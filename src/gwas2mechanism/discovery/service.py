"""Discovery orchestration and outputs.

Outputs (``runs/<phenotype>/<run>/discovery/``):
    phenotype_resolution.json
    candidate_studies.parquet
    candidate_files.parquet
    study_selection.parquet
    selection_report.tsv / selection_report.html
    selected_analyses.json
"""

from __future__ import annotations

import html
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import polars as pl

from gwas2mechanism.config import Config
from gwas2mechanism.discovery.ancestry import apply_sample_size, classify_study
from gwas2mechanism.discovery.base import DiscoveryBackend
from gwas2mechanism.discovery.gwas_catalog import GWASCatalogBackend
from gwas2mechanism.discovery.gwaspoker import GWASPokerBackend, gwaspoker_available
from gwas2mechanism.discovery.models import CandidateFile, CandidateStudy, PhenotypeResolution
from gwas2mechanism.discovery.phenotype import (
    PhenotypeResolver,
    match_study,
    ols_descendants_factory,
    ols_search_factory,
)
from gwas2mechanism.discovery.ranking import (
    SELECTION_EXPORT_COLUMNS,
    build_candidate_table,
    score_candidates,
    select_studies,
)
from gwas2mechanism.paths import RunLayout
from gwas2mechanism.utils.io import (
    atomic_write_json,
    atomic_write_text,
    ensure_dir,
    write_parquet,
    write_tsv,
)

log = logging.getLogger(__name__)


def make_backend(config: Config) -> list[DiscoveryBackend]:
    """Ordered backends: GWASPoker first when requested/available, Catalog always last."""
    cache = config.resources_path
    catalog = GWASCatalogBackend(config.discovery, cache)
    choice = config.discovery.backend
    if choice == "gwas_catalog" or config.discovery.gwas_catalog.mode == "bulk":
        return [catalog]
    if choice in ("auto", "gwaspoker") and gwaspoker_available():
        try:
            # GWASPoker is the fast enumeration pass. File checks happen only
            # after preliminary ranking below; asking it to inspect thousands
            # of candidates here defeats bounded discovery.
            return [
                GWASPokerBackend(
                    max_studies=min(config.discovery.gwas_catalog.max_studies, 200),
                    check_files=False,
                ),
                catalog,
            ]
        except Exception as exc:  # pragma: no cover - depends on optional package
            log.warning("GWASPoker unavailable (%s); using GWAS Catalog", exc)
    elif choice == "gwaspoker":
        log.warning("GWASPoker requested but not installed; using GWAS Catalog backend")
    return [catalog]


def resolve_phenotype(phenotype: str, config: Config, backend: DiscoveryBackend, offline: bool = False) -> PhenotypeResolution:
    cfg = config.phenotype_resolution
    search = None if offline else ols_search_factory(cfg.ols_url, config.discovery.timeout_s)
    desc = None if offline else ols_descendants_factory(cfg.ols_url, config.discovery.timeout_s)
    try:
        extra = backend.ontology_terms(phenotype)
    except Exception as exc:
        log.warning("Backend trait lookup failed: %s", exc)
        extra = []
    return PhenotypeResolver(cfg, search, desc, extra_terms=extra).resolve(phenotype)


def _search_with_fallback(backends: list[DiscoveryBackend], resolution: PhenotypeResolution) -> tuple[list[CandidateStudy], DiscoveryBackend]:
    last_error: Exception | None = None
    for backend in backends:
        try:
            studies = backend.search_studies(resolution)
            log.info("%s returned %d candidate studies", backend.name, len(studies))
            if studies or backend is backends[-1]:
                return studies, backend
        except Exception as exc:
            last_error = exc
            log.warning("Discovery backend %s failed: %s", backend.name, exc)
    if last_error:
        raise last_error
    return [], backends[-1]


def discover(
    phenotype: str,
    config: Config,
    layout: RunLayout,
    populations: list[str] | None = None,
    offline: bool = False,
) -> pl.DataFrame:
    ensure_dir(layout.discovery)
    backends = make_backend(config)
    try:
        resolution = resolve_phenotype(phenotype, config, backends[0], offline=offline)
        atomic_write_json(layout.phenotype_resolution, resolution.to_dict())
        log.info(
            "Phenotype %r -> %r (%s; terms=%s)",
            phenotype,
            resolution.phenotype_canonical,
            resolution.resolution_type,
            sorted(resolution.term_ids),
        )
        studies, backend = _search_with_fallback(backends, resolution)

        for study in studies:
            classify_study(study, config.populations.mapping)
            apply_sample_size(study)
            match_study(study, resolution, config.phenotype_resolution.fuzzy_min_ratio)
        studies = [s for s in studies if s.match_score > 0 or config.phenotype_resolution.allow_fuzzy]

        # Preliminary (metadata-only) ranking decides which studies to list/probe.
        prelim = score_candidates(build_candidate_table(studies, {}), config)
        files_by_study: dict[str, list[CandidateFile]] = {}
        if config.discovery.list_files and not prelim.is_empty():
            top = (
                prelim.filter(pl.col("usable"))
                .sort("selection_score", descending=True)
                .group_by("analysis_label", maintain_order=True)
                .head(max(config.discovery.probe_top_n * 3, 1))
            )
            chosen = set(top["study_accession"].to_list())
            by_acc = {s.study_accession: s for s in studies}
            targets = [by_acc[a] for a in chosen]
            with ThreadPoolExecutor(max_workers=config.discovery.concurrency) as pool:
                for study, files in zip(
                    targets, pool.map(backend.list_files, targets), strict=True
                ):
                    files_by_study[study.study_accession] = files
            if config.discovery.probe_headers:
                probe_ids = set(
                    top.group_by("analysis_label", maintain_order=True)
                    .head(config.discovery.probe_top_n)["study_accession"]
                    .to_list()
                )
                to_probe = [
                    f
                    for acc in probe_ids
                    for f in _probe_targets(files_by_study.get(acc, []), config.discovery.prefer_harmonised)
                ]
                with ThreadPoolExecutor(max_workers=config.discovery.concurrency) as pool:
                    list(pool.map(lambda f: backend.probe(f, config.discovery.probe_bytes), to_probe))

        table = build_candidate_table(studies, files_by_study)
        scored = score_candidates(table, config)
        selection = select_studies(scored, config, populations)
    finally:
        for b in backends:
            b.close()

    write_parquet(_studies_frame(studies), layout.candidate_studies)
    write_parquet(_files_frame(files_by_study), layout.candidate_files)
    write_parquet(selection, layout.study_selection)
    export = selection.select([c for c in SELECTION_EXPORT_COLUMNS if c in selection.columns]) if not selection.is_empty() else selection
    write_tsv(export, layout.selection_report_tsv)
    atomic_write_text(layout.selection_report_html, selection_html(resolution, export))
    atomic_write_json(layout.selected_analyses, selected_analyses(selection))
    return selection


def _probe_targets(files: list[CandidateFile], prefer_harmonised: bool) -> list[CandidateFile]:
    data = [f for f in files if f.kind in ("harmonised", "raw")]
    if prefer_harmonised:
        data.sort(key=lambda f: not f.is_harmonised)
    return [f for f in files if f.kind == "meta_yaml"][:1] + data[:1]


def _studies_frame(studies: list[CandidateStudy]) -> pl.DataFrame:
    rows = [s.to_row() for s in studies]
    for row in rows:
        row["ancestries"] = str(row["ancestries"])
        for key in ("mapped_trait_ids", "mapped_trait_labels", "background_trait_labels", "cohorts"):
            row[key] = "; ".join(row[key])
    if not rows:
        return pl.DataFrame({"study_accession": []}, schema={"study_accession": pl.String})
    return pl.DataFrame(rows, infer_schema_length=None)


def _files_frame(files_by_study: dict[str, list[CandidateFile]]) -> pl.DataFrame:
    rows = [f.to_row() for files in files_by_study.values() for f in files]
    for row in rows:
        row["header_columns"] = "\t".join(row["header_columns"])
    if not rows:
        return pl.DataFrame({"study_accession": [], "url": []}, schema={"study_accession": pl.String, "url": pl.String})
    return pl.DataFrame(rows, infer_schema_length=None)


def selected_analyses(selection: pl.DataFrame) -> list[dict[str, Any]]:
    if selection.is_empty():
        return []
    keep = [
        "analysis_id",
        "analysis_label",
        "population",
        "study_accession",
        "is_primary",
        "reported_trait",
        "reported_ancestry",
        "population_approximate",
        "n_total",
        "n_cases",
        "n_controls",
        "n_eff",
        "trait_type",
        "selected_file",
        "summary_stats_url",
        "harmonised_available",
        "build",
    ]
    return selection.filter(pl.col("selected")).select([c for c in keep if c in selection.columns]).to_dicts()


def selection_html(resolution: PhenotypeResolution, table: pl.DataFrame) -> str:
    head = "".join(f"<th>{html.escape(c)}</th>" for c in table.columns)
    body = []
    for row in table.iter_rows():
        cls = ' class="sel"' if table.columns and "selected" in table.columns and row[table.columns.index("selected")] else ""
        body.append("<tr" + cls + ">" + "".join(f"<td>{html.escape('' if v is None else str(v))}</td>" for v in row) + "</tr>")
    terms = ", ".join(f"{html.escape(t.label)} ({t.term_id})" for t in resolution.terms) or "none (free text)"
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>Study selection</title>
<style>body{{font-family:system-ui,sans-serif;margin:16px}}table{{border-collapse:collapse;font-size:12px}}
td,th{{border:1px solid #ccc;padding:3px 6px}}tr.sel{{background:#e8f5e9}}th{{background:#f3f3f3;position:sticky;top:0}}</style>
</head><body><h1>GWAS study selection</h1>
<p><b>Phenotype input:</b> {html.escape(resolution.phenotype_input)} &nbsp; <b>Canonical:</b> {html.escape(resolution.phenotype_canonical)}
&nbsp; <b>Resolution:</b> {html.escape(resolution.resolution_type)} &nbsp; <b>Terms:</b> {terms}</p>
<p>Every candidate is listed with the components of its selection score and the reason it was selected or excluded.
Selected rows are highlighted. Sample sizes are taken from Catalog metadata and are never imputed.</p>
<div style="overflow:auto"><table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>
</body></html>"""
