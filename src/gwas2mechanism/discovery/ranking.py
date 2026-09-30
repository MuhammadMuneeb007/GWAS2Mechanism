"""Explicit, interpretable study ranking and per-population selection.

``selection_score = sum_k weight_k * component_k - sum_j penalty_j`` with every
component, weight and penalty exported in ``study_selection.parquet``.
Unknown values score 0 and are flagged - they are never imputed. Sample sizes
are never fabricated.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

import polars as pl

from gwas2mechanism.config import Config
from gwas2mechanism.constants import COMBINED_PUBLISHED, MULTI, SUPERPOPULATIONS
from gwas2mechanism.discovery.models import CandidateFile, CandidateStudy
from gwas2mechanism.discovery.phenotype import is_stratified_subset


def effective_n(n_cases: float | None, n_controls: float | None) -> float | None:
    """N_eff = 4 / (1/N_cases + 1/N_controls); None unless both are positive."""
    if not n_cases or not n_controls or n_cases <= 0 or n_controls <= 0:
        return None
    return 4.0 / (1.0 / n_cases + 1.0 / n_controls)


def _file_summary(files: list[CandidateFile]) -> dict[str, object]:
    data = [f for f in files if f.kind in ("harmonised", "raw")]
    harmonised = [f for f in data if f.is_harmonised]
    meta = [f for f in files if f.kind == "meta_yaml"]
    best = (harmonised or data or [None])[0]
    grch38 = bool(harmonised) or any(f.build == "GRCh38" for f in data + meta)
    probed = [f for f in (harmonised + data) if f.probe_status == "header_parsed"]
    ref = probed[0] if probed else None
    return {
        "n_files": len(data),
        "harmonised_available": bool(harmonised),
        "grch38_available": grch38,
        "build": "GRCh38" if grch38 else (next((f.build for f in data + meta if f.build), None)),
        "selected_file": best.url if best else None,
        "file_type": next((f.file_type for f in meta + data if f.file_type), None),
        "columns_verified": ref is not None,
        "has_beta": ref.has_beta if ref else None,
        "has_or": ref.has_or if ref else None,
        "has_se": ref.has_se if ref else None,
        "has_eaf": ref.has_eaf if ref else None,
        "has_n_column": ref.has_n if ref else None,
        "has_alleles": ref.has_alleles if ref else None,
        "probe_bytes": sum(f.probe_bytes for f in files),
    }


def build_candidate_table(studies: Iterable[CandidateStudy], files_by_study: dict[str, list[CandidateFile]]) -> pl.DataFrame:
    rows = []
    for s in studies:
        row = {
            "study_accession": s.study_accession,
            "reported_trait": s.reported_trait,
            "mapped_trait": "; ".join(s.mapped_trait_labels) or None,
            "mapped_trait_ids": "; ".join(s.mapped_trait_ids) or None,
            "n_mapped_traits": len(s.mapped_trait_ids),
            "background_traits": "; ".join(s.background_trait_labels) or None,
            "population": s.population,
            "analysis_label": s.analysis_label,
            "reported_ancestry": s.reported_ancestry,
            "population_approximate": s.population_approximate,
            "population_reason": s.population_reason,
            "n_total": s.n_total,
            "n_cases": s.n_cases,
            "n_controls": s.n_controls,
            "n_eff": effective_n(s.n_cases, s.n_controls),
            "n_source": s.n_source,
            "trait_type": s.trait_type,
            "initial_sample_size": s.initial_sample_size_text,
            "pubmed_id": s.pubmed_id,
            "first_author": s.first_author,
            "publication_date": s.publication_date,
            "cohorts": "; ".join(s.cohorts) or None,
            "gxe": s.gxe,
            "summary_stats_available": bool(s.summary_stats_available),
            "summary_stats_url": s.summary_stats_url,
            "phenotype_match_type": s.match_type,
            "phenotype_match_score": s.match_score,
            "phenotype_match_reason": s.match_reason,
            "stratified_subset": is_stratified_subset(s.reported_trait),
            "backend": s.backend,
        }
        row.update(_file_summary(files_by_study.get(s.study_accession, [])))
        rows.append(row)
    schema_overrides = {
        "n_total": pl.Float64,
        "n_cases": pl.Float64,
        "n_controls": pl.Float64,
        "n_eff": pl.Float64,
        "has_beta": pl.Boolean,
        "has_or": pl.Boolean,
        "has_se": pl.Boolean,
        "has_eaf": pl.Boolean,
        "has_n_column": pl.Boolean,
        "has_alleles": pl.Boolean,
    }
    if not rows:
        return pl.DataFrame(schema={"study_accession": pl.String, **schema_overrides})
    return pl.DataFrame(rows, schema_overrides=schema_overrides, infer_schema_length=None)


def score_candidates(table: pl.DataFrame, config: Config) -> pl.DataFrame:
    if table.is_empty():
        return table
    w = config.ranking.weights
    pen = config.ranking.penalties
    disc = config.discovery
    n_for_rank = pl.coalesce(pl.col("n_eff"), pl.col("n_total"))
    scored = table.with_columns(
        comp_phenotype_match=pl.col("phenotype_match_score"),
        comp_ancestry_suitability=pl.when(pl.col("population").is_in(list(SUPERPOPULATIONS)))
        .then(pl.when(pl.col("population_approximate")).then(0.7).otherwise(1.0))
        .when(pl.col("population") == MULTI)
        .then(0.5)
        .otherwise(0.0),
        comp_summary_stats_available=pl.col("summary_stats_available").cast(pl.Float64),
        comp_harmonised_available=pl.col("harmonised_available").fill_null(False).cast(pl.Float64),
        comp_grch38_available=pl.col("grch38_available").fill_null(False).cast(pl.Float64),
        comp_sample_size=(n_for_rank.log10() / 7.0).clip(0.0, 1.0).fill_null(0.0),
        comp_effect_columns=(
            (pl.col("has_beta").fill_null(False) | pl.col("has_or").fill_null(False))
            & pl.col("has_se").fill_null(False)
        ).cast(pl.Float64),
        comp_eaf_available=pl.col("has_eaf").fill_null(False).cast(pl.Float64),
        comp_n_column_available=pl.col("has_n_column").fill_null(False).cast(pl.Float64),
        comp_metadata_completeness=pl.mean_horizontal(
            pl.col("pubmed_id").is_not_null().cast(pl.Float64),
            pl.col("n_total").is_not_null().cast(pl.Float64),
            pl.col("reported_ancestry").is_not_null().cast(pl.Float64),
            pl.col("cohorts").is_not_null().cast(pl.Float64),
            pl.col("mapped_trait").is_not_null().cast(pl.Float64),
        ),
        penalty_stratified_subset=pl.col("stratified_subset").cast(pl.Float64) * pen.get("stratified_subset", 0.0),
        penalty_multi_trait_mapping=(pl.col("n_mapped_traits") > 1).cast(pl.Float64) * pen.get("multi_trait_mapping", 0.0),
    )
    components = [c for c in scored.columns if c.startswith("comp_")]
    score = pl.sum_horizontal(
        [pl.col(c) * w.get(c.removeprefix("comp_"), 0.0) for c in components]
    ) - pl.col("penalty_stratified_subset") - pl.col("penalty_multi_trait_mapping")
    max_score = sum(w.get(c.removeprefix("comp_"), 0.0) for c in components)
    data_quality = pl.mean_horizontal(
        [pl.col(c) for c in ("comp_harmonised_available", "comp_grch38_available", "comp_effect_columns", "comp_eaf_available", "comp_n_column_available")]
    )

    min_match = config.ranking.min_phenotype_match
    exclusion = (
        pl.when(~pl.col("summary_stats_available")).then(pl.lit("no full summary statistics"))
        .when(pl.col("phenotype_match_score") < (0.0 if config.phenotype_resolution.allow_fuzzy else min_match))
        .then(pl.lit("phenotype match below threshold"))
        .when(pl.col("phenotype_match_score") <= 0).then(pl.lit("phenotype does not match"))
        .when(pl.col("analysis_label").is_null()).then(pl.lit("population UNKNOWN/UNMAPPED - no matching LD reference"))
        .when(pl.lit(disc.exclude_gxe) & pl.col("gxe")).then(pl.lit("GxE/GxG interaction study"))
        .when(pl.lit(disc.exclude_background_traits) & pl.col("background_traits").is_not_null())
        .then(pl.lit("conditional/background-trait analysis"))
        .when(pl.col("has_alleles") == False)  # noqa: E712 - explicit tri-state check
        .then(pl.lit("summary statistics lack allele columns"))
        .when(pl.lit(config.harmonization.require_grch38 and config.harmonization.liftover == "none") & ~pl.col("grch38_available"))
        .then(pl.lit("not GRCh38 and liftover disabled"))
        .otherwise(pl.lit(None, dtype=pl.String))
    )
    return scored.with_columns(
        selection_score=score,
        selection_score_max=pl.lit(max_score),
        data_quality_score=data_quality,
        exclusion_reason=exclusion,
    ).with_columns(usable=pl.col("exclusion_reason").is_null())


def select_studies(
    scored: pl.DataFrame,
    config: Config,
    populations: list[str] | None = None,
) -> pl.DataFrame:
    """Pick the strongest usable study per analysis label (or all, with all_studies)."""
    if scored.is_empty():
        return scored.with_columns(selected=pl.lit(False), selection_reason=pl.lit(None, pl.String), analysis_id=pl.lit(None, pl.String), is_primary=pl.lit(False))
    wanted = list(populations or config.populations.allowed)
    labels = [p for p in wanted if p in SUPERPOPULATIONS]
    if config.combined.published_combined:
        labels.append(COMBINED_PUBLISHED)
    ranked = scored.sort(
        ["analysis_label", "usable", "selection_score", "n_eff", "n_total"],
        descending=[False, True, True, True, True],
        nulls_last=True,
    ).with_columns(
        rank_in_population=pl.col("selection_score").rank("ordinal", descending=True).over(["analysis_label", "usable"])
    )
    eligible = pl.col("usable") & pl.col("analysis_label").is_in(labels)
    if config.discovery.all_studies:
        selected = eligible
    else:
        selected = eligible & (pl.col("rank_in_population") == 1)
    out = ranked.with_columns(selected=selected).with_columns(
        is_primary=pl.col("selected") & (pl.col("rank_in_population") == 1),
    )
    out = out.with_columns(
        analysis_id=pl.when(pl.col("selected") & pl.col("is_primary"))
        .then(pl.col("analysis_label"))
        .when(pl.col("selected"))
        .then(pl.concat_str([pl.col("analysis_label"), pl.col("study_accession")], separator="__"))
        .otherwise(pl.lit(None, pl.String)),
        selection_reason=pl.when(pl.col("selected") & pl.col("is_primary"))
        .then(
            pl.format(
                "highest selection score ({}/{}) among usable {} studies; match={}",
                pl.col("selection_score").round(3),
                pl.col("selection_score_max").round(1),
                pl.col("analysis_label"),
                pl.col("phenotype_match_type"),
            )
        )
        .when(pl.col("selected"))
        .then(pl.lit("additional usable study (--all-studies)"))
        .when(~pl.col("usable"))
        .then(pl.concat_str([pl.lit("excluded: "), pl.col("exclusion_reason")]))
        .when(~pl.col("analysis_label").is_in(labels))
        .then(pl.lit("population not requested"))
        .otherwise(pl.format("usable but ranked #{} for {}", pl.col("rank_in_population"), pl.col("analysis_label"))),
    )
    return out


SELECTION_EXPORT_COLUMNS = [
    "study_accession",
    "analysis_id",
    "selected",
    "is_primary",
    "selection_reason",
    "selection_score",
    "data_quality_score",
    "phenotype_match_score",
    "phenotype_match_type",
    "phenotype_match_reason",
    "reported_trait",
    "mapped_trait",
    "population",
    "analysis_label",
    "reported_ancestry",
    "population_approximate",
    "n_total",
    "n_cases",
    "n_controls",
    "n_eff",
    "n_source",
    "summary_stats_available",
    "harmonised_available",
    "build",
    "has_beta",
    "has_or",
    "has_se",
    "has_eaf",
    "has_n_column",
    "columns_verified",
    "exclusion_reason",
    "pubmed_id",
    "selected_file",
]


def score_value(n: float | None) -> float:
    """Helper used in docs/tests: the sample-size component for a given N."""
    if not n or n <= 0:
        return 0.0
    return max(0.0, min(1.0, math.log10(n) / 7.0))
