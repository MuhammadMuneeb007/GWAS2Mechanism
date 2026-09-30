"""Official GWAS Catalog backend - the canonical provenance source.

Two modes:

* ``rest`` (default): REST API v2 (``/efo-traits``, ``/studies?efo_id=``,
  ``/studies?disease_trait=``). Small metadata requests only.
* ``bulk``: the official download TSVs (studies + ancestries), cached once
  and scanned lazily with Polars. Also used for offline/air-gapped runs and
  the synthetic smoke test (paths may be local).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import httpx
import polars as pl

from gwas2mechanism.config import DiscoveryConfig
from gwas2mechanism.discovery.ancestry import parse_ancestry_string
from gwas2mechanism.discovery.base import DiscoveryBackend, gcst_directory
from gwas2mechanism.discovery.models import (
    AncestryRecord,
    CandidateStudy,
    OntologyTerm,
    PhenotypeResolution,
)
from gwas2mechanism.discovery.phenotype import canonical_term_id, normalize_text
from gwas2mechanism.utils.download import USER_AGENT, Downloader, is_local, local_path

log = logging.getLogger(__name__)


class GWASCatalogBackend(DiscoveryBackend):
    name = "gwas_catalog"

    def __init__(self, config: DiscoveryConfig, cache_dir: Path):
        self.config = config
        self.cfg = config.gwas_catalog
        self.cache_dir = cache_dir
        self._client: httpx.Client | None = None

    # ------------------------------------------------------------------ http
    @property
    def client(self) -> httpx.Client:
        if self._client is None:
            self._client = httpx.Client(
                timeout=self.config.timeout_s,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
                follow_redirects=True,
            )
        return self._client

    def close(self) -> None:
        if self._client is not None:
            self._client.close()

    def _get(self, path: str, params: dict[str, Any]) -> dict[str, Any]:
        response = self.client.get(f"{self.cfg.rest_url}{path}", params=params)
        if response.status_code == 404:
            return {}
        response.raise_for_status()
        return response.json()

    def _paged(self, path: str, params: dict[str, Any], key: str) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        page = 0
        while len(out) < self.cfg.max_studies:
            payload = self._get(path, {**params, "size": self.cfg.page_size, "page": page})
            items = payload.get("_embedded", {}).get(key, [])
            out.extend(items)
            info = payload.get("page", {})
            if not items or page + 1 >= int(info.get("totalPages", 1) or 1):
                break
            page += 1
        return out

    # ------------------------------------------------------------ ontology
    def ontology_terms(self, phenotype: str) -> list[OntologyTerm]:
        if self.cfg.mode == "bulk":
            return self._bulk_terms(phenotype)
        try:
            items = self._paged("/efo-traits", {"trait": phenotype}, "efo_traits")
        except httpx.HTTPError as exc:
            log.warning("GWAS Catalog trait search failed: %s", exc)
            return []
        return [
            OntologyTerm(
                term_id=canonical_term_id(x["efo_id"]),
                label=x.get("efo_trait") or x["efo_id"],
                iri=x.get("uri"),
                source="gwas_catalog",
            )
            for x in items
            if x.get("efo_id")
        ]

    # -------------------------------------------------------------- studies
    def search_studies(self, resolution: PhenotypeResolution) -> list[CandidateStudy]:
        if self.cfg.mode == "bulk":
            return self._bulk_studies(resolution)
        found: dict[str, dict[str, Any]] = {}
        term_ids = sorted(resolution.term_ids | set(resolution.descendant_ids))
        for term_id in term_ids:
            for item in self._paged("/studies", {"efo_id": term_id}, "studies"):
                found.setdefault(item["accession_id"], item)
        # Reported-trait query catches studies whose ontology mapping differs.
        texts = dict.fromkeys([resolution.phenotype_input.strip(), resolution.phenotype_canonical.strip()])
        for text in texts:
            for item in self._paged("/studies", {"disease_trait": text}, "studies"):
                found.setdefault(item["accession_id"], item)
        return [self._study_from_rest(item) for item in found.values()]

    def _study_from_rest(self, item: dict[str, Any]) -> CandidateStudy:
        ancestries: list[AncestryRecord] = []
        for text in item.get("discovery_ancestry") or []:
            ancestries += parse_ancestry_string(text, "initial")
        for text in item.get("replication_ancestry") or []:
            ancestries += parse_ancestry_string(text, "replication")
        url = item.get("full_summary_stats")
        if not url and item.get("full_summary_stats_available"):
            url = gcst_directory(self.cfg.sumstats_base, item["accession_id"])
        return CandidateStudy(
            study_accession=item["accession_id"],
            reported_trait=item.get("disease_trait"),
            mapped_trait_ids=[canonical_term_id(t["efo_id"]) for t in item.get("efo_traits") or []],
            mapped_trait_labels=[t.get("efo_trait", "") for t in item.get("efo_traits") or []],
            background_trait_labels=[t.get("efo_trait", "") for t in item.get("bg_efo_traits") or []],
            pubmed_id=str(item["pubmed_id"]) if item.get("pubmed_id") else None,
            cohorts=list(item.get("cohort") or []),
            initial_sample_size_text=item.get("initial_sample_size"),
            ancestries=ancestries,
            gxe=bool(item.get("gxe") or item.get("gxg")),
            snp_count=item.get("snp_count"),
            summary_stats_available=item.get("full_summary_stats_available"),
            summary_stats_url=url,
            backend="gwas_catalog_rest",
        )

    # ------------------------------------------------------------------ bulk
    def _bulk_file(self, source: str, name: str) -> Path:
        if is_local(source):
            return local_path(source)
        dest = self.cache_dir / "gwas_catalog" / name
        Downloader().fetch(source, dest)
        return dest

    def _bulk_frames(self) -> tuple[pl.LazyFrame, pl.LazyFrame]:
        studies = self._bulk_file(self.cfg.studies_source, "studies.tsv")
        ancestries = self._bulk_file(self.cfg.ancestries_source, "ancestries.tsv")
        opts = {"separator": "\t", "quote_char": None, "infer_schema": False, "truncate_ragged_lines": True}
        return pl.scan_csv(studies, **opts), pl.scan_csv(ancestries, **opts)

    def _bulk_terms(self, phenotype: str) -> list[OntologyTerm]:
        studies, _ = self._bulk_frames()
        query = normalize_text(phenotype)
        pairs = (
            studies.select(
                pl.col("MAPPED_TRAIT").str.split(", ").alias("label"),
                pl.col("MAPPED_TRAIT_URI").str.split(", ").alias("uri"),
            )
            .filter(pl.col("label").list.len() == pl.col("uri").list.len())
            .explode(["label", "uri"])
            .filter(pl.col("label").str.to_lowercase().str.contains(query, literal=True))
            .unique()
            .collect()
        )
        return [
            OntologyTerm(canonical_term_id(r["uri"]), r["label"], iri=r["uri"], source="bulk")
            for r in pairs.iter_rows(named=True)
        ]

    def _bulk_studies(self, resolution: PhenotypeResolution) -> list[CandidateStudy]:
        studies, ancestries = self._bulk_frames()
        ids = sorted(resolution.term_ids | set(resolution.descendant_ids))
        syn = list(resolution.synonyms)
        # Cheap vectorised PRE-filter only; match_study() makes the decision.
        prefilter = pl.lit(False)
        for term_id in ids:
            prefilter = prefilter | pl.col("MAPPED_TRAIT_URI").str.contains(term_id, literal=True)
        for text in syn:
            prefilter = prefilter | pl.col("DISEASE/TRAIT").str.to_lowercase().str.contains(text, literal=True)
            prefilter = prefilter | pl.col("MAPPED_TRAIT").str.to_lowercase().str.contains(text, literal=True)
        hits = studies.filter(prefilter).collect()
        if hits.is_empty():
            return []
        acc = hits["STUDY ACCESSION"].to_list()
        anc = ancestries.filter(pl.col("STUDY ACCESSION").is_in(acc)).collect()
        anc_by: dict[str, list[AncestryRecord]] = {}
        for row in anc.iter_rows(named=True):
            stage = (row.get("STAGE") or "unknown").strip().lower()
            n_text = (row.get("NUMBER OF INDIVIDUALS") or "").replace(",", "").strip()
            categories = [c.strip() for c in (row.get("BROAD ANCESTRAL CATEGORY") or "NR").split(",")]
            records = [
                AncestryRecord(stage=stage, category=c, country=row.get("COUNTRY OF RECRUITMENT")) for c in categories if c
            ]
            if records and n_text.isdigit():
                records[0].n = int(n_text)
            anc_by.setdefault(row["STUDY ACCESSION"], []).extend(records)

        out = []
        for row in hits.iter_rows(named=True):
            accession = row["STUDY ACCESSION"]
            uris = [u for u in (row.get("MAPPED_TRAIT_URI") or "").split(", ") if u]
            labels = [x for x in (row.get("MAPPED_TRAIT") or "").split(", ") if x]
            has_ss = (row.get("FULL SUMMARY STATISTICS") or "").strip().lower() == "yes"
            location = (row.get("SUMMARY STATS LOCATION") or "").strip() or None
            if has_ss and not location:
                location = gcst_directory(self.cfg.sumstats_base, accession)
            if location and is_local(self.cfg.sumstats_base) and not is_local(location):
                location = gcst_directory(self.cfg.sumstats_base, accession)
            out.append(
                CandidateStudy(
                    study_accession=accession,
                    reported_trait=row.get("DISEASE/TRAIT"),
                    mapped_trait_ids=[canonical_term_id(u) for u in uris],
                    mapped_trait_labels=labels,
                    background_trait_labels=[
                        x for x in (row.get("MAPPED BACKGROUND TRAIT") or "").split(", ") if x
                    ],
                    pubmed_id=row.get("PUBMED ID"),
                    first_author=row.get("FIRST AUTHOR"),
                    publication_date=row.get("DATE"),
                    cohorts=[c for c in (row.get("COHORT") or "").split("|") if c],
                    initial_sample_size_text=row.get("INITIAL SAMPLE SIZE"),
                    ancestries=anc_by.get(accession, []),
                    gxe=(row.get("GXE") or "").strip().lower() == "yes",
                    summary_stats_available=has_ss,
                    summary_stats_url=location,
                    backend="gwas_catalog_bulk",
                )
            )
        return out
