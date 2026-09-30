"""Adapter for GWASPoker (GWASPokerforPRS2) - optional fast discovery backend.

GWASPoker is used for what it is good at: EFO-driven study search, ancestry
metadata, and cheap file-availability / GWAS-SSF checks. Its results are
converted into GWAS2Mechanism models and then pass through the *same*
phenotype matching, ancestry classification and ranking as every other backend,
so no GWASPoker heuristics leak into selection. Any failure falls back to the
official GWAS Catalog backend (see ``discovery.service``).

Pinned: MuhammadMuneeb007/GWASPokerforPRS2@11538899ad0d8b3962a6d172d644fae43295a34b
"""

from __future__ import annotations

import logging
from typing import Any

from gwas2mechanism.discovery.ancestry import parse_ancestry_string
from gwas2mechanism.discovery.base import DiscoveryBackend
from gwas2mechanism.discovery.models import (
    CandidateFile,
    CandidateStudy,
    OntologyTerm,
    PhenotypeResolution,
)
from gwas2mechanism.discovery.phenotype import canonical_term_id

log = logging.getLogger(__name__)

#: Only non-model sample-size sources are accepted (no LLM-extracted counts).
ACCEPTED_SAMPLE_SOURCES = {"structured", "api", "rest", "text", "regex", "ssf"}


def gwaspoker_available() -> bool:
    try:
        import gwaspoker.catalog.discovery  # noqa: F401
    except Exception:
        return False
    return True


class GWASPokerBackend(DiscoveryBackend):
    name = "gwaspoker"

    def __init__(self, max_studies: int = 500, check_files: bool = True):
        from gwaspoker.catalog.discovery import DiscoveryService

        self._service = DiscoveryService(enable_llm=False)
        self.max_studies = max_studies
        self.check_files = check_files
        self._results: dict[str, Any] = {}

    def close(self) -> None:
        self._service.close()

    def ontology_terms(self, phenotype: str) -> list[OntologyTerm]:
        terms = self._service.catalog.search_traits(phenotype, limit=25)
        return [
            OntologyTerm(canonical_term_id(t.efo_id), t.label, iri=t.uri, source="gwaspoker")
            for t in terms
            if getattr(t, "efo_id", None)
        ]

    def search_studies(self, resolution: PhenotypeResolution) -> list[CandidateStudy]:
        results = self._service.search(
            resolution.phenotype_canonical,
            population=None,
            limit=self.max_studies,
            summary_stats_only=False,
            resolve_samples=True,
        )
        if self.check_files and results:
            results = self._service.check_files(results)
        studies = []
        for result in results:
            self._results[result.study.study_accession] = result
            studies.append(self._convert(result))
        return studies

    def _convert(self, result: Any) -> CandidateStudy:
        study = result.study
        ancestries = []
        for block in study.ancestries:
            groups = list(block.ancestral_groups) or ["NR"]
            records = [
                record
                for group in groups
                for record in parse_ancestry_string(group, block.stage)
            ]
            if block.number_of_individuals and len(records) == 1 and records[0].n is None:
                records[0].n = block.number_of_individuals
            if block.countries_of_recruitment:
                for r in records:
                    r.country = ", ".join(block.countries_of_recruitment)
            ancestries += records
        counts = study.samples
        candidate = CandidateStudy(
            study_accession=study.study_accession,
            reported_trait=study.reported_trait,
            mapped_trait_ids=[canonical_term_id(t.efo_id) for t in study.mapped_traits if t.efo_id],
            mapped_trait_labels=[t.label for t in study.mapped_traits],
            background_trait_labels=[t.label for t in study.background_traits],
            pubmed_id=study.pubmed_id,
            first_author=study.first_author,
            publication_date=study.publication_date,
            cohorts=list(study.cohorts),
            initial_sample_size_text=study.initial_sample_description,
            ancestries=ancestries,
            snp_count=study.snp_count,
            summary_stats_available=study.summary_statistics_available,
            summary_stats_url=study.summary_statistics_location,
            backend="gwaspoker",
        )
        if counts.cases is not None and str(getattr(counts.cases_source, "value", "")).lower() in ACCEPTED_SAMPLE_SOURCES:
            candidate.n_cases = counts.cases
        if counts.controls is not None and str(getattr(counts.controls_source, "value", "")).lower() in ACCEPTED_SAMPLE_SOURCES:
            candidate.n_controls = counts.controls
        return candidate

    def list_files(self, study: CandidateStudy) -> list[CandidateFile]:
        files = super().list_files(study)
        result = self._results.get(study.study_accession)
        if result is not None and result.ssf_status:
            for f in files:
                if f.kind in ("raw", "harmonised") and f.file_type is None:
                    f.file_type = result.ssf_status
        return files
