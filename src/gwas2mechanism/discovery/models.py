"""Data models shared by discovery backends, the resolver and the ranker."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class OntologyTerm:
    term_id: str  # e.g. MONDO_0005277 (underscore form, as used by the GWAS Catalog)
    label: str
    synonyms: tuple[str, ...] = ()
    iri: str | None = None
    source: str = "unknown"  # ols | gwas_catalog | bulk | user

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class PhenotypeResolution:
    phenotype_input: str
    phenotype_canonical: str
    phenotype_slug: str
    terms: list[OntologyTerm]
    synonyms: list[str]
    descendant_ids: list[str] = field(default_factory=list)
    rejected_terms: list[dict[str, Any]] = field(default_factory=list)
    resolution_type: str = "free_text"  # exact_label | exact_synonym | user_pinned | fuzzy | free_text
    source: str = "none"

    @property
    def term_ids(self) -> set[str]:
        return {t.term_id for t in self.terms}

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["term_ids"] = sorted(self.term_ids)
        return out


@dataclass
class AncestryRecord:
    stage: str  # initial | replication | unknown
    category: str  # original GWAS Catalog broad ancestral category
    n: int | None = None
    country: str | None = None
    population: str = "UNKNOWN"  # EUR/AFR/EAS/SAS/AMR/UNMAPPED/UNKNOWN
    approximate: bool = False


@dataclass
class CandidateStudy:
    study_accession: str
    reported_trait: str | None = None
    mapped_trait_ids: list[str] = field(default_factory=list)
    mapped_trait_labels: list[str] = field(default_factory=list)
    background_trait_labels: list[str] = field(default_factory=list)
    pubmed_id: str | None = None
    first_author: str | None = None
    publication_date: str | None = None
    cohorts: list[str] = field(default_factory=list)
    initial_sample_size_text: str | None = None
    ancestries: list[AncestryRecord] = field(default_factory=list)
    gxe: bool = False
    snp_count: int | None = None
    summary_stats_available: bool | None = None
    summary_stats_url: str | None = None
    backend: str = "gwas_catalog"
    # filled by classification / resolution / ranking
    population: str = "UNKNOWN"
    reported_ancestry: str | None = None
    population_approximate: bool = False
    population_reason: str = ""
    analysis_label: str | None = None
    n_total: int | None = None
    n_cases: int | None = None
    n_controls: int | None = None
    n_source: str | None = None
    trait_type: str | None = None  # binary | quantitative | None
    match_type: str = "none"
    match_score: float = 0.0
    match_reason: str = ""

    def to_row(self) -> dict[str, Any]:
        row = asdict(self)
        row["ancestries"] = [asdict(a) for a in self.ancestries]
        return row


@dataclass
class CandidateFile:
    study_accession: str
    url: str
    filename: str
    kind: str  # harmonised | raw | meta_yaml | md5 | other
    size: int | None = None
    build: str | None = None
    is_harmonised: bool = False
    file_type: str | None = None  # GWAS-SSF / pre-GWAS-SSF / unknown
    header_columns: list[str] = field(default_factory=list)
    has_beta: bool | None = None
    has_or: bool | None = None
    has_se: bool | None = None
    has_eaf: bool | None = None
    has_n: bool | None = None
    has_alleles: bool | None = None
    probe_status: str = "not_probed"
    probe_bytes: int = 0

    def to_row(self) -> dict[str, Any]:
        return asdict(self)
