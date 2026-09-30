"""Framework-wide constants. Nothing here is phenotype-specific."""

from __future__ import annotations

#: 1000 Genomes super-populations with dedicated LD references.
SUPERPOPULATIONS: tuple[str, ...] = ("EUR", "AFR", "EAS", "SAS", "AMR")

#: A study whose discovery sample spans several super-populations.
MULTI = "MULTI"
#: The ancestry could not be established from metadata.
UNKNOWN_POPULATION = "UNKNOWN"
#: A population the Catalog reports but 1000 Genomes does not represent.
UNMAPPED_POPULATION = "UNMAPPED"

#: Analysis labels for combined analyses.
COMBINED_PUBLISHED = "COMBINED_PUBLISHED"
COMBINED_META = "COMBINED_META"
MULTIANCESTRY = "MULTIANCESTRY"

ALL_POPULATION_LABELS: tuple[str, ...] = (
    *SUPERPOPULATIONS,
    MULTI,
    UNKNOWN_POPULATION,
    UNMAPPED_POPULATION,
)

#: Canonical internal summary-statistics schema (order preserved on export).
SUMSTATS_COLUMNS: tuple[str, ...] = (
    "STUDY_ID",
    "TRAIT",
    "POPULATION",
    "CHR",
    "POS",
    "RSID",
    "REF",
    "ALT",
    "EA",
    "NEA",
    "BETA",
    "SE",
    "OR",
    "Z",
    "P",
    "NEG_LOG10_P",
    "EAF",
    "N",
    "N_CASES",
    "N_CONTROLS",
    "N_EFF",
    "BUILD",
    "SOURCE_FILE",
    "VARIANT_ID",
)

#: Ordered VEP functional classes (kept from the legacy Step05 logic).
FUNCTIONAL_CATEGORIES: tuple[str, ...] = (
    "SPLICING",
    "CODING",
    "UTR",
    "REGULATORY",
    "NONCODING_RNA",
    "INTRONIC",
    "UPSTREAM",
    "DOWNSTREAM",
    "INTERGENIC",
    "OTHER",
)

GTEX_QTL_TYPES: dict[str, str] = {"eqtl": "eQTL", "sqtl": "sQTL", "apaqtl": "apaQTL"}

WORDING = {
    "variant": "candidate causal variant",
    "coloc": "colocalisation evidence (fine-mapping PIP screen)",
    "splice": "predicted splice effect",
    "regulatory": "candidate regulatory mechanism",
}
