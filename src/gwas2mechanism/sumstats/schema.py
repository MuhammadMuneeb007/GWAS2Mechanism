"""Canonical summary-statistics schema and input column recognition."""

from __future__ import annotations

import polars as pl

from gwas2mechanism.constants import SUMSTATS_COLUMNS

CANONICAL_DTYPES: dict[str, pl.DataType] = {
    "STUDY_ID": pl.String,
    "TRAIT": pl.String,
    "POPULATION": pl.String,
    "CHR": pl.String,
    "POS": pl.Int64,
    "RSID": pl.String,
    "REF": pl.String,
    "ALT": pl.String,
    "EA": pl.String,
    "NEA": pl.String,
    "BETA": pl.Float64,
    "SE": pl.Float64,
    "OR": pl.Float64,
    "Z": pl.Float64,
    "P": pl.Float64,
    "NEG_LOG10_P": pl.Float64,
    "EAF": pl.Float64,
    "N": pl.Float64,
    "N_CASES": pl.Float64,
    "N_CONTROLS": pl.Float64,
    "N_EFF": pl.Float64,
    "BUILD": pl.String,
    "SOURCE_FILE": pl.String,
    "VARIANT_ID": pl.String,
}
assert tuple(CANONICAL_DTYPES) == SUMSTATS_COLUMNS

# Input header aliases (lower-case) -> canonical field. Order = priority:
# harmonised GWAS Catalog columns (hm_*) win over submitted ones.
COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "CHR": ("hm_chrom", "chromosome", "chrom", "chr", "#chrom", "#chr", "chromosome_name", "ch"),
    "POS": ("hm_pos", "base_pair_location", "pos", "bp", "position", "base_pair", "genpos", "bp_hg38", "pos_b38"),
    "RSID": ("hm_rsid", "rsid", "rs_id", "snp", "markername", "variant_id", "snpid", "rsids", "id"),
    "EA": ("hm_effect_allele", "effect_allele", "ea", "a1", "alt", "allele1", "tested_allele", "effectallele"),
    "NEA": ("hm_other_allele", "other_allele", "nea", "a2", "ref", "allele0", "allele2", "non_effect_allele", "otherallele"),
    "BETA": ("hm_beta", "beta", "b", "effect", "log_odds", "logor", "effect_size", "beta_fixed"),
    "OR": ("hm_odds_ratio", "odds_ratio", "or"),
    "OR_L95": ("hm_ci_lower", "ci_lower", "or_95l", "lower_ci", "or_l95"),
    "OR_U95": ("hm_ci_upper", "ci_upper", "or_95u", "upper_ci", "or_u95"),
    "SE": ("standard_error", "se", "stderr", "sebeta", "se_beta", "standard_error_of_beta", "se_fixed"),
    "Z": ("z", "zscore", "z_score", "zstat", "z_stat"),
    "P": ("p_value", "p", "pval", "pvalue", "p-value", "p.value", "p_bolt_lmm", "p_fixed", "p_value_nominal"),
    "NEG_LOG10_P": ("neg_log_10_p_value", "neg_log10_p", "log10p", "mlog10p", "minus_log10_pval", "lp"),
    "EAF": ("hm_effect_allele_frequency", "effect_allele_frequency", "eaf", "af", "freq", "af_alt", "a1freq", "maf_effect", "frq", "freq1"),
    "N": ("n", "n_total", "sample_size", "neff", "n_samples", "totalsamplesize", "nobs"),
    "N_CASES": ("n_cases", "ncase", "n_case", "cases", "num_cases"),
    "N_CONTROLS": ("n_controls", "ncontrol", "n_control", "controls", "num_controls"),
    "INFO": ("info", "imputation_score", "r2", "rsq", "impinfo"),
}

REQUIRED_FIELDS = ("CHR", "POS", "EA", "NEA")


def map_columns(header: list[str]) -> dict[str, str]:
    """Return ``{canonical: source_column}`` for a raw header (first alias wins)."""
    lowered = {h.strip().lower(): h for h in header}
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            source = lowered.get(alias)
            if source is not None and source not in used:
                mapping[canonical] = source
                used.add(source)
                break
    return mapping


def capabilities(mapping: dict[str, str]) -> dict[str, bool]:
    return {
        "has_alleles": "EA" in mapping and "NEA" in mapping,
        "has_beta": "BETA" in mapping,
        "has_or": "OR" in mapping,
        "has_se": "SE" in mapping or ("OR_L95" in mapping and "OR_U95" in mapping),
        "has_eaf": "EAF" in mapping,
        "has_n": "N" in mapping,
        "has_p": "P" in mapping or "NEG_LOG10_P" in mapping,
    }


def empty_sumstats() -> pl.DataFrame:
    return pl.DataFrame(schema=CANONICAL_DTYPES)


def validate(frame: pl.DataFrame | pl.LazyFrame) -> list[str]:
    schema = frame.collect_schema() if isinstance(frame, pl.LazyFrame) else frame.schema
    problems = [f"missing column {c}" for c in SUMSTATS_COLUMNS if c not in schema]
    for column, dtype in CANONICAL_DTYPES.items():
        if column in schema and schema[column] != dtype:
            problems.append(f"{column} has dtype {schema[column]}, expected {dtype}")
    return problems
