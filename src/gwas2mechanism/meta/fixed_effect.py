"""Inverse-variance fixed-effect meta-analysis of ancestry-stratified GWAS.

For each variant with estimates from k populations (after allele alignment):

    w_i      = 1 / SE_i^2
    beta_meta = sum(w_i beta_i) / sum(w_i)
    SE_meta   = sqrt(1 / sum(w_i))
    Z_meta    = beta_meta / SE_meta
    P_meta    = 2 * norm.sf(|Z_meta|)
    Q         = sum(w_i (beta_i - beta_meta)^2)       (Cochran's Q, df = k - 1)
    I^2       = max(0, (Q - df) / Q)

Alleles are aligned to a common effect allele (the lexicographically larger
allele of the unordered pair) before combining; incompatible allele pairs never
share a key and are therefore never combined. Heterogeneity is *retained and
reported*, never hidden.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from gwas2mechanism.constants import COMBINED_META, SUMSTATS_COLUMNS
from gwas2mechanism.sumstats.schema import CANONICAL_DTYPES
from gwas2mechanism.utils.variants import is_palindromic

log = logging.getLogger(__name__)
LN10 = math.log(10.0)


class IncompatibleScalesError(ValueError):
    pass


def align_for_meta(lf: pl.LazyFrame, population: str, drop_palindromic_maf: float | None) -> pl.LazyFrame:
    """Orient each study to effect allele = A2 of the sorted allele pair."""
    parts = pl.col("VARIANT_ID").str.split_exact(":", 3)
    a2 = parts.struct.field("field_3")
    flip = pl.col("EA") != a2
    out = lf.select(
        "VARIANT_ID",
        "CHR",
        "POS",
        parts.struct.field("field_2").alias("A1"),
        a2.alias("A2"),
        pl.when(flip).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
        pl.col("SE"),
        pl.when(flip).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")).alias("EAF"),
        pl.col("N"),
        pl.col("N_CASES"),
        pl.col("N_CONTROLS"),
        pl.col("RSID"),
        pl.lit(population).alias("POP"),
    ).filter(pl.col("SE").is_finite() & (pl.col("SE") > 0) & pl.col("BETA").is_finite())
    if drop_palindromic_maf is not None:
        maf = pl.min_horizontal(pl.col("EAF"), 1.0 - pl.col("EAF"))
        ambiguous = is_palindromic("A1", "A2") & (pl.col("EAF").is_null() | (maf > drop_palindromic_maf))
        out = out.filter(~ambiguous)
    return out


def inverse_variance_meta(long: pl.DataFrame, populations: list[str]) -> pl.DataFrame:
    """Vectorised per-variant IVW meta-analysis on a long (variant x population) frame."""
    agg = (
        long.with_columns(W=1.0 / pl.col("SE") ** 2)
        .group_by("VARIANT_ID")
        .agg(
            pl.col("CHR").first(),
            pl.col("POS").first(),
            pl.col("A1").first(),
            pl.col("A2").first(),
            pl.col("RSID").drop_nulls().first(),
            pl.len().alias("N_POPULATIONS"),
            pl.col("POP").sort().str.join(",").alias("POPULATIONS"),
            pl.col("W").sum().alias("SUM_W"),
            (pl.col("W") * pl.col("BETA")).sum().alias("SUM_WB"),
            (pl.col("W") * pl.col("BETA") ** 2).sum().alias("SUM_WB2"),
            (pl.col("W") * pl.col("EAF")).sum().alias("SUM_W_EAF"),
            pl.col("W").filter(pl.col("EAF").is_not_null()).sum().alias("SUM_W_EAF_OBS"),
            pl.col("N").sum().alias("N"),
            pl.col("N_CASES").sum().alias("N_CASES"),
            pl.col("N_CONTROLS").sum().alias("N_CONTROLS"),
            *[
                pl.col("BETA").filter(pl.col("POP") == p).first().alias(f"BETA_{p}")
                for p in populations
            ],
            *[pl.col("SE").filter(pl.col("POP") == p).first().alias(f"SE_{p}") for p in populations],
        )
    )
    agg = agg.with_columns(
        BETA=pl.col("SUM_WB") / pl.col("SUM_W"),
        SE=(1.0 / pl.col("SUM_W")).sqrt(),
        Q=(pl.col("SUM_WB2") - pl.col("SUM_WB") ** 2 / pl.col("SUM_W")).clip(lower_bound=0.0),
        EAF=pl.when(pl.col("SUM_W_EAF_OBS") > 0).then(pl.col("SUM_W_EAF") / pl.col("SUM_W_EAF_OBS")),
    ).with_columns(Z=pl.col("BETA") / pl.col("SE"))

    z = agg["Z"].to_numpy().astype(np.float64)
    q = agg["Q"].to_numpy().astype(np.float64)
    k = agg["N_POPULATIONS"].to_numpy().astype(np.float64)
    df = k - 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 2.0 * stats.norm.sf(np.abs(z))
        nlp = -(math.log(2.0) + stats.norm.logsf(np.abs(z))) / LN10
        q_p = np.where(df > 0, stats.chi2.sf(q, np.maximum(df, 1.0)), np.nan)
        i2 = np.where((df > 0) & (q > 0), np.clip((q - df) / q, 0.0, 1.0), np.where(df > 0, 0.0, np.nan))
    direction = pl.concat_str(
        [
            pl.when(pl.col(f"BETA_{p}").is_null())
            .then(pl.lit("?"))
            .when(pl.col(f"BETA_{p}") > 0)
            .then(pl.lit("+"))
            .when(pl.col(f"BETA_{p}") < 0)
            .then(pl.lit("-"))
            .otherwise(pl.lit("0"))
            for p in populations
        ]
    )
    return agg.with_columns(
        P=pl.Series(p),
        NEG_LOG10_P=pl.Series(nlp),
        Q_DF=pl.Series(df),
        Q_P=pl.Series(q_p),
        I2=pl.Series(i2),
        DIRECTION=direction,
    ).drop("SUM_W", "SUM_WB", "SUM_WB2", "SUM_W_EAF", "SUM_W_EAF_OBS")


def run_meta(
    sumstats: Mapping[str, Path],
    trait_types: Mapping[str, str | None],
    out_path: Path,
    trait: str,
    drop_palindromic_maf: float | None = 0.4,
) -> pl.DataFrame:
    if len(sumstats) < 2:
        raise ValueError("fixed-effect meta-analysis requires >= 2 ancestry-specific datasets")
    types = {t for t in trait_types.values() if t}
    if len({("binary" if t == "binary" else "other") for t in types}) > 1:
        raise IncompatibleScalesError(
            f"populations report different effect scales ({sorted(types)}); refusing to combine"
        )
    populations = sorted(sumstats)
    frames = [
        align_for_meta(pl.scan_parquet(path), pop, drop_palindromic_maf) for pop, path in sumstats.items()
    ]
    long = pl.concat(frames, how="vertical_relaxed").collect()
    meta = inverse_variance_meta(long, populations)
    out = meta.with_columns(
        STUDY_ID=pl.lit(COMBINED_META),
        TRAIT=pl.lit(trait),
        POPULATION=pl.lit(COMBINED_META),
        REF=pl.lit(None, pl.String),
        ALT=pl.lit(None, pl.String),
        EA=pl.col("A2"),
        NEA=pl.col("A1"),
        OR=pl.lit(None, pl.Float64),
        N_EFF=pl.when((pl.col("N_CASES") > 0) & (pl.col("N_CONTROLS") > 0)).then(
            4.0 / (1.0 / pl.col("N_CASES") + 1.0 / pl.col("N_CONTROLS"))
        ),
        BUILD=pl.lit("GRCh38"),
        SOURCE_FILE=pl.lit("fixed_effect_ivw:" + ",".join(populations)),
    )
    extra = ["N_POPULATIONS", "POPULATIONS", "DIRECTION", "Q", "Q_DF", "Q_P", "I2"] + [
        c for c in out.columns if c.startswith(("BETA_", "SE_"))
    ]
    out = out.select(
        [pl.col(c).cast(CANONICAL_DTYPES[c]) for c in SUMSTATS_COLUMNS] + extra
    ).sort(["CHR", "POS"])
    from gwas2mechanism.utils.io import write_parquet

    write_parquet(out, out_path)
    log.info("COMBINED_META: %d variants from %s", out.height, populations)
    return out
