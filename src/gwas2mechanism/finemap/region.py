"""Regional extraction and GWAS <-> reference-panel allele harmonisation.

Vectorised version of the legacy Step03.2/Step04 matching: variants are joined
on (CHR, POS) and kept only when the unordered GWAS allele pair equals the
panel's {REF, ALT}. Effects are then oriented to the panel ALT allele:

    EA == ALT -> BETA_ALT =  BETA, Z_ALT =  Z, EAF_ALT = EAF
    EA == REF -> BETA_ALT = -BETA, Z_ALT = -Z, EAF_ALT = 1 - EAF

Strand-flipped or allele-mismatched variants are dropped (not guessed).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gwas2mechanism.utils.variants import allele_pair_key


def extract_region(sumstats: Path, chrom: str, start: int, end: int) -> pl.DataFrame:
    return (
        pl.scan_parquet(sumstats)
        .filter((pl.col("CHR") == chrom) & pl.col("POS").is_between(start, end))
        .collect()
    )


def match_to_reference(gwas: pl.DataFrame, reference: pl.DataFrame) -> pl.DataFrame:
    """Return matched, ALT-oriented variants keyed by the panel ``VARIANT_ID``."""
    if gwas.is_empty() or reference.is_empty():
        return pl.DataFrame(
            schema={
                "VARIANT_ID": pl.String, "CHR": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String,
                "RSID": pl.String, "EA": pl.String, "NEA": pl.String, "BETA": pl.Float64, "SE": pl.Float64,
                "P": pl.Float64, "EAF": pl.Float64, "N": pl.Float64, "BETA_ALT": pl.Float64, "Z_ALT": pl.Float64,
                "EAF_ALT": pl.Float64, "ORIENTATION": pl.String,
            }
        )
    ref = reference.select(
        pl.col("CHR"), pl.col("POS"), pl.col("ID").alias("VARIANT_ID"), pl.col("REF"), pl.col("ALT"),
    ).with_columns(PAIR=allele_pair_key("REF", "ALT"))
    g = gwas.drop([c for c in ("REF", "ALT", "VARIANT_ID") if c in gwas.columns]).with_columns(
        PAIR=allele_pair_key("EA", "NEA")
    )
    joined = g.join(ref, on=["CHR", "POS", "PAIR"], how="inner")
    flip = pl.col("EA") == pl.col("REF")
    z = pl.col("BETA") / pl.col("SE")
    out = joined.with_columns(
        BETA_ALT=pl.when(flip).then(-pl.col("BETA")).otherwise(pl.col("BETA")),
        Z_ALT=pl.when(flip).then(-z).otherwise(z),
        EAF_ALT=pl.when(flip).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")),
        ORIENTATION=pl.when(flip).then(pl.lit("EA=REF_FLIPPED")).otherwise(pl.lit("EA=ALT")),
    ).filter(pl.col("Z_ALT").is_finite())
    out = out.sort("P").unique(subset=["VARIANT_ID"], keep="first", maintain_order=True)
    keep = [
        "VARIANT_ID", "CHR", "POS", "REF", "ALT", "RSID", "EA", "NEA", "BETA", "SE", "P", "EAF", "N",
        "N_CASES", "N_CONTROLS", "N_EFF", "BETA_ALT", "Z_ALT", "EAF_ALT", "ORIENTATION",
    ]
    return out.select([c for c in keep if c in out.columns]).sort("POS")
