"""PIP-based colocalisation screen (explicitly not formal coloc.susie)."""

from __future__ import annotations

import numpy as np
import polars as pl


def variant_clpp(gwas_pip: np.ndarray | float, qtl_pip: np.ndarray | float) -> np.ndarray:
    return np.asarray(gwas_pip, dtype=np.float64) * np.asarray(qtl_pip, dtype=np.float64)


def locus_clpp(values: np.ndarray | list[float]) -> float:
    array = np.clip(np.asarray(values, dtype=np.float64), 0.0, 1.0)
    if array.size == 0:
        return 0.0
    if np.any(array == 1.0):
        return 1.0
    return float(-np.expm1(np.log1p(-array).sum()))


def screen(gwas: pl.DataFrame, qtl: pl.DataFrame, cross_ancestry: pl.DataFrame | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Exact normalized-allele match and vCLPP/LCLPP aggregation."""
    pip_col = "PIP"
    left = gwas.select("VARIANT_ID", "LOCUS_ID", "POPULATION", pl.col(pip_col).alias("GWAS_PIP"))
    if cross_ancestry is not None and not cross_ancestry.is_empty():
        left = left.join(cross_ancestry.select("VARIANT_ID", "CROSS_ANCESTRY_PIP").unique("VARIANT_ID"), on="VARIANT_ID", how="left")
    overlap = left.join(qtl, on="VARIANT_ID", how="inner").with_columns((pl.col("GWAS_PIP") * pl.col("QTL_PIP")).alias("vCLPP"))
    if overlap.is_empty():
        return overlap, pl.DataFrame(schema={"LOCUS_ID": pl.String, "POPULATION": pl.String, "GENE": pl.String, "TISSUE": pl.String, "QTL_TYPE": pl.String, "LCLPP": pl.Float64})
    pairs = overlap.group_by("LOCUS_ID", "POPULATION", "GENE", "TISSUE", "QTL_TYPE").agg(
        pl.col("vCLPP").map_batches(lambda s: pl.Series([locus_clpp(s.to_numpy())]), returns_scalar=True).alias("LCLPP"),
        pl.max("vCLPP").alias("MAX_vCLPP"),
        pl.len().alias("N_MATCHING_VARIANTS"),
        pl.col("VARIANT_ID").sort_by("vCLPP", descending=True).first().alias("TOP_VARIANT"),
    )
    return overlap, pairs
