"""Keep distinct genomic statements that the legacy experiments validated."""

from __future__ import annotations

import polars as pl


def classify_lncrna_context(overlaps: pl.DataFrame) -> pl.DataFrame:
    """Add non-equivalent evidence flags without collapsing them into one label."""
    required = {"VARIANT_ID", "OVERLAPS_PROTEIN_CODING_CDS", "OVERLAPS_LNCRNA_GENE_BODY"}
    missing = required - set(overlaps.columns)
    if missing:
        raise ValueError(f"lncRNA overlap table lacks {sorted(missing)}")
    return overlaps.with_columns(
        (~pl.col("OVERLAPS_PROTEIN_CODING_CDS")).alias("OUTSIDE_PROTEIN_CODING_CDS"),
        pl.col("OVERLAPS_LNCRNA_GENE_BODY").alias("WITHIN_LNCRNA_GENE_BODY"),
    )
