"""LD clumping adapters and deterministic locus construction.

PLINK2 remains the production clumping engine.  The small in-memory clumper is
intended for fixtures and for references represented by an explicit LD matrix;
it never pretends that physical distance is LD.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import polars as pl

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


def greedy_ld_clump(
    variants: pl.DataFrame,
    ld: np.ndarray,
    *,
    p1: float = 5e-8,
    p2: float = 1.0,
    r2: float = 0.1,
    kb: int = 1000,
) -> pl.DataFrame:
    """Greedily clump a small, position-aligned table using a real LD matrix."""
    required = {"CHR", "POS", "P", "VARIANT_ID"}
    missing = required - set(variants.columns)
    if missing:
        raise ValueError(f"Clumping input lacks {sorted(missing)}")
    if ld.shape != (variants.height, variants.height):
        raise ValueError("LD dimensions must match the variant table")
    order = np.argsort(variants["P"].to_numpy().astype(float), kind="stable")
    chrom = variants["CHR"].cast(pl.String).to_numpy()
    pos = variants["POS"].to_numpy().astype(np.int64)
    p = variants["P"].to_numpy().astype(float)
    available = np.isfinite(p) & (p <= p2)
    rows: list[dict[str, object]] = []
    for lead in order:
        lead = int(lead)
        if not available[lead] or p[lead] > p1:
            continue
        in_window = (chrom == chrom[lead]) & (np.abs(pos - pos[lead]) <= kb * 1000)
        members = available & in_window & (np.square(ld[lead]) >= r2)
        members[lead] = True
        member_idx = np.flatnonzero(members).astype(np.int64)
        rows.append(
            {
                "LEAD_VARIANT": variants["VARIANT_ID"][lead],
                "CHR": str(chrom[lead]),
                "LEAD_POS": int(pos[lead]),
                "LEAD_P": float(p[lead]),
                "N_MEMBERS": int(member_idx.size),
                "MEMBERS": variants["VARIANT_ID"].gather(pl.Series(member_idx)).to_list(),
            }
        )
        available[members] = False
    return pl.DataFrame(rows) if rows else pl.DataFrame(
        schema={
            "LEAD_VARIANT": pl.String,
            "CHR": pl.String,
            "LEAD_POS": pl.Int64,
            "LEAD_P": pl.Float64,
            "N_MEMBERS": pl.Int64,
            "MEMBERS": pl.List(pl.String),
        }
    )


def construct_loci(
    leads: pl.DataFrame,
    *,
    population: str,
    window_bp: int = 1_000_000,
    merge_gap_bp: int = 0,
) -> pl.DataFrame:
    """Expand lead positions and merge overlapping intervals by chromosome."""
    if leads.is_empty():
        return pl.DataFrame(
            schema={
                "LOCUS_ID": pl.String,
                "POPULATION": pl.String,
                "CHR": pl.String,
                "START": pl.Int64,
                "END": pl.Int64,
                "LEAD_VARIANTS": pl.List(pl.String),
                "LEAD_P": pl.Float64,
            }
        )
    pos_col = "LEAD_POS" if "LEAD_POS" in leads.columns else "POS"
    id_col = "LEAD_VARIANT" if "LEAD_VARIANT" in leads.columns else "VARIANT_ID"
    p_col = "LEAD_P" if "LEAD_P" in leads.columns else "P"
    source = leads.select(
        pl.col("CHR").cast(pl.String),
        pl.col(pos_col).cast(pl.Int64).alias("POS"),
        pl.col(id_col).cast(pl.String).alias("VARIANT_ID"),
        pl.col(p_col).cast(pl.Float64).alias("P"),
    ).sort(["CHR", "POS"])
    out: list[dict[str, object]] = []
    for chrom in source["CHR"].unique(maintain_order=True):
        records = source.filter(pl.col("CHR") == chrom).iter_rows(named=True)
        current: dict[str, object] | None = None
        for record in records:
            start = max(1, int(record["POS"]) - window_bp)
            end = int(record["POS"]) + window_bp
            if current is None or start > int(current["END"]) + merge_gap_bp:
                if current is not None:
                    out.append(current)
                current = {
                    "POPULATION": population,
                    "CHR": str(chrom),
                    "START": start,
                    "END": end,
                    "LEAD_VARIANTS": [record["VARIANT_ID"]],
                    "LEAD_P": float(record["P"]),
                }
            else:
                current["END"] = max(int(current["END"]), end)
                current["LEAD_VARIANTS"].append(record["VARIANT_ID"])  # type: ignore[union-attr]
                current["LEAD_P"] = min(float(current["LEAD_P"]), float(record["P"]))
        if current is not None:
            out.append(current)
    result = pl.DataFrame(out).sort(["CHR", "START"])
    return result.with_columns(
        pl.concat_str(
            [pl.lit(population), pl.lit(":"), pl.col("CHR"), pl.lit(":"), pl.col("START"), pl.lit("-"), pl.col("END")]
        ).alias("LOCUS_ID")
    ).select("LOCUS_ID", "POPULATION", "CHR", "START", "END", "LEAD_VARIANTS", "LEAD_P")


def union_loci(frames: Iterable[pl.DataFrame], merge_gap_bp: int = 0) -> pl.DataFrame:
    """Union significant intervals across populations without using pooled LD."""
    present = [frame for frame in frames if not frame.is_empty()]
    if not present:
        return construct_loci(pl.DataFrame(), population="combined")
    leads = pl.concat(present, how="diagonal_relaxed").explode(
        "LEAD_VARIANTS", empty_as_null=True
    ).select(
        "CHR",
        ((pl.col("START") + pl.col("END")) // 2).alias("LEAD_POS"),
        pl.col("LEAD_VARIANTS").alias("LEAD_VARIANT"),
        "LEAD_P",
    )
    return construct_loci(
        leads,
        population="combined",
        window_bp=0,
        merge_gap_bp=merge_gap_bp,
    )


def run_plink_clump(
    sumstats_tsv: Path,
    reference_prefix: Path,
    output_prefix: Path,
    tools: ToolResolver,
    *,
    p1: float,
    p2: float,
    r2: float,
    kb: int,
    threads: int,
) -> Path:
    """Run PLINK2 clumping against one explicitly selected population panel."""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            *tools.require("plink2"),
            "--pfile",
            str(reference_prefix),
            "--clump",
            str(sumstats_tsv),
            "--clump-id-field",
            "VARIANT_ID",
            "--clump-p-field",
            "P",
            "--clump-p1",
            str(p1),
            "--clump-p2",
            str(p2),
            "--clump-r2",
            str(r2),
            "--clump-kb",
            str(kb),
            "--threads",
            str(threads),
            "--out",
            str(output_prefix),
        ],
        log_file=output_prefix.with_suffix(".log"),
    )
    result = output_prefix.with_suffix(".clumps")
    if not result.exists():
        raise RuntimeError(f"PLINK2 did not create {result}")
    return result
