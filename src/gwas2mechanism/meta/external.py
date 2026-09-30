"""Optional external meta-analysis tools: METAL and MR-MEGA.

Both are *optional*. The fast mode never requires them. Inputs are written from
the harmonised Parquet files; outputs are parsed back into Polars.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

import polars as pl

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver

log = logging.getLogger(__name__)


def _export(path: Path, out: Path) -> Path:
    (
        pl.scan_parquet(path)
        .select("VARIANT_ID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "EAF", "N")
        .with_columns(pl.col("N").fill_null(0))
        .sink_csv(out, separator="\t")
    )
    return out


def run_metal(sumstats: Mapping[str, Path], workdir: Path, tools: ToolResolver) -> pl.DataFrame | None:
    cmd = tools.command("metal")
    if cmd is None:
        log.info("METAL not installed; skipping external fixed-effect meta-analysis")
        return None
    workdir.mkdir(parents=True, exist_ok=True)
    lines = [
        "SCHEME STDERR",
        "AVERAGEFREQ ON",
        "MINMAXFREQ ON",
        "MARKER VARIANT_ID",
        "ALLELE EA NEA",
        "EFFECT BETA",
        "STDERR SE",
        "PVALUE P",
        "FREQLABEL EAF",
        "WEIGHT N",
    ]
    for pop, path in sumstats.items():
        lines.append(f"PROCESS {_export(path, workdir / f'{pop}.tsv')}")
    lines += [f"OUTFILE {workdir / 'metal'} .tbl", "ANALYZE HETEROGENEITY", "QUIT"]
    script = workdir / "metal.txt"
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run([*cmd, str(script)], log_file=workdir / "metal.log", cwd=workdir)
    table = workdir / "metal1.tbl"
    return pl.read_csv(table, separator="\t") if table.exists() else None


def run_mrmega(sumstats: Mapping[str, Path], workdir: Path, tools: ToolResolver, n_pcs: int = 1) -> pl.DataFrame | None:
    """MR-MEGA meta-regression (requires >= n_pcs + 3 populations)."""
    cmd = tools.command("mrmega")
    if cmd is None or len(sumstats) < n_pcs + 3:
        log.info("MR-MEGA skipped (tool missing or too few populations)")
        return None
    workdir.mkdir(parents=True, exist_ok=True)
    listing = []
    for pop, path in sumstats.items():
        out = workdir / f"{pop}.mrmega.tsv"
        (
            pl.scan_parquet(path)
            .filter(pl.col("EAF").is_not_null() & pl.col("N").is_not_null())
            .select(
                pl.col("VARIANT_ID").alias("MARKERNAME"),
                pl.col("EA"),
                pl.col("NEA"),
                pl.col("EAF"),
                pl.col("BETA"),
                pl.col("SE"),
                pl.col("N"),
                pl.col("CHR").alias("CHROMOSOME"),
                pl.col("POS").alias("POSITION"),
            )
            .sink_csv(out, separator="\t")
        )
        listing.append(str(out))
    in_file = workdir / "mr-mega.in"
    in_file.write_text("\n".join(listing) + "\n", encoding="utf-8")
    run([*cmd, "-i", str(in_file), "--pc", str(n_pcs), "-o", str(workdir / "mrmega")], log_file=workdir / "mrmega.log")
    result = workdir / "mrmega.result"
    return pl.read_csv(result, separator="\t", null_values=["NA"]) if result.exists() else None
