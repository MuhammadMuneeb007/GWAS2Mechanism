"""Optional formal coloc.susie policy and execution adapter."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


@dataclass(frozen=True)
class FormalColocDecision:
    run: bool
    reason: str


def decide_formal_coloc(*, enabled: bool, dense_gwas: Path | None, dense_qtl: Path | None, qtl_size_bytes: int | None, max_download_mb: int) -> FormalColocDecision:
    if not enabled:
        return FormalColocDecision(False, "disabled by configuration")
    if dense_gwas is None or not dense_gwas.exists():
        return FormalColocDecision(False, "dense regional GWAS statistics unavailable")
    if dense_qtl is None or not dense_qtl.exists():
        if qtl_size_bytes and qtl_size_bytes > max_download_mb * 1024 * 1024:
            return FormalColocDecision(False, "dense QTL resource exceeds configured download limit")
        return FormalColocDecision(False, "dense regional QTL statistics unavailable")
    return FormalColocDecision(True, "sufficient dense regional statistics available")


def run_coloc_susie(
    gwas_fit_rds: Path,
    qtl_fit_rds: Path,
    output_tsv: Path,
    *,
    tools: ToolResolver,
) -> pl.DataFrame:
    """Run ``coloc::coloc.susie`` on two precomputed susieR fit objects.

    This deliberately accepts RDS fits rather than silently rebuilding SuSiE
    models with potentially inconsistent variant order or LD.  Callers must
    first pass :func:`decide_formal_coloc` and align both fits to the same
    variant keys.
    """
    for path in (gwas_fit_rds, qtl_fit_rds):
        if not path.is_file():
            raise FileNotFoundError(path)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    script = output_tsv.with_suffix(".coloc.R")
    script.write_text(
        """args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(coloc))
gwas_fit <- readRDS(args[1])
qtl_fit <- readRDS(args[2])
result <- coloc.susie(gwas_fit, qtl_fit)
write.table(result$summary, file=args[3], sep="\\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
""",
        encoding="utf-8",
    )
    run(
        [*tools.require("rscript"), str(script), str(gwas_fit_rds), str(qtl_fit_rds), str(output_tsv)],
        log_file=output_tsv.with_suffix(".log"),
    )
    return pl.read_csv(output_tsv, separator="\t")
