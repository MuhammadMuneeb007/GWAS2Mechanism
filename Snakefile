"""GWAS2Mechanism restartable workflow entry point."""

from pathlib import Path

RUN_DIR = str(Path(config.get("run_dir", "runs/unconfigured")).resolve())
PHENOTYPE = config.get("phenotype")
if not PHENOTYPE:
    raise ValueError("Snakemake requires --config phenotype=<trait> run_dir=<path>")

include: "workflow/rules/discover.smk"
include: "workflow/rules/download.smk"
include: "workflow/rules/harmonize.smk"
include: "workflow/rules/references.smk"
include: "workflow/rules/clump.smk"
include: "workflow/rules/finemap.smk"
include: "workflow/rules/multiancestry.smk"
include: "workflow/rules/annotate.smk"
include: "workflow/rules/qtl.smk"
include: "workflow/rules/splice.smk"
include: "workflow/rules/coloc.smk"
include: "workflow/rules/report.smk"

rule all:
    input:
        f"{RUN_DIR}/report/report.html",
        f"{RUN_DIR}/report/report.md",
        f"{RUN_DIR}/report/summary.tsv",
