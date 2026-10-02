"""Read-only resource inspection and post-setup validation."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from gwas2mechanism.config import Config
from gwas2mechanism.constants import SUPERPOPULATIONS
from gwas2mechanism.provenance import ResourceManifest
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.utils.io import read_json
from gwas2mechanism.utils.tools import ToolResolver


@dataclass(frozen=True)
class StatusRow:
    resource: str
    status: str
    detail: str = ""


def _file_ok(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _ratio_status(done: int, total: int) -> str:
    if done == total and total:
        return "COMPLETE"
    if done:
        return "PARTIAL"
    return "MISSING"


def inspect_resources(config: Config) -> dict[str, object]:
    root = config.resources_path.resolve()
    tools = ToolResolver(config.tools)
    panel = ReferencePanel(config.reference, root, tools)
    chromosomes = config.genome.chromosomes
    genome = root / "genome" / "GRCh38"
    fasta = genome / "GRCh38.primary_assembly.fa"
    gtf = genome / "gencode.v50.annotation.gtf"
    rows = [
        StatusRow("GRCh38 FASTA", "COMPLETE" if _file_ok(fasta) else "MISSING", str(fasta)),
        StatusRow("GENCODE 50 annotation", "COMPLETE" if _file_ok(gtf) else "MISSING", str(gtf)),
        StatusRow("GRCh38 FASTA index", "COMPLETE" if _file_ok(Path(f"{fasta}.fai")) else "MISSING"),
    ]
    raw_done = sum(panel.raw_chromosome_complete(chrom) for chrom in chromosomes)
    all_done = sum(panel.all_chromosome_complete(chrom) for chrom in chromosomes)
    rows.append(StatusRow("1000G raw", _ratio_status(raw_done, len(chromosomes)), f"{raw_done}/{len(chromosomes)}"))
    rows.append(StatusRow("1000G ALL PGEN (intermediate)", _ratio_status(all_done, len(chromosomes)), f"{all_done}/{len(chromosomes)}"))
    for population in SUPERPOPULATIONS:
        done = sum(panel.population_chromosome_complete(population, chrom) for chrom in chromosomes)
        rows.append(StatusRow(population, _ratio_status(done, len(chromosomes)), f"{done}/{len(chromosomes)}"))

    gtex_root = root / "gtex"
    gtex_files = [path for path in gtex_root.rglob("*") if path.is_file() and not path.name.endswith((".part", ".json"))] if gtex_root.exists() else []
    for label, token in (("GTEx SuSiE eQTL", "eqtl"), ("GTEx SuSiE sQTL", "sqtl")):
        matching = [path for path in gtex_files if token in path.name.lower()]
        rows.append(StatusRow(label, "COMPLETE" if matching else "MISSING", f"{len(matching)} files" if matching else ""))

    vep = root / "vep"
    vep_files = [path for path in vep.rglob("*") if path.is_file() and path.stat().st_size > 0 and path.name not in {"install.log", ".complete.json"}] if vep.exists() else []
    rows.append(StatusRow("VEP cache", "COMPLETE" if vep_files else "MISSING", str(vep)))
    rows.append(StatusRow("SpliceAI", "COMPLETE" if tools.available("spliceai") else "MISSING", tools.version("spliceai") or ""))
    pangolin_db = root / "splicing" / "pangolin" / "gencode.annotation.db"
    rows.append(StatusRow("Pangolin annotation DB", "COMPLETE" if _file_ok(pangolin_db) and tools.available("pangolin") else "MISSING", str(pangolin_db)))

    counts = {name: sum(row.status == name for row in rows) for name in ("COMPLETE", "PARTIAL", "MISSING")}
    total = len(rows)
    return {
        "resource_cache": str(root),
        "rows": rows,
        "total": total,
        "complete": counts["COMPLETE"],
        "partial": counts["PARTIAL"],
        "missing": counts["MISSING"],
        "percentage_complete": round(100 * counts["COMPLETE"] / total, 1) if total else 100.0,
    }


def validate_resources(config: Config, *, reference: bool = True, gtex: bool = True, vep: bool = True, splice: bool = True, populations: list[str] | None = None) -> list[str]:
    root = config.resources_path.resolve()
    tools = ToolResolver(config.tools)
    errors: list[str] = []
    genome = root / "genome" / "GRCh38"
    if reference or splice:
        for path in (genome / "GRCh38.primary_assembly.fa", genome / "GRCh38.primary_assembly.fa.fai", genome / "gencode.v50.annotation.gtf"):
            if not _file_ok(path):
                errors.append(f"missing or empty: {path}")
    if reference:
        panel = ReferencePanel(config.reference, root, tools)
        for chrom in config.genome.chromosomes:
            if not panel.raw_chromosome_complete(chrom):
                errors.append(f"1000G raw chr{chrom} incomplete")
            if not panel.all_chromosome_complete(chrom):
                errors.append(f"1000G ALL intermediate chr{chrom} incomplete")
            for population in populations or config.populations.allowed:
                if not panel.population_chromosome_complete(population, chrom):
                    errors.append(f"1000G {population} chr{chrom} incomplete")
    status = inspect_resources(config)
    by_name = {row.resource: row.status for row in status["rows"]}
    if gtex:
        for name in ("GTEx SuSiE eQTL", "GTEx SuSiE sQTL"):
            if by_name.get(name) != "COMPLETE":
                errors.append(f"{name} incomplete")
        manifest = ResourceManifest(root / "manifest.tsv").read()
        if not manifest.is_empty():
            objects = manifest.filter(pl.col("resource").str.starts_with("gtex/") & ~pl.col("resource").str.ends_with("/susie"))
            for row in objects.iter_rows(named=True):
                path = Path(row["local_path"])
                if not _file_ok(path):
                    errors.append(f"GTEx object missing or empty: {path}")
                if path.name.endswith((".tar", ".tar.gz", ".tgz")):
                    stem = path.name.removesuffix(".gz").removesuffix(".tar").removesuffix(".tgz")
                    marker = path.parent / "extracted" / stem / ".complete.json"
                    if not _file_ok(marker):
                        errors.append(f"GTEx archive not successfully extracted: {path}")
    if vep and by_name.get("VEP cache") != "COMPLETE":
        errors.append("VEP cache incomplete")
    elif vep:
        marker = root / "vep" / ".complete.json"
        installed = tools.version("vep")
        try:
            cached = read_json(marker).get("version")
        except (OSError, ValueError):
            cached = None
        if installed and cached and str(installed).split(".", 1)[0] != str(cached).split(".", 1)[0]:
            errors.append(f"VEP cache/software release mismatch: cache {cached}, executable {installed}")
    if splice:
        if by_name.get("SpliceAI") != "COMPLETE":
            errors.append("SpliceAI executable unavailable")
        if by_name.get("Pangolin annotation DB") != "COMPLETE":
            errors.append("Pangolin executable or annotation DB unavailable")
        for name in ("spliceai", "pangolin"):
            command = tools.command(name)
            if command:
                check = subprocess.run([*command, "--help"], capture_output=True, text=True, timeout=120)
                if check.returncode:
                    errors.append(f"{name} executable check failed")
    rscript = tools.command("rscript")
    if rscript:
        check = subprocess.run([*rscript, "-e", "stopifnot(all(vapply(c('susieR','coloc','data.table'), requireNamespace, logical(1), quietly=TRUE)))"], capture_output=True, text=True, timeout=120)
        if check.returncode:
            errors.append("R packages susieR, coloc and/or data.table unavailable")
    else:
        errors.append("Rscript unavailable")
    return errors


def human_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024 or unit == "TiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024
    return f"{amount:.1f} TiB"
