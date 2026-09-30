"""SpliceAI input preparation, execution and VCF parsing."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from gwas2mechanism.annotation.vep import write_vep_vcf
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


def _text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def run_spliceai(variants: pl.DataFrame, output_vcf: Path, fasta: Path, annotation: str, tools: ToolResolver, distance: int = 50, mask: int = 0) -> Path:
    unique = variants.unique("VARIANT_ID", keep="first")
    input_vcf = output_vcf.with_name("spliceai_input.vcf")
    write_vep_vcf(unique, input_vcf)
    run([*tools.require("spliceai"), "-I", str(input_vcf), "-O", str(output_vcf), "-R", str(fasta), "-A", annotation, "-D", str(distance), "-M", str(mask)], log_file=output_vcf.with_suffix(".log"))
    return output_vcf


def parse_spliceai_info(info: str) -> list[dict[str, object]]:
    field = next((part.split("=", 1)[1] for part in info.split(";") if part.startswith("SpliceAI=")), "")
    rows = []
    for item in field.split(","):
        parts = item.split("|")
        if len(parts) < 10:
            continue
        values = [float(value) if value not in {"", "."} else None for value in parts[2:6]]
        rows.append({"ALLELE": parts[0], "GENE": parts[1], "DS_AG": values[0], "DS_AL": values[1], "DS_DG": values[2], "DS_DL": values[3], "DP_AG": parts[6], "DP_AL": parts[7], "DP_DG": parts[8], "DP_DL": parts[9], "SPLICEAI_MAX": max((value for value in values if value is not None), default=None)})
    return rows


def parse_spliceai_vcf(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows: list[dict[str, object]] = []
    with _text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, variant_id, ref, alt, _qual, _filter, info, *_ = line.rstrip().split("\t")
            key = variant_id if variant_id not in {"", "."} else f"{chrom}:{pos}:{ref}:{alt}"
            for record in parse_spliceai_info(info):
                rows.append({"VARIANT_ID": key, **record})
    all_predictions = pl.DataFrame(rows) if rows else pl.DataFrame(schema={"VARIANT_ID": pl.String, "GENE": pl.String, "SPLICEAI_MAX": pl.Float64})
    summary = all_predictions.group_by("VARIANT_ID").agg(pl.max("SPLICEAI_MAX"), pl.col("GENE").drop_nulls().unique().sort()) if rows else all_predictions
    return all_predictions, summary


def classify_score(score: float | None) -> str:
    if score is None:
        return "UNSCORED"
    if score >= 0.8:
        return "HIGH_PRECISION"
    if score >= 0.5:
        return "RECOMMENDED"
    if score >= 0.2:
        return "HIGH_RECALL"
    return "LOW"
