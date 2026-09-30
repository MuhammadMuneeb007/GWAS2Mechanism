"""Pangolin adapter kept independent from SpliceAI."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


def run_pangolin(
    variants: pl.DataFrame,
    output_csv: Path,
    fasta: Path,
    annotation_db: Path,
    tools: ToolResolver,
    distance: int = 500,
    mask: bool = False,
) -> Path:
    """Run the upstream Pangolin CLI: variants, FASTA, annotation DB, prefix."""
    if not annotation_db.exists():
        raise FileNotFoundError(
            f"Pangolin annotation database is missing: {annotation_db}; run `gwas2m setup --splice`"
        )
    input_csv = output_csv.with_name("pangolin_input.csv")
    input_csv.parent.mkdir(parents=True, exist_ok=True)
    variants.select(
        pl.col("CHR").alias("CHROM"), "POS", "REF", "ALT", pl.col("VARIANT_ID").alias("REFERENCE_ID")
    ).unique("REFERENCE_ID").write_csv(input_csv)
    prefix = output_csv.with_suffix("")
    command = [
        *tools.require("pangolin"),
        "-m",
        str(mask),
        "-d",
        str(distance),
        "-c",
        "CHROM,POS,REF,ALT",
        str(input_csv),
        str(fasta),
        str(annotation_db),
        str(prefix),
    ]
    run(command, log_file=output_csv.with_suffix(".log"))
    produced = prefix.with_suffix(".csv")
    if produced != output_csv and produced.exists():
        produced.replace(output_csv)
    if not output_csv.exists():
        raise RuntimeError(f"Pangolin did not create {output_csv}")
    return output_csv


def parse_pangolin_table(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split(",")
        if "Pangolin" not in header or "REFERENCE_ID" not in header:
            raise ValueError(f"Unexpected Pangolin CSV header: {header}")
        pangolin_index = header.index("Pangolin")
        reference_index = header.index("REFERENCE_ID")
        split_count = len(header) - 1
        for line in handle:
            fields = line.rstrip("\n").split(",", split_count)
            if len(fields) <= max(pangolin_index, reference_index):
                continue
            variant_id = fields[reference_index]
            for prediction in fields[pangolin_index].split(","):
                parts = prediction.split("|")
                if len(parts) < 3:
                    continue
                def score(value: str) -> float | None:
                    try:
                        return float(value.rsplit(":", 1)[-1])
                    except ValueError:
                        return None
                gain = score(parts[1])
                loss = score(parts[2])
                maximum = max((abs(value) for value in (gain, loss) if value is not None), default=None)
                rows.append({"VARIANT_ID": variant_id, "GENE": parts[0], "PANGOLIN_GAIN": gain, "PANGOLIN_LOSS": loss, "PANGOLIN_SCORE": maximum})
    frame = pl.DataFrame(rows) if rows else pl.DataFrame(schema={"VARIANT_ID": pl.String, "GENE": pl.String, "PANGOLIN_SCORE": pl.Float64})
    summary = frame.group_by("VARIANT_ID").agg(pl.max("PANGOLIN_SCORE"), pl.col("GENE").drop_nulls().unique().sort()) if rows else frame
    return frame, summary


def integrate_splicing(spliceai: pl.DataFrame, pangolin: pl.DataFrame) -> pl.DataFrame:
    return spliceai.join(pangolin, on="VARIANT_ID", how="full", coalesce=True)
