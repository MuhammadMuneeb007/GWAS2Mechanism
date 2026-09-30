"""Prioritised-variant VEP adapter and consequence classification."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver

_SPLICE = {"splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"}
_CODING = {"stop_gained", "stop_lost", "start_lost", "frameshift_variant", "missense_variant", "inframe_insertion", "inframe_deletion", "synonymous_variant", "coding_sequence_variant"}
_UTR = {"5_prime_UTR_variant", "3_prime_UTR_variant"}
_REGULATORY = {"regulatory_region_variant", "TF_binding_site_variant", "mature_miRNA_variant"}
_NONCODING = {"non_coding_transcript_exon_variant", "non_coding_transcript_variant", "NMD_transcript_variant"}


def functional_category(consequence: str | None) -> str:
    terms = set((consequence or "").split("&"))
    if terms & _SPLICE:
        return "SPLICING"
    if terms & _CODING:
        return "CODING"
    if terms & _UTR:
        return "UTR"
    if terms & _REGULATORY:
        return "REGULATORY"
    if terms & _NONCODING:
        return "NONCODING_RNA"
    if "intron_variant" in terms:
        return "INTRONIC"
    if "upstream_gene_variant" in terms:
        return "UPSTREAM"
    if "downstream_gene_variant" in terms:
        return "DOWNSTREAM"
    if "intergenic_variant" in terms:
        return "INTERGENIC"
    return "OTHER"


def select_prioritised_variants(
    population_results: list[pl.DataFrame],
    cross_ancestry: pl.DataFrame | None = None,
    *,
    min_pip: float = 0.01,
) -> pl.DataFrame:
    """Return the deduplicated union requested by fast mode."""
    selected: list[pl.DataFrame] = []
    for frame in population_results:
        if frame.is_empty():
            continue
        cs = pl.col("CREDIBLE_SETS").list.len() > 0 if "CREDIBLE_SETS" in frame.columns else pl.lit(False)
        top = pl.col("PIP") == pl.col("PIP").max().over("LOCUS_ID")
        selected.append(frame.filter((pl.col("PIP") >= min_pip) | cs | top))
    if cross_ancestry is not None and not cross_ancestry.is_empty():
        pip_col = "CROSS_ANCESTRY_PIP"
        cs = pl.col("CREDIBLE_SETS").list.len() > 0 if "CREDIBLE_SETS" in cross_ancestry.columns else pl.lit(False)
        selected.append(cross_ancestry.filter((pl.col(pip_col) >= min_pip) | cs))
    if not selected:
        return pl.DataFrame(schema={"VARIANT_ID": pl.String, "CHR": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String})
    return pl.concat(selected, how="diagonal_relaxed").unique("VARIANT_ID", keep="first")


def write_vep_vcf(variants: pl.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = variants.select("CHR", "POS", "VARIANT_ID", "REF", "ALT").unique("VARIANT_ID")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, variant_id, ref, alt in rows.iter_rows():
            handle.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\tPASS\t.\n")
    return path


def run_vep(input_vcf: Path, output_tsv: Path, tools: ToolResolver, cache_dir: Path, extra_args: list[str] | None = None) -> Path:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            *tools.require("vep"), "--input_file", str(input_vcf), "--output_file", str(output_tsv),
            "--format", "vcf", "--tab", "--cache", "--offline", "--dir_cache", str(cache_dir),
            "--assembly", "GRCh38", "--everything", "--force_overwrite", *(extra_args or []),
        ],
        log_file=output_tsv.with_suffix(".log"),
    )
    return output_tsv


def parse_vep_tsv(path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    frame = pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=10000, null_values=["-", "NA", ""])
    rename = {"#Uploaded_variation": "VARIANT_ID", "Consequence": "CONSEQUENCE", "SYMBOL": "GENE", "Gene": "GENE_ID"}
    frame = frame.rename({key: value for key, value in rename.items() if key in frame.columns})
    frame = frame.with_columns(
        pl.col("CONSEQUENCE").map_elements(functional_category, return_dtype=pl.String).alias("FUNCTIONAL_CATEGORY")
    )
    ranks = {"SPLICING": 0, "CODING": 1, "UTR": 2, "REGULATORY": 3, "NONCODING_RNA": 4, "INTRONIC": 5, "UPSTREAM": 6, "DOWNSTREAM": 7, "INTERGENIC": 8, "OTHER": 9}
    best = frame.with_columns(pl.col("FUNCTIONAL_CATEGORY").replace_strict(ranks, default=99).alias("_rank")).sort(["VARIANT_ID", "_rank"]).unique("VARIANT_ID", keep="first").drop("_rank")
    return frame, best
