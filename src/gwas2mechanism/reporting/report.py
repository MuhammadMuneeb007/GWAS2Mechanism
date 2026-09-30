"""Self-contained Markdown/HTML scientific report."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from gwas2mechanism.utils.io import atomic_write_text, write_tsv


def _load(path: Path) -> pl.DataFrame:
    return pl.read_parquet(path) if path.exists() else pl.DataFrame()


def _markdown_table(frame: pl.DataFrame, n: int = 20) -> str:
    if frame.is_empty():
        return "_No results available._"
    sample = frame.head(n)
    columns = sample.columns
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in sample.iter_rows():
        lines.append("| " + " | ".join(str(value).replace("|", "\\|") if value is not None else "" for value in row) + " |")
    return "\n".join(lines)


def generate_report(run_dir: Path, top_n: int = 50) -> tuple[Path, Path, Path]:
    manifest_path = run_dir / "manifest.yaml"
    manifest: dict[str, Any] = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    evidence = _load(run_dir / "evidence" / "evidence.parquet")
    comparison = _load(run_dir / "evidence" / "cross_population.parquet")
    selections = _load(run_dir / "discovery" / "study_selection.parquet")
    qtl = _load(run_dir / "qtl" / "all_tissues" / "variant_gene_tissue.parquet")
    failures = _load(run_dir / "finemapping" / "failed_loci.parquet")
    title = f"GWAS2Mechanism: {manifest.get('canonical_phenotype', manifest.get('phenotype', 'unknown phenotype'))}"
    sections = [
        f"# {title}",
        "",
        "> Results prioritize candidate causal variants and candidate molecular mechanisms; they do not establish experimental causality.",
        "",
        "## Selected GWAS studies",
        _markdown_table(selections, top_n),
        "",
        "## Population-specific and cross-ancestry evidence",
        _markdown_table(comparison, top_n),
        "",
        "## Top mechanistic candidates",
        _markdown_table(evidence.sort("PRIORITY_SCORE", descending=True) if "PRIORITY_SCORE" in evidence.columns else evidence, top_n),
        "",
        "## GTEx all-tissue evidence",
        f"All available tissues were eligible. Matched rows: {qtl.height}; tissues represented: {qtl['TISSUE'].n_unique() if 'TISSUE' in qtl.columns else 0}.",
        _markdown_table(qtl, top_n),
        "",
        "## Failed or skipped analyses",
        _markdown_table(failures, top_n),
        "",
        "## Scientific limitations",
        "1. Aggregated summary statistics cannot be retrospectively split by ancestry.",
        "2. LD must match cohort ancestry; 1000 Genomes may remain an imperfect proxy.",
        "3. Smaller studies can lack power; absence of evidence is not biological absence.",
        "4. PIP-based shared-variant screening is not formal colocalisation.",
        "5. VEP, SpliceAI and Pangolin outputs are predictions, not validation.",
        "6. Gene overlap and tissue QTL evidence do not prove mediation or causal tissue.",
        "7. Cross-population heterogeneity must be interpreted alongside meta-analysis.",
        "",
        "## Software and resources",
        "```yaml",
        yaml.safe_dump({"software": manifest.get("software", {}), "resources": manifest.get("resources", {})}, sort_keys=False).rstrip(),
        "```",
    ]
    markdown = "\n".join(sections) + "\n"
    report_dir = run_dir / "report"
    md_path = atomic_write_text(report_dir / "report.md", markdown)
    body = "<br>".join(html.escape(markdown).splitlines())
    html_path = atomic_write_text(report_dir / "report.html", f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title><style>body{{font:15px system-ui;max-width:1200px;margin:2rem auto;line-height:1.45}}code{{white-space:pre-wrap}}</style></head><body>{body}</body></html>\n")
    summary = pl.DataFrame([{"phenotype": manifest.get("phenotype"), "run_id": manifest.get("run_id"), "evidence_rows": evidence.height, "qtl_tissues": qtl["TISSUE"].n_unique() if "TISSUE" in qtl.columns else 0, "failed_loci": failures.height}])
    summary_path = write_tsv(summary, report_dir / "summary.tsv")
    return html_path, md_path, summary_path
