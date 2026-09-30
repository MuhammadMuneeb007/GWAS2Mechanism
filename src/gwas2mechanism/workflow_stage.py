"""Implementation behind modular Snakemake stage rules.

Every marker is written only after expected scientific outputs exist. Expensive
work is therefore safely restartable and failed loci are recorded rather than
silently converted into successful results.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import polars as pl

from gwas2mechanism.annotation.vep import (
    parse_vep_tsv,
    run_vep,
    select_prioritised_variants,
    write_vep_vcf,
)
from gwas2mechanism.clumping.loci import construct_loci, run_plink_clump, union_loci
from gwas2mechanism.coloc.screen import screen
from gwas2mechanism.config import Config, load_config
from gwas2mechanism.discovery.service import discover
from gwas2mechanism.evidence.integrate import (
    add_priority_components,
    cross_population_comparison,
    integrate_evidence,
)
from gwas2mechanism.finemap.ld import LDProvider
from gwas2mechanism.finemap.multisusie import MultiSuSiEInputs, run_multisusie
from gwas2mechanism.finemap.region import extract_region, match_to_reference
from gwas2mechanism.finemap.susie import FineMapResult, run_susier
from gwas2mechanism.meta.fixed_effect import run_meta
from gwas2mechanism.paths import RunLayout
from gwas2mechanism.qtl.gtex import integrate_local_archive
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.reporting.report import generate_report
from gwas2mechanism.setup_resources import setup_resources
from gwas2mechanism.splicing.pangolin import integrate_splicing, parse_pangolin_table, run_pangolin
from gwas2mechanism.splicing.spliceai import parse_spliceai_vcf, run_spliceai
from gwas2mechanism.sumstats.download import download_analysis
from gwas2mechanism.sumstats.harmonize import harmonize_analysis
from gwas2mechanism.utils.io import (
    atomic_write_json,
    atomic_write_text,
    concat_parquets,
    read_json,
    write_parquet,
    write_tsv,
)
from gwas2mechanism.utils.runtime import resolve_threads
from gwas2mechanism.utils.tools import ToolResolver
from gwas2mechanism.utils.variants import split_variant_id

log = logging.getLogger(__name__)


class StageRunner:
    def __init__(self, run_dir: Path, phenotype: str):
        self.layout = RunLayout(run_dir.resolve())
        self.phenotype = phenotype
        self.config: Config = load_config([self.layout.config])
        self.tools = ToolResolver(self.config.tools)
        self.threads = resolve_threads(self.config.performance.threads)

    @property
    def analyses(self) -> list[dict[str, Any]]:
        return read_json(self.layout.selected_analyses) if self.layout.selected_analyses.exists() else []

    def mark(self, stage: str) -> None:
        atomic_write_text(self.layout.stage_marker(stage), "OK\n")

    def discover(self) -> None:
        discover(self.phenotype, self.config, self.layout, self.config.populations.allowed if self.config.populations.mode == "explicit" else None)

    def download(self) -> None:
        if not self.analyses:
            raise RuntimeError("No selected studies; inspect discovery/selection_report.tsv")
        records = []
        for analysis in self.analyses:
            record = download_analysis(analysis, self.layout, self.config, self.tools)
            atomic_write_json(self.layout.gwas_dir(analysis["analysis_label"]) / "download.json", record)
            records.append(record)
        atomic_write_json(self.layout.root / "gwas" / "downloads.json", records)

    def harmonize(self) -> None:
        usable: dict[str, Path] = {}
        types: dict[str, str | None] = {}
        fasta = self.config.resources_path / "genome" / "GRCh38" / "GRCh38.primary_assembly.fa"
        if not fasta.exists():
            fasta = None
        for analysis in self.analyses:
            label = analysis["analysis_label"]
            download = read_json(self.layout.gwas_dir(label) / "download.json")
            harmonize_analysis(analysis, download, self.layout.sumstats(label), self.layout.significant(label), self.layout.gwas_qc(label), self.config, self.phenotype, self.threads, fasta)
            if label in self.config.populations.allowed:
                usable[label] = self.layout.sumstats(label)
                types[label] = analysis.get("trait_type")
        if self.config.combined.fixed_effect_meta and len(usable) >= self.config.combined.min_populations_for_meta:
            run_meta(usable, types, self.layout.meta_sumstats, self.phenotype, self.config.harmonization.drop_palindromic_ambiguous_maf)

    def references(self) -> None:
        populations = sorted({a["analysis_label"] for a in self.analyses if a["analysis_label"] in self.config.populations.allowed})
        setup_resources(self.config, reference=True, gtex=False, vep=False, splice=False, populations=populations, full=False, dry_run=False)

    @staticmethod
    def _read_plink_leads(path: Path) -> pl.DataFrame:
        frame = pl.read_csv(path, separator="\t", infer_schema_length=10000, null_values=["NA", "."])
        rename = {"#CHROM": "CHR", "POS": "LEAD_POS", "ID": "LEAD_VARIANT", "P": "LEAD_P"}
        frame = frame.rename({key: value for key, value in rename.items() if key in frame.columns})
        needed = {"CHR", "LEAD_POS", "LEAD_VARIANT", "LEAD_P"}
        if not needed <= set(frame.columns):
            raise ValueError(f"Unexpected PLINK2 clump columns: {frame.columns}")
        return frame.select("CHR", "LEAD_POS", "LEAD_VARIANT", "LEAD_P").with_columns(pl.col("CHR").cast(pl.String))

    def clump(self) -> None:
        panel = ReferencePanel(self.config.reference, self.config.resources_path, self.tools)
        frames = []
        for analysis in self.analyses:
            population = analysis["analysis_label"]
            if population not in self.config.populations.allowed or not self.layout.significant(population).exists():
                continue
            significant = pl.read_parquet(self.layout.significant(population))
            leads = []
            for chrom in significant["CHR"].unique().to_list():
                subset = significant.filter(pl.col("CHR") == chrom).select("VARIANT_ID", "P")
                if subset.is_empty():
                    continue
                input_path = self.layout.loci_dir(population) / f"chr{chrom}.clump.tsv"
                input_path.parent.mkdir(parents=True, exist_ok=True)
                subset.write_csv(input_path, separator="\t")
                prefix = self.layout.loci_dir(population) / f"chr{chrom}"
                result = run_plink_clump(input_path, panel.prefix(population, str(chrom)), prefix, self.tools, p1=self.config.clumping.p1, p2=self.config.clumping.p2, r2=self.config.clumping.r2, kb=self.config.clumping.kb, threads=self.threads)
                leads.append(self._read_plink_leads(result))
            lead_frame = pl.concat(leads, how="diagonal_relaxed") if leads else pl.DataFrame()
            loci = construct_loci(lead_frame, population=population, window_bp=self.config.loci.window_bp)
            write_parquet(loci, self.layout.loci(population))
            frames.append(loci)
        if self.layout.meta_sumstats.exists():
            # Meta-analysis signals contribute to the union. They are not clumped
            # with a fake pooled LD panel: fixed windows around significant hits are used.
            meta = pl.read_parquet(self.layout.meta_sumstats).filter(pl.col("P") < self.config.gwas.p_threshold)
            if not meta.is_empty():
                meta_leads = meta.select("CHR", pl.col("POS").alias("LEAD_POS"), pl.col("VARIANT_ID").alias("LEAD_VARIANT"), pl.col("P").alias("LEAD_P"))
                frames.append(construct_loci(meta_leads, population="COMBINED_META", window_bp=self.config.loci.window_bp))
        write_parquet(union_loci(frames, self.config.loci.combined_merge_bp), self.layout.combined_loci)

    def _finemap_one(self, population: str, locus: dict[str, Any], panel: ReferencePanel, ld_provider: LDProvider) -> FineMapResult:
        regional = extract_region(self.layout.sumstats(population), locus["CHR"], locus["START"], locus["END"])
        reference = panel.region_variants(population, locus["CHR"], locus["START"], locus["END"])
        matched = match_to_reference(regional, reference)
        if matched.height < self.config.finemapping.min_variants:
            raise ValueError(f"only {matched.height} variants match the population reference")
        workdir = self.layout.finemap_locus_dir(population, locus["LOCUS_ID"].replace(":", "_"))
        ld_result = ld_provider.compute(population, locus["CHR"], matched["VARIANT_ID"].to_list(), workdir, self.config.finemapping.ld_ridge, self.config.finemapping.max_ld_variants)
        order = pl.DataFrame({"VARIANT_ID": ld_result.ids, "_ORDER": range(len(ld_result.ids))})
        matched = order.join(matched, on="VARIANT_ID", how="inner").sort("_ORDER").drop("_ORDER")
        if matched.height != len(ld_result.ids):
            raise ValueError("GWAS/LD variant order mismatch")
        variants = matched.drop("Z", strict=False).rename({"Z_ALT": "Z"})
        n_col = "N_EFF" if variants["N_EFF"].drop_nulls().len() else "N"
        sample_size = float(variants[n_col].drop_nulls().median())
        return run_susier(variants, ld_result.R, sample_size=sample_size, locus_id=locus["LOCUS_ID"], population=population, workdir=workdir, tools=self.tools, coverage=self.config.finemapping.credible_set_coverage, max_effects=self.config.finemapping.default_L)

    def finemap(self) -> None:
        panel = ReferencePanel(self.config.reference, self.config.resources_path, self.tools)
        provider = LDProvider(panel, self.config.reference.ld_engine, self.config.reference.maf, self.config.reference.geno, self.tools.command("plink2"), self.threads)
        all_failed = []
        for population in self.config.populations.allowed:
            if not self.layout.loci(population).exists():
                continue
            results = []
            for locus in pl.read_parquet(self.layout.loci(population)).iter_rows(named=True):
                try:
                    results.append(self._finemap_one(population, locus, panel, provider))
                except Exception as exc:  # record per-locus failures; continue independent jobs
                    all_failed.append({"LOCUS_ID": locus["LOCUS_ID"], "POPULATION": population, "REASON": str(exc)})
            out = self.layout.finemap_dir(population)
            for name, attr in (("finemapped_variants.parquet", "variants"), ("credible_sets.parquet", "credible_sets"), ("finemap_summary.parquet", "summary"), ("finemap_diagnostics.parquet", "diagnostics")):
                write_parquet(pl.concat([getattr(result, attr) for result in results], how="diagonal_relaxed") if results else pl.DataFrame(), out / name)
        write_parquet(pl.DataFrame(all_failed) if all_failed else pl.DataFrame(schema={"LOCUS_ID": pl.String, "POPULATION": pl.String, "REASON": pl.String}), self.layout.root / "finemapping" / "failed_loci.parquet")

    def multiancestry(self) -> None:
        populations = [population for population in self.config.populations.allowed if self.layout.sumstats(population).exists()]
        variants_out: list[pl.DataFrame] = []
        summaries: list[pl.DataFrame] = []
        diagnostics: list[dict[str, Any]] = []
        panel = ReferencePanel(self.config.reference, self.config.resources_path, self.tools)
        provider = LDProvider(panel, self.config.reference.ld_engine, self.config.reference.maf, self.config.reference.geno, self.tools.command("plink2"), self.threads)
        if len(populations) >= 2 and self.layout.combined_loci.exists():
            for locus in pl.read_parquet(self.layout.combined_loci).iter_rows(named=True):
                try:
                    matched_by_pop: dict[str, pl.DataFrame] = {}
                    for population in populations:
                        regional = extract_region(self.layout.sumstats(population), locus["CHR"], locus["START"], locus["END"])
                        reference = panel.region_variants(population, locus["CHR"], locus["START"], locus["END"])
                        matched_by_pop[population] = match_to_reference(regional, reference)
                    common = set.intersection(*(set(frame["VARIANT_ID"].to_list()) for frame in matched_by_pop.values()))
                    if len(common) < self.config.finemapping.min_variants:
                        raise ValueError(f"only {len(common)} variants are shared across population references")
                    ordered_ids = sorted(common)
                    z_list = []
                    ld_list = []
                    n_list = []
                    for population in populations:
                        workdir = self.layout.multiancestry_dir / "loci" / locus["LOCUS_ID"].replace(":", "_") / population
                        ld_result = provider.compute(population, locus["CHR"], ordered_ids, workdir, self.config.finemapping.ld_ridge, self.config.finemapping.max_ld_variants)
                        if set(ld_result.ids) != common:
                            raise ValueError(f"{population} LD filtering changed the shared variant set")
                        ordered = pl.DataFrame({"VARIANT_ID": ld_result.ids, "_ORDER": range(len(ld_result.ids))}).join(matched_by_pop[population], on="VARIANT_ID", how="inner").sort("_ORDER")
                        z_list.append(ordered["Z_ALT"].to_numpy())
                        ld_list.append(ld_result.R)
                        n_source = ordered["N_EFF"].drop_nulls() if "N_EFF" in ordered.columns and ordered["N_EFF"].drop_nulls().len() else ordered["N"].drop_nulls()
                        n_list.append(float(n_source.median()))
                    inputs = MultiSuSiEInputs(populations, z_list, ld_list, n_list, ld_result.ids)
                    variant_result, summary = run_multisusie(inputs, coverage=self.config.finemapping.credible_set_coverage, rho=self.config.finemapping.multisusie_rho)
                    variants_out.append(variant_result.with_columns(pl.lit(locus["LOCUS_ID"]).alias("LOCUS_ID")))
                    summaries.append(summary.with_columns(pl.lit(locus["LOCUS_ID"]).alias("LOCUS_ID")))
                    diagnostics.append({"LOCUS_ID": locus["LOCUS_ID"], "STATUS": "OK", "POPULATIONS": ",".join(populations)})
                except Exception as exc:
                    diagnostics.append({"LOCUS_ID": locus["LOCUS_ID"], "STATUS": "FAILED", "REASON": str(exc)})
        else:
            diagnostics.append({"STATUS": "SKIPPED", "REASON": "fewer than two compatible population datasets"})
        write_parquet(pl.concat(variants_out, how="diagonal_relaxed") if variants_out else pl.DataFrame(), self.layout.multiancestry_dir / "variants.parquet")
        # MultiSuSiE API variants differ in credible-set representation; raw
        # credible sets are emitted by future adapters only when unambiguous.
        write_parquet(pl.DataFrame(), self.layout.multiancestry_dir / "credible_sets.parquet")
        write_parquet(pl.concat(summaries, how="diagonal_relaxed") if summaries else pl.DataFrame(), self.layout.multiancestry_dir / "locus_summary.parquet")
        write_parquet(pl.DataFrame(diagnostics), self.layout.multiancestry_dir / "diagnostics.parquet")

    def annotate(self) -> None:
        frames = [pl.read_parquet(path) for path in self.layout.root.glob("finemapping/*/finemapped_variants.parquet")]
        cross_path = self.layout.multiancestry_dir / "variants.parquet"
        cross = pl.read_parquet(cross_path) if cross_path.exists() else None
        prioritised = select_prioritised_variants(frames, cross, min_pip=self.config.annotation.min_pip)
        if prioritised.is_empty():
            raise RuntimeError("No prioritised variants to annotate")
        vcf = write_vep_vcf(prioritised, self.layout.annotation / "prioritised.vcf")
        raw = run_vep(vcf, self.layout.annotation / "vep.tsv", self.tools, self.config.resources_path / "vep", self.config.annotation.vep_extra_args)
        all_consequences, best = parse_vep_tsv(raw)
        write_parquet(all_consequences, self.layout.annotation / "vep_all_consequences.parquet")
        write_parquet(best, self.layout.annotation / "vep_best_consequence.parquet")
        write_parquet(prioritised.join(best, on="VARIANT_ID", how="left"), self.layout.annotation / "variant_summary.parquet")

    def qtl(self) -> None:
        variants = pl.read_parquet(self.layout.annotation / "variant_summary.parquet")
        if self.config.qtl.local_archive_dir:
            root = Path(self.config.qtl.local_archive_dir)
        else:
            gtex_root = self.config.resources_path / "gtex"
            if self.config.qtl.release == "latest":
                releases = sorted(gtex_root.glob("v*"), key=lambda path: int(path.name[1:]) if path.name[1:].isdigit() else -1)
                if not releases:
                    raise FileNotFoundError("No cached GTEx release; run `gwas2m setup --gtex`")
                root = releases[-1]
            else:
                root = gtex_root / self.config.qtl.release
        matches = integrate_local_archive(variants, root, tissues=self.config.qtl.tissues, qtl_types=tuple(value.replace("eqtl", "eQTL").replace("sqtl", "sQTL").replace("apaqtl", "apaQTL") for value in self.config.qtl.types))
        write_parquet(matches.filter(pl.col("QTL_TYPE") == "eQTL"), self.layout.qtl / "eqtl_variant_matches.parquet")
        write_parquet(matches.filter(pl.col("QTL_TYPE") == "sQTL"), self.layout.qtl / "sqtl_variant_matches.parquet")
        write_parquet(matches, self.layout.qtl / "variant_gene_tissue.parquet")
        write_parquet(matches.group_by("GENE", "TISSUE", "QTL_TYPE").agg(pl.max("QTL_PIP"), pl.len().alias("N_VARIANTS")), self.layout.qtl / "gene_tissue_summary.parquet")
        write_parquet(matches.group_by("TISSUE", "QTL_TYPE").agg(pl.len().alias("N_VARIANTS")), self.layout.qtl / "locus_summary.parquet")

    def splice(self) -> None:
        variants = pl.read_parquet(self.layout.annotation / "variant_summary.parquet")
        if not {"CHR", "POS", "REF", "ALT"} <= set(variants.columns):
            variants = split_variant_id(variants)
        fasta = self.config.resources_path / "genome" / "GRCh38" / "GRCh38.primary_assembly.fa.gz"
        spliceai_summary = pl.DataFrame(schema={"VARIANT_ID": pl.String, "SPLICEAI_MAX": pl.Float64})
        pangolin_summary = pl.DataFrame(schema={"VARIANT_ID": pl.String, "PANGOLIN_SCORE": pl.Float64})
        if self.config.splicing.spliceai:
            out = run_spliceai(variants, self.layout.splicing / "spliceai" / "predictions.vcf", fasta, self.config.splicing.spliceai_annotation, self.tools, self.config.splicing.spliceai_distance, self.config.splicing.spliceai_mask)
            all_predictions, spliceai_summary = parse_spliceai_vcf(out)
            write_parquet(all_predictions, self.layout.splicing / "spliceai" / "all_predictions.parquet")
            write_parquet(spliceai_summary, self.layout.splicing / "spliceai" / "variant_summary.parquet")
            write_parquet(variants.join(spliceai_summary, on="VARIANT_ID", how="anti"), self.layout.splicing / "spliceai" / "unscored.parquet")
        if self.config.splicing.pangolin:
            annotation_db = self.config.resources_path / "splicing" / "pangolin" / "gencode.annotation.db"
            out = run_pangolin(
                variants,
                self.layout.splicing / "pangolin" / "predictions.csv",
                fasta,
                annotation_db,
                self.tools,
                self.config.splicing.pangolin_distance,
                self.config.splicing.pangolin_mask,
            )
            all_predictions, pangolin_summary = parse_pangolin_table(out)
            write_parquet(all_predictions, self.layout.splicing / "pangolin" / "all_predictions.parquet")
            write_parquet(pangolin_summary, self.layout.splicing / "pangolin" / "variant_summary.parquet")
        write_parquet(integrate_splicing(spliceai_summary, pangolin_summary), self.layout.splicing / "integrated_splicing.parquet")

    def coloc(self) -> None:
        finemap = concat_parquets(self.layout.root.glob("finemapping/*/finemapped_variants.parquet"))
        qtl = pl.read_parquet(self.layout.qtl / "variant_gene_tissue.parquet")
        cross_path = self.layout.multiancestry_dir / "variants.parquet"
        cross = pl.read_parquet(cross_path) if cross_path.exists() else None
        variant, pairs = screen(finemap, qtl, cross)
        write_parquet(variant, self.layout.coloc / "variant_colocalization.parquet")
        write_parquet(pairs, self.layout.coloc / "gene_tissue_colocalization.parquet")
        gene = pairs.group_by("GENE").agg(pl.max("LCLPP"), pl.len().alias("N_PAIRS"))
        write_parquet(gene, self.layout.coloc / "gene_summary.parquet")
        write_parquet(pairs.group_by("LOCUS_ID").agg(pl.max("LCLPP"), pl.len().alias("N_PAIRS")), self.layout.coloc / "locus_summary.parquet")
        write_parquet(pairs.sort("LCLPP", descending=True).head(100), self.layout.coloc / "top_candidates.parquet")
        annotation = pl.read_parquet(self.layout.annotation / "vep_best_consequence.parquet")
        splicing = pl.read_parquet(self.layout.splicing / "integrated_splicing.parquet")
        evidence = integrate_evidence(finemap, phenotype=self.phenotype, annotation=annotation, qtl=qtl, coloc=variant, splicing=splicing, cross_ancestry=cross)
        comparison = cross_population_comparison(finemap, self.config.evidence.support_pip)
        evidence = add_priority_components(evidence.join(comparison, on="VARIANT_ID", how="left"), self.config.evidence.priority_score.weights)
        write_parquet(evidence, self.layout.evidence / "evidence.parquet")
        write_parquet(comparison, self.layout.evidence / "cross_population.parquet")
        write_tsv(evidence, self.layout.evidence / "evidence.tsv")

    def report(self) -> None:
        generate_report(self.layout.root, self.config.report.top_n)

    def run(self, stage: str) -> None:
        getattr(self, stage)()
        self.mark(stage)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["discover", "download", "harmonize", "references", "clump", "finemap", "multiancestry", "annotate", "qtl", "splice", "coloc", "report"])
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--phenotype", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    StageRunner(args.run_dir, args.phenotype).run(args.stage)


if __name__ == "__main__":
    main()
