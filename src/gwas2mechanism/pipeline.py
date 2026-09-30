"""Run orchestration helpers, including the download-free synthetic smoke test."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from gwas2mechanism.annotation.vep import functional_category, select_prioritised_variants
from gwas2mechanism.clumping.loci import construct_loci, greedy_ld_clump, union_loci
from gwas2mechanism.coloc.screen import screen
from gwas2mechanism.config import Config, save_config
from gwas2mechanism.evidence.integrate import (
    add_priority_components,
    cross_population_comparison,
    integrate_evidence,
)
from gwas2mechanism.finemap.susie import run_numpy_smoke_susie
from gwas2mechanism.paths import RunLayout, create_run
from gwas2mechanism.provenance import RunManifest
from gwas2mechanism.reporting.report import generate_report
from gwas2mechanism.utils.io import atomic_write_json, write_parquet, write_tsv


def _synthetic_sumstats(population: str, shift: float) -> pl.DataFrame:
    pos = np.arange(100_000, 100_000 + 20 * 1000, 1000, dtype=np.int64)
    beta = np.linspace(0.05, 0.2, pos.size) + shift
    se = np.full(pos.size, 0.025)
    beta[9] = 0.36 + shift
    z = beta / se
    p = 2 * stats.norm.sf(np.abs(z))
    ref = ["A" if i % 2 == 0 else "C" for i in range(pos.size)]
    alt = ["G" if i % 2 == 0 else "T" for i in range(pos.size)]
    variant_id = [f"1:{p0}:{min(r, a)}:{max(r, a)}" for p0, r, a in zip(pos, ref, alt, strict=True)]
    return pl.DataFrame(
        {
            "STUDY_ID": [f"SYNTH_{population}"] * pos.size,
            "TRAIT": ["synthetic phenotype"] * pos.size,
            "POPULATION": [population] * pos.size,
            "CHR": ["1"] * pos.size,
            "POS": pos,
            "RSID": [f"rs{i + 1}" for i in range(pos.size)],
            "REF": ref,
            "ALT": alt,
            "EA": alt,
            "NEA": ref,
            "BETA": beta,
            "SE": se,
            "OR": np.exp(beta),
            "Z": z,
            "P": p,
            "NEG_LOG10_P": -np.log10(p),
            "EAF": np.linspace(0.08, 0.45, pos.size),
            "N": [50_000 if population == "EUR" else 20_000] * pos.size,
            "N_CASES": [10_000 if population == "EUR" else 4_000] * pos.size,
            "N_CONTROLS": [40_000 if population == "EUR" else 16_000] * pos.size,
            "N_EFF": [32_000 if population == "EUR" else 12_800] * pos.size,
            "BUILD": ["GRCh38"] * pos.size,
            "SOURCE_FILE": ["synthetic_fixture"] * pos.size,
            "VARIANT_ID": variant_id,
        }
    )


def _synthetic_ld(n: int, decay: float) -> np.ndarray:
    index = np.arange(n)
    return np.power(decay, np.abs(index[:, None] - index[None, :]), dtype=np.float64)


def run_synthetic_smoke(config: Config, phenotype: str, run_id: str = "synthetic-smoke") -> RunLayout:
    """Execute every internal layer without network access or external tools."""
    layout = create_run(Path(config.run_root), phenotype, run_id)
    save_config(config, layout.config)
    manifest = RunManifest(layout.manifest)
    manifest.initialise(phenotype=phenotype, run_id=run_id, command_line=sys.argv, parameters=config.dump())
    manifest.update("canonical_phenotype", phenotype)
    manifest.update("mode", "synthetic_smoke")

    selection = pl.DataFrame(
        [
            {"study_accession": "SYNTH_EUR", "reported_trait": phenotype, "population": "EUR", "reported_ancestry": "European", "n_total": 50_000, "n_cases": 10_000, "n_controls": 40_000, "summary_stats_available": True, "harmonised_available": True, "build": "GRCh38", "selected": True, "selection_reason": "synthetic fixture"},
            {"study_accession": "SYNTH_AFR", "reported_trait": phenotype, "population": "AFR", "reported_ancestry": "African", "n_total": 20_000, "n_cases": 4_000, "n_controls": 16_000, "summary_stats_available": True, "harmonised_available": True, "build": "GRCh38", "selected": True, "selection_reason": "synthetic fixture"},
        ]
    )
    write_parquet(selection, layout.study_selection)
    write_parquet(selection, layout.candidate_studies)
    write_parquet(pl.DataFrame(schema={"study_accession": pl.String, "url": pl.String}), layout.candidate_files)
    write_tsv(selection, layout.selection_report_tsv)
    atomic_write_json(layout.phenotype_resolution, {"phenotype_input": phenotype, "phenotype_canonical": phenotype, "phenotype_slug": phenotype.replace(" ", "_"), "source": "synthetic"})

    sumstats: dict[str, pl.DataFrame] = {}
    finemaps: list[pl.DataFrame] = []
    loci_frames: list[pl.DataFrame] = []
    cs_frames: list[pl.DataFrame] = []
    summaries: list[pl.DataFrame] = []
    diagnostics: list[pl.DataFrame] = []
    for population, shift, decay in (("EUR", 0.0, 0.85), ("AFR", -0.015, 0.65)):
        frame = _synthetic_sumstats(population, shift)
        sumstats[population] = frame
        write_parquet(frame, layout.sumstats(population))
        write_parquet(frame.filter(pl.col("P") < config.gwas.p_threshold), layout.significant(population))
        atomic_write_json(layout.gwas_qc(population), {"population": population, "n_variants": frame.height, "build": "GRCh38", "fixture": True})
        ld = _synthetic_ld(frame.height, decay)
        clumps = greedy_ld_clump(frame, ld, p1=config.clumping.p1, p2=config.clumping.p2, r2=config.clumping.r2, kb=config.clumping.kb)
        loci = construct_loci(clumps, population=population, window_bp=50_000)
        write_parquet(loci, layout.loci(population))
        loci_frames.append(loci)
        locus_id = loci["LOCUS_ID"][0]
        result = run_numpy_smoke_susie(frame, ld, locus_id=locus_id, population=population, coverage=config.finemapping.credible_set_coverage, max_effects=2)
        write_parquet(result.variants, layout.finemap_dir(population) / "finemapped_variants.parquet")
        write_parquet(result.credible_sets, layout.finemap_dir(population) / "credible_sets.parquet")
        write_parquet(result.summary, layout.finemap_dir(population) / "finemap_summary.parquet")
        write_parquet(result.diagnostics, layout.finemap_dir(population) / "finemap_diagnostics.parquet")
        finemaps.append(result.variants)
        cs_frames.append(result.credible_sets)
        summaries.append(result.summary)
        diagnostics.append(result.diagnostics)
    write_parquet(union_loci(loci_frames), layout.combined_loci)
    write_parquet(pl.concat(summaries), layout.root / "finemapping" / "finemap_summary.parquet")
    write_parquet(pl.concat(diagnostics), layout.root / "finemapping" / "finemap_diagnostics.parquet")
    write_parquet(pl.DataFrame(schema={"LOCUS_ID": pl.String, "REASON": pl.String}), layout.root / "finemapping" / "failed_loci.parquet")

    all_finemap = pl.concat(finemaps, how="diagonal_relaxed")
    cross = all_finemap.group_by("VARIANT_ID").agg((1.0 - (1.0 - pl.col("PIP")).product()).alias("CROSS_ANCESTRY_PIP"))
    write_parquet(cross, layout.multiancestry_dir / "variants.parquet")
    write_parquet(pl.concat(cs_frames, how="diagonal_relaxed"), layout.multiancestry_dir / "credible_sets.parquet")
    write_parquet(pl.DataFrame([{"METHOD": "synthetic_consensus", "STATUS": "fixture_only", "N_POPULATIONS": 2}]), layout.multiancestry_dir / "locus_summary.parquet")
    write_parquet(pl.DataFrame([{"NOTE": "Synthetic smoke test; production uses MultiSuSiE with separate LD matrices"}]), layout.multiancestry_dir / "diagnostics.parquet")

    prioritised = select_prioritised_variants(finemaps, cross, min_pip=config.annotation.min_pip)
    consequences = ["missense_variant" if i % 7 == 0 else "intron_variant" for i in range(prioritised.height)]
    annotation = prioritised.select("VARIANT_ID").with_columns(
        pl.Series("CONSEQUENCE", consequences),
        pl.Series("FUNCTIONAL_CATEGORY", [functional_category(value) for value in consequences]),
        pl.Series("GENE", [f"GENE{(i % 4) + 1}" for i in range(prioritised.height)]),
    )
    write_parquet(annotation, layout.annotation / "vep_all_consequences.parquet")
    write_parquet(annotation, layout.annotation / "vep_best_consequence.parquet")
    write_parquet(annotation, layout.annotation / "variant_summary.parquet")

    qtl_rows = []
    tissues = ["Artery_Aorta", "Brain_Cortex", "Liver"]
    for i, variant_id in enumerate(prioritised["VARIANT_ID"].to_list()[:9]):
        for qtl_type in ("eQTL", "sQTL"):
            qtl_rows.append({"VARIANT_ID": variant_id, "TISSUE": tissues[i % len(tissues)], "QTL_TYPE": qtl_type, "GENE_ID": f"ENSG{i % 4 + 1}", "GENE": f"GENE{i % 4 + 1}", "QTL_PHENOTYPE": f"QTL{i}", "QTL_PIP": float(0.1 + 0.08 * (i % 5))})
    qtl = pl.DataFrame(qtl_rows)
    write_parquet(qtl.filter(pl.col("QTL_TYPE") == "eQTL"), layout.qtl / "eqtl_variant_matches.parquet")
    write_parquet(qtl.filter(pl.col("QTL_TYPE") == "sQTL"), layout.qtl / "sqtl_variant_matches.parquet")
    write_parquet(qtl, layout.qtl / "variant_gene_tissue.parquet")
    write_parquet(qtl.group_by("GENE", "TISSUE", "QTL_TYPE").agg(pl.max("QTL_PIP"), pl.len().alias("N_VARIANTS")), layout.qtl / "gene_tissue_summary.parquet")
    write_parquet(qtl.group_by("TISSUE", "QTL_TYPE").agg(pl.len().alias("N_VARIANTS")), layout.qtl / "locus_summary.parquet")

    splice = prioritised.select("VARIANT_ID").with_columns(pl.Series("SPLICEAI_MAX", np.linspace(0.01, 0.9, prioritised.height)), pl.Series("PANGOLIN_SCORE", np.linspace(0.02, 0.7, prioritised.height)))
    write_parquet(splice.select("VARIANT_ID", "SPLICEAI_MAX"), layout.splicing / "spliceai" / "all_predictions.parquet")
    write_parquet(splice.select("VARIANT_ID", "SPLICEAI_MAX"), layout.splicing / "spliceai" / "variant_summary.parquet")
    write_parquet(pl.DataFrame(schema={"VARIANT_ID": pl.String, "REASON": pl.String}), layout.splicing / "spliceai" / "unscored.parquet")
    write_parquet(splice.select("VARIANT_ID", "PANGOLIN_SCORE"), layout.splicing / "pangolin" / "all_predictions.parquet")
    write_parquet(splice.select("VARIANT_ID", "PANGOLIN_SCORE"), layout.splicing / "pangolin" / "variant_summary.parquet")
    write_parquet(splice, layout.splicing / "integrated_splicing.parquet")

    variant_coloc, gene_tissue = screen(all_finemap, qtl, cross)
    write_parquet(variant_coloc, layout.coloc / "variant_colocalization.parquet")
    write_parquet(gene_tissue, layout.coloc / "gene_tissue_colocalization.parquet")
    gene_summary = gene_tissue.group_by("GENE").agg(pl.max("LCLPP"), pl.len().alias("N_PAIRS"))
    write_parquet(gene_summary, layout.coloc / "gene_summary.parquet")
    write_parquet(gene_tissue.group_by("LOCUS_ID").agg(pl.max("LCLPP"), pl.len().alias("N_PAIRS")), layout.coloc / "locus_summary.parquet")
    write_parquet(gene_tissue.sort("LCLPP", descending=True).head(100), layout.coloc / "top_candidates.parquet")

    evidence = integrate_evidence(all_finemap, phenotype=phenotype, annotation=annotation, qtl=qtl, coloc=variant_coloc, splicing=splice, cross_ancestry=cross)
    comparison = cross_population_comparison(all_finemap, config.evidence.support_pip)
    evidence = add_priority_components(evidence.join(comparison, on="VARIANT_ID", how="left"), config.evidence.priority_score.weights)
    write_parquet(evidence, layout.evidence / "evidence.parquet")
    write_parquet(comparison, layout.evidence / "cross_population.parquet")
    write_tsv(evidence, layout.evidence / "evidence.tsv")
    manifest.update("study_accessions", ["SYNTH_EUR", "SYNTH_AFR"])
    manifest.update("population_assignments", {"SYNTH_EUR": "EUR", "SYNTH_AFR": "AFR"})
    manifest.update("sample_sizes", {"EUR": 50_000, "AFR": 20_000})
    manifest.update("gtex_version", "synthetic_fixture")
    manifest.update("synthetic", True)
    generate_report(layout.root, config.report.top_n)
    return layout
