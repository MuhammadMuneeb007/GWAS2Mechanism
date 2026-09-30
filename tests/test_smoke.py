from __future__ import annotations

from pathlib import Path

import polars as pl
from typer.testing import CliRunner

from gwas2mechanism.cli import app
from gwas2mechanism.config import load_config
from gwas2mechanism.pipeline import run_synthetic_smoke


def test_complete_synthetic_smoke(tmp_path: Path) -> None:
    cfg = load_config(overrides={"run_root": str(tmp_path / "runs")})
    layout = run_synthetic_smoke(cfg, "asthma", run_id="test")
    assert (layout.report / "report.html").exists()
    assert (layout.evidence / "evidence.parquet").exists()
    qtl = pl.read_parquet(layout.qtl / "variant_gene_tissue.parquet")
    assert set(qtl["TISSUE"]) == {"Artery_Aorta", "Brain_Cortex", "Liver"}


def test_cli_help_and_two_phenotypes_are_accepted(tmp_path: Path) -> None:
    runner = CliRunner()
    assert runner.invoke(app, ["--help"]).exit_code == 0
    for phenotype in ("migraine", "asthma"):
        result = runner.invoke(app, ["run", "--phenotype", phenotype, "--synthetic", "--set", f"run_root={tmp_path / phenotype}"])
        assert result.exit_code == 0, result.output
