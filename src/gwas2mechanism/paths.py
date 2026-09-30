"""Run directory layout - the single source of truth for output locations."""

from __future__ import annotations

import datetime as _dt
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from gwas2mechanism.constants import COMBINED_META


def slugify(text: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_")
    return slug or "phenotype"


def new_run_id() -> str:
    return _dt.datetime.now().strftime("%Y%m%d-%H%M%S")


@dataclass(frozen=True)
class RunLayout:
    root: Path

    # -- top level -------------------------------------------------------
    @property
    def manifest(self) -> Path:
        return self.root / "manifest.yaml"

    @property
    def config(self) -> Path:
        return self.root / "config.yaml"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def stages(self) -> Path:
        return self.root / ".stages"

    def stage_marker(self, stage: str, key: str | None = None) -> Path:
        name = f"{stage}.{key}.done" if key else f"{stage}.done"
        return self.stages / name

    # -- discovery -------------------------------------------------------
    @property
    def discovery(self) -> Path:
        return self.root / "discovery"

    @property
    def phenotype_resolution(self) -> Path:
        return self.discovery / "phenotype_resolution.json"

    @property
    def candidate_studies(self) -> Path:
        return self.discovery / "candidate_studies.parquet"

    @property
    def candidate_files(self) -> Path:
        return self.discovery / "candidate_files.parquet"

    @property
    def study_selection(self) -> Path:
        return self.discovery / "study_selection.parquet"

    @property
    def selection_report_tsv(self) -> Path:
        return self.discovery / "selection_report.tsv"

    @property
    def selection_report_html(self) -> Path:
        return self.discovery / "selection_report.html"

    @property
    def selected_analyses(self) -> Path:
        """JSON list of analyses (population label -> study) chosen for the run."""
        return self.discovery / "selected_analyses.json"

    # -- GWAS --------------------------------------------------------------
    def gwas_dir(self, analysis: str) -> Path:
        return self.root / "gwas" / analysis

    def sumstats(self, analysis: str) -> Path:
        return self.gwas_dir(analysis) / "sumstats.parquet"

    def significant(self, analysis: str) -> Path:
        return self.gwas_dir(analysis) / "significant.parquet"

    def gwas_qc(self, analysis: str) -> Path:
        return self.gwas_dir(analysis) / "qc.json"

    @property
    def meta_sumstats(self) -> Path:
        return self.gwas_dir(COMBINED_META) / "meta_sumstats.parquet"

    # -- loci ---------------------------------------------------------------
    def loci_dir(self, analysis: str) -> Path:
        return self.root / "loci" / analysis

    def loci(self, analysis: str) -> Path:
        return self.loci_dir(analysis) / "loci.parquet"

    @property
    def combined_loci(self) -> Path:
        return self.root / "loci" / "combined" / "loci.parquet"

    # -- fine-mapping ---------------------------------------------------------
    def finemap_dir(self, analysis: str) -> Path:
        return self.root / "finemapping" / analysis

    def finemap_locus_dir(self, analysis: str, locus_id: str) -> Path:
        return self.finemap_dir(analysis) / "loci" / locus_id

    @property
    def multiancestry_dir(self) -> Path:
        return self.root / "finemapping" / "multiancestry"

    # -- downstream ---------------------------------------------------------
    @property
    def annotation(self) -> Path:
        return self.root / "annotation"

    @property
    def qtl(self) -> Path:
        return self.root / "qtl" / "all_tissues"

    @property
    def splicing(self) -> Path:
        return self.root / "splicing"

    @property
    def coloc(self) -> Path:
        return self.root / "colocalization"

    @property
    def evidence(self) -> Path:
        return self.root / "evidence"

    @property
    def report(self) -> Path:
        return self.root / "report"


def create_run(run_root: Path, phenotype: str, run_id: str | None = None) -> RunLayout:
    base = run_root / slugify(phenotype)
    layout = RunLayout(base / (run_id or new_run_id()))
    layout.root.mkdir(parents=True, exist_ok=True)
    update_latest(base, layout.root)
    return layout


def update_latest(phenotype_dir: Path, run_dir: Path) -> Path:
    """Point ``<phenotype>/latest`` at ``run_dir`` (symlink; text-file fallback)."""
    latest = phenotype_dir / "latest"
    try:
        if latest.is_symlink() or latest.is_file():
            latest.unlink()
        elif latest.is_dir():
            shutil.rmtree(latest)
        os.symlink(run_dir.name, latest, target_is_directory=True)
    except OSError:
        # Windows without developer mode: record the run id instead.
        latest.write_text(run_dir.name + "\n", encoding="utf-8")
    return latest


def resolve_run_dir(path: Path) -> Path:
    """Accept a run directory or a ``latest`` pointer (symlink or text file)."""
    if path.is_file() and path.name == "latest":
        return (path.parent / path.read_text(encoding="utf-8").strip()).resolve()
    return path.resolve()
