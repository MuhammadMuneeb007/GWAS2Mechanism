"""Typer command-line interface for ``gwas2m``."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from gwas2mechanism import __version__
from gwas2mechanism.config import bundled_path, load_config, parse_set_overrides, save_config
from gwas2mechanism.discovery.service import discover as discover_studies
from gwas2mechanism.paths import create_run, resolve_run_dir
from gwas2mechanism.pipeline import run_synthetic_smoke
from gwas2mechanism.provenance import ResourceManifest, RunManifest
from gwas2mechanism.reporting.report import generate_report
from gwas2mechanism.utils.runtime import resolve_threads
from gwas2mechanism.utils.tools import ToolResolver

app = typer.Typer(help="Phenotype-agnostic, population-aware GWAS-to-mechanism framework.", no_args_is_help=True)
console = Console()


def _config(paths: list[Path] | None, sets: list[str] | None, phenotype: str | None = None):
    overrides = parse_set_overrides(sets)
    if phenotype:
        overrides["phenotype"] = phenotype
    return load_config(paths, overrides)


@app.callback()
def main(version: Annotated[bool, typer.Option("--version", help="Show version and exit.")] = False) -> None:
    if version:
        console.print(__version__)
        raise typer.Exit()


@app.command()
def doctor(config: Annotated[list[Path] | None, typer.Option("--config", exists=True)] = None) -> None:
    """Check Python packages, external tools, and configuration."""
    cfg = _config(config, None)
    tools = ToolResolver(cfg.tools)
    table = Table("Component", "Status", "Version/path")
    table.add_row("configuration", "OK", str(bundled_path("config", "default.yaml")))
    for name in cfg.tools:
        command = tools.command(name)
        table.add_row(name, "OK" if command else "optional/missing", tools.version(name) or (" ".join(command) if command else "-"))
    console.print(table)
    console.print("[green]Core configuration is valid.[/green] Missing external tools are required only for their stages.")


@app.command()
def discover(
    phenotype: Annotated[str, typer.Option("--phenotype", help="Phenotype name or ontology label.")],
    populations: Annotated[list[str] | None, typer.Option("--populations")] = None,
    all_studies: Annotated[bool, typer.Option("--all-studies")] = False,
    config: Annotated[list[Path] | None, typer.Option("--config", exists=True)] = None,
    set_value: Annotated[list[str] | None, typer.Option("--set")] = None,
    run_root: Annotated[Path | None, typer.Option("--run-root")] = None,
    offline: Annotated[bool, typer.Option("--offline", help="Disable ontology network lookup; Catalog discovery may still need a local bulk source.")] = False,
) -> None:
    """Resolve a phenotype, enumerate candidates, rank, and select studies."""
    overrides = parse_set_overrides(set_value)
    overrides["phenotype"] = phenotype
    overrides.setdefault("discovery", {})["all_studies"] = all_studies
    if run_root:
        overrides["run_root"] = str(run_root)
    cfg = load_config(config, overrides)
    layout = create_run(Path(cfg.run_root), phenotype)
    save_config(cfg, layout.config)
    manifest = RunManifest(layout.manifest)
    manifest.initialise(phenotype=phenotype, run_id=layout.root.name, command_line=sys.argv, parameters=cfg.dump())
    result = discover_studies(phenotype, cfg, layout, populations, offline)
    manifest.update("selected_studies", result.filter(result["selected"]).to_dicts() if "selected" in result.columns else [])
    console.print(f"Discovery written to [bold]{layout.discovery}[/bold]")


def run_pipeline(
    phenotype: Annotated[str, typer.Option("--phenotype")],
    populations: Annotated[list[str] | None, typer.Option("--populations")] = None,
    tissues: Annotated[str, typer.Option("--tissues")] = "all",
    mode: Annotated[str, typer.Option("--mode")] = "fast",
    threads: Annotated[int | None, typer.Option("--threads", min=1)] = None,
    combined: Annotated[bool, typer.Option("--combined/--no-combined")] = True,
    config: Annotated[list[Path] | None, typer.Option("--config", exists=True)] = None,
    set_value: Annotated[list[str] | None, typer.Option("--set")] = None,
    synthetic: Annotated[bool, typer.Option("--synthetic", help="Run the tiny download-free end-to-end fixture.")] = False,
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Build and print the Snakemake command only.")] = False,
) -> None:
    """Run or resume the complete restartable workflow."""
    overrides = parse_set_overrides(set_value)
    overrides.update({"phenotype": phenotype, "mode": mode})
    overrides.setdefault("qtl", {})["tissues"] = tissues
    overrides.setdefault("combined", {})["fixed_effect_meta"] = combined
    # ``auto`` is a CLI sentinel, not a biological population label.
    if populations and [value.upper() for value in populations] != ["AUTO"]:
        overrides.setdefault("populations", {})["mode"] = "explicit"
        overrides["populations"]["allowed"] = [value.upper() for value in populations]
    if threads:
        overrides.setdefault("performance", {})["threads"] = threads
    cfg = load_config(config, overrides)
    if synthetic:
        layout = run_synthetic_smoke(cfg, phenotype)
        console.print(f"Synthetic end-to-end run completed: [bold]{layout.root}[/bold]")
        return
    layout = create_run(Path(cfg.run_root), phenotype)
    save_config(cfg, layout.config)
    command = [*ToolResolver(cfg.tools).require("snakemake"), "--snakefile", str(bundled_path("Snakefile")), "--directory", str(Path.cwd()), "--configfile", str(layout.config), "--config", f"run_dir={layout.root}", f"phenotype={phenotype}", "--cores", str(resolve_threads(cfg.performance.threads)), "--rerun-incomplete", "--printshellcmds"]
    if dry_run:
        command.append("--dry-run")
    console.print(" ".join(map(str, command)))
    completed = subprocess.run(command, check=False)
    raise typer.Exit(completed.returncode)


app.command("run")(run_pipeline)


@app.command()
def resume(run: Annotated[Path, typer.Option("--run", exists=True)], threads: Annotated[int | None, typer.Option("--threads", min=1)] = None) -> None:
    """Resume a run from verified Snakemake outputs."""
    run_dir = resolve_run_dir(run)
    cfg = load_config([run_dir / "config.yaml"])
    command = [*ToolResolver(cfg.tools).require("snakemake"), "--snakefile", str(bundled_path("Snakefile")), "--configfile", str(run_dir / "config.yaml"), "--config", f"run_dir={run_dir}", f"phenotype={cfg.phenotype}", "--cores", str(threads or resolve_threads(cfg.performance.threads)), "--rerun-incomplete"]
    raise typer.Exit(subprocess.run(command, check=False).returncode)


@app.command("report")
def report_command(run: Annotated[Path, typer.Option("--run", exists=True)]) -> None:
    html_path, markdown_path, summary_path = generate_report(resolve_run_dir(run))
    console.print(f"Report: {html_path}\nMarkdown: {markdown_path}\nSummary: {summary_path}")


@app.command()
def resources(config: Annotated[list[Path] | None, typer.Option("--config", exists=True)] = None) -> None:
    cfg = _config(config, None)
    manifest = ResourceManifest(cfg.resources_path / "manifest.tsv").read()
    console.print(manifest if not manifest.is_empty() else f"No cached resources recorded under {cfg.resources_path}")


@app.command()
def setup(
    reference: Annotated[bool, typer.Option("--reference")] = False,
    gtex: Annotated[bool, typer.Option("--gtex")] = False,
    vep: Annotated[bool, typer.Option("--vep")] = False,
    splice: Annotated[bool, typer.Option("--splice")] = False,
    all_resources: Annotated[bool, typer.Option("--all")] = False,
    full: Annotated[bool, typer.Option("--full")] = False,
    populations: Annotated[list[str] | None, typer.Option("--populations")] = None,
    config: Annotated[list[Path] | None, typer.Option("--config", exists=True)] = None,
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Print the resource plan without downloading.")] = False,
) -> None:
    """Download, validate and cache versioned scientific resources."""
    from gwas2mechanism.setup_resources import setup_resources

    cfg = _config(config, None)
    chosen = all_resources or full
    plan = setup_resources(cfg, reference=reference or chosen, gtex=gtex or chosen, vep=vep or chosen, splice=splice or chosen, populations=populations or cfg.populations.allowed, full=full, dry_run=dry_run)
    console.print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    app()
