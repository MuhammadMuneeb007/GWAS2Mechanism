"""Generate and optionally submit an idempotent SLURM resource-setup DAG."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any

from gwas2mechanism.config import Config
from gwas2mechanism.constants import SUPERPOPULATIONS


def _quote(value: object) -> str:
    return shlex.quote(str(value))


def _header(name: str, logs: Path, *, cpus: int, memory: str, partition: str | None = None, array: str | None = None) -> list[str]:
    lines = ["#!/usr/bin/env bash", f"#SBATCH --job-name=g2m-{name}", f"#SBATCH --cpus-per-task={cpus}", f"#SBATCH --mem={memory}", f"#SBATCH --output={logs}/{name}.%A_%a.out" if array else f"#SBATCH --output={logs}/{name}.%j.out", f"#SBATCH --error={logs}/{name}.%A_%a.err" if array else f"#SBATCH --error={logs}/{name}.%j.err"]
    if array:
        lines.append(f"#SBATCH --array={array}")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    lines.extend(["set -Eeuo pipefail", ""])
    return lines


def _preamble(root: Path, resource: str, source: str) -> list[str]:
    resource_root = (root / ".gwas2m" / "resources").resolve()
    return [
        f"cd {_quote(root)}",
        f"echo {_quote(f'Resource: {resource}')}",
        "echo \"Task ID: ${SLURM_ARRAY_TASK_ID:-none}\"",
        "echo \"Hostname: $(hostname)\"",
        "echo \"CPU allocation: ${SLURM_CPUS_PER_TASK:-1}\"",
        f"echo {_quote(f'Resource root: {resource_root}')}",
        f"echo {_quote(f'Source: {source}')}",
    ]


def generate_slurm_jobs(config: Config, *, reference: bool, gtex: bool, vep: bool, splice: bool, populations: list[str], partition: str | None = None) -> dict[str, Any]:
    invalid = set(populations) - set(SUPERPOPULATIONS)
    if invalid:
        raise ValueError(f"Unsupported populations: {sorted(invalid)}")
    if partition and not re.fullmatch(r"[A-Za-z0-9_.-]+", partition):
        raise ValueError("SLURM partition contains unsafe characters")
    root = Path(os.environ.get("GWAS2M_PROJECT_ROOT", Path.cwd())).resolve()
    jobs = root / ".gwas2m" / "setup_jobs"
    logs = root / ".gwas2m" / "setup_logs"
    jobs.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    executable = root / ".gwas2m" / "envs" / "gwas2mechanism" / "bin" / "gwas2m"
    command = _quote(executable)
    chromosomes = config.genome.chromosomes
    setup_cfg = config.performance.setup

    scripts: dict[str, Path] = {}

    def write(key: str, filename: str, lines: list[str]) -> None:
        path = jobs / filename
        path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        path.chmod(0o750)
        scripts[key] = path

    lines = _header("environment", logs, cpus=1, memory="2G", partition=partition)
    lines += _preamble(root, "environment check", "project-local conda environments")
    lines += [f"{command} setup-task --task environment-check", "echo 'Environment check complete.'"]
    write("environment", "00_environment_check.sh", lines)

    if reference or splice:
        lines = _header("genome", logs, cpus=2, memory="8G", partition=partition)
        lines += _preamble(root, "GENCODE GRCh38", "GENCODE release 50")
        lines += [f"{command} setup-task --task genome", "echo 'GENCODE genome setup complete.'"]
        write("genome", "01_genome_gencode.sh", lines)

    if reference:
        chrom_values = " ".join(_quote(value) for value in chromosomes)
        lines = _header("1000g_download", logs, cpus=2, memory="8G", partition=partition, array=f"1-{len(chromosomes)}%{setup_cfg.slurm_download_array_limit}")
        lines += _preamble(root, "1000 Genomes raw chromosome", config.reference.base_url)
        lines += [f"CHROMS=({chrom_values})", "CHR=${CHROMS[$((SLURM_ARRAY_TASK_ID-1))]}", "echo \"Chromosome: $CHR\"", f"{command} setup-task --task reference-download --chromosome \"$CHR\"", "echo \"1000G raw chromosome $CHR complete.\""]
        write("reference_download", "02_1000g_download.sh", lines)

        lines = _header("1000g_convert", logs, cpus=4, memory="32G", partition=partition, array=f"1-{len(chromosomes)}%{setup_cfg.slurm_convert_array_limit}")
        lines += _preamble(root, "1000 Genomes ALL PGEN intermediate", "verified raw chromosome VCF")
        lines += [f"CHROMS=({chrom_values})", "CHR=${CHROMS[$((SLURM_ARRAY_TASK_ID-1))]}", "echo \"Chromosome: $CHR\"", f"{command} setup-task --task reference-convert --chromosome \"$CHR\"", "echo \"ALL PGEN chromosome $CHR complete.\""]
        write("reference_convert", "03_1000g_convert_all.sh", lines)

        tasks = len(chromosomes) * len(populations)
        pop_values = " ".join(_quote(value) for value in populations)
        lines = _header("1000g_population", logs, cpus=2, memory="16G", partition=partition, array=f"1-{tasks}%{setup_cfg.slurm_prepare_array_limit}")
        lines += _preamble(root, "ancestry-specific 1000 Genomes PGEN", "ALL PGEN intermediate")
        lines += [f"CHROMS=({chrom_values})", f"POPS=({pop_values})", f"N_CHR={len(chromosomes)}", "INDEX=$((SLURM_ARRAY_TASK_ID-1))", "POP=${POPS[$((INDEX/N_CHR))]}", "CHR=${CHROMS[$((INDEX%N_CHR))]}", "echo \"Chromosome: $CHR\"", "echo \"Ancestry: $POP\"", f"{command} setup-task --task reference-population --chromosome \"$CHR\" --population \"$POP\"", "echo \"1000G $POP chromosome $CHR complete.\""]
        write("reference_population", "04_1000g_prepare_populations.sh", lines)

    if gtex:
        lines = _header("gtex", logs, cpus=4, memory="16G", partition=partition)
        lines += _preamble(root, "Adult GTEx compact SuSiE eQTL/sQTL", config.qtl.bucket_download_url)
        lines += [f"{command} setup-task --task gtex", "echo 'GTEx setup complete.'"]
        write("gtex", "05_gtex.sh", lines)

    if vep:
        lines = _header("vep", logs, cpus=4, memory="24G", partition=partition)
        lines += _preamble(root, "Ensembl VEP cache", "Ensembl VEP installer")
        lines += [f"{command} setup-task --task vep", "echo 'VEP cache setup complete.'"]
        write("vep", "06_vep.sh", lines)

    if splice:
        lines = _header("splicing", logs, cpus=4, memory="24G", partition=partition)
        lines += _preamble(root, "SpliceAI and Pangolin", "GENCODE release 50")
        lines += [f"{command} setup-task --task splice", "echo 'Splicing resource setup complete.'"]
        write("splice", "07_splicing.sh", lines)

    pop_args = " ".join(f"--populations {_quote(pop)}" for pop in populations)
    component_args = " ".join(f"--{name}" for name, enabled in (("reference", reference), ("gtex", gtex), ("vep", vep), ("splice", splice)) if enabled)
    lines = _header("verify", logs, cpus=2, memory="8G", partition=partition)
    lines += _preamble(root, "resource validation", "local resource cache")
    lines += [f"{command} setup-task --task verify {component_args} {pop_args}", "echo 'Resource validation complete.'"]
    write("verify", "08_verify.sh", lines)

    submit = ["#!/usr/bin/env bash", "set -Eeuo pipefail", f"cd {_quote(root)}", "submit() { sbatch --parsable \"$@\" | cut -d';' -f1; }", f"ENV_JOB=$(submit {_quote(scripts['environment'])})", "echo \"environment: $ENV_JOB\""]
    terminal: list[str] = []
    if "genome" in scripts:
        submit += [f"GENOME_JOB=$(submit --dependency=afterok:$ENV_JOB {_quote(scripts['genome'])})", "echo \"genome: $GENOME_JOB\""]
        terminal.append("$GENOME_JOB")
    if reference:
        submit += [f"RAW_JOB=$(submit --dependency=afterok:$ENV_JOB {_quote(scripts['reference_download'])})", f"ALL_JOB=$(submit --dependency=afterok:$RAW_JOB {_quote(scripts['reference_convert'])})", f"POP_JOB=$(submit --dependency=afterok:$ALL_JOB {_quote(scripts['reference_population'])})", "echo \"1000G raw/ALL/population: $RAW_JOB $ALL_JOB $POP_JOB\""]
        terminal.append("$POP_JOB")
    if gtex:
        submit += [f"GTEX_JOB=$(submit --dependency=afterok:$ENV_JOB {_quote(scripts['gtex'])})", "echo \"GTEx: $GTEX_JOB\""]
        terminal.append("$GTEX_JOB")
    if vep:
        submit += [f"VEP_JOB=$(submit --dependency=afterok:$ENV_JOB {_quote(scripts['vep'])})", "echo \"VEP: $VEP_JOB\""]
        terminal.append("$VEP_JOB")
    if splice:
        dependency = "$ENV_JOB:$GENOME_JOB" if "genome" in scripts else "$ENV_JOB"
        submit += [f"SPLICE_JOB=$(submit --dependency=afterok:{dependency} {_quote(scripts['splice'])})", "echo \"splicing: $SPLICE_JOB\""]
        terminal.append("$SPLICE_JOB")
    dependency = ":".join(terminal) or "$ENV_JOB"
    submit += [f"VERIFY_JOB=$(submit --dependency=afterok:{dependency} {_quote(scripts['verify'])})", "echo \"validation: $VERIFY_JOB\"", "echo 'SLURM setup DAG submitted.'"]
    submit_path = jobs / "submit_all.sh"
    submit_path.write_text("\n".join(submit) + "\n", encoding="utf-8", newline="\n")
    submit_path.chmod(0o750)
    return {"jobs_dir": str(jobs), "logs_dir": str(logs), "scripts": {key: str(path) for key, path in scripts.items()}, "submit_script": str(submit_path), "populations": populations, "partition": partition}


def submit_slurm_jobs(plan: dict[str, Any]) -> str:
    if shutil.which("sbatch") is None:
        raise FileNotFoundError("sbatch was not found; jobs were generated but not submitted")
    completed = subprocess.run(["bash", str(plan["submit_script"])], check=True, text=True, capture_output=True)
    return completed.stdout
