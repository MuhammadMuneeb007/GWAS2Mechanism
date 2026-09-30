"""Optional SuSiEx command-line adapter."""

from __future__ import annotations

from pathlib import Path

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


def run_susiex(config_file: Path, output_prefix: Path, tools: ToolResolver, threads: int) -> Path:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    run(
        [*tools.require("susiex"), "--config", str(config_file), "--out", str(output_prefix), "--threads", str(threads)],
        log_file=output_prefix.with_suffix(".log"),
    )
    result = output_prefix.with_suffix(".snp")
    if not result.exists():
        raise RuntimeError(f"SuSiEx did not create {result}")
    return result
