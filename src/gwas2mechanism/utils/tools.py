"""Resolve external tool commands.

Each tool entry in the configuration may be:

* ``auto``          -> the executable on ``PATH``; otherwise the dedicated conda
                       environment created by ``gwas2m setup`` via ``mamba run``;
* a string          -> an executable name or path;
* a list of strings -> a full command prefix (e.g. ``["mamba", "run", "-n", "x", "vep"]``
                       or ``["python", "mock.py"]`` in tests).
"""

from __future__ import annotations

import functools
import json
import logging
import re
import shutil
import subprocess

log = logging.getLogger(__name__)

#: tool key -> (executable name, isolated env name or None)
TOOL_SPECS: dict[str, tuple[str, str | None]] = {
    "plink2": ("plink2", None),
    "bcftools": ("bcftools", None),
    "samtools": ("samtools", None),
    "tabix": ("tabix", None),
    "bgzip": ("bgzip", None),
    "aria2c": ("aria2c", None),
    "rscript": ("Rscript", None),
    "vep": ("vep", "gwas2mechanism-vep"),
    "vep_install": ("vep_install", "gwas2mechanism-vep"),
    "spliceai": ("spliceai", "gwas2mechanism-spliceai"),
    "pangolin": ("pangolin", "gwas2mechanism-pangolin"),
    "susiex": ("SuSiEx", None),
    "metal": ("metal", None),
    "mrmega": ("MR-MEGA", None),
    "snakemake": ("snakemake", None),
}

VERSION_ARGS: dict[str, list[str]] = {
    "plink2": ["--version"],
    "bcftools": ["--version"],
    "samtools": ["--version"],
    "rscript": ["--version"],
    "vep": ["--help"],
    "snakemake": ["--version"],
    "aria2c": ["--version"],
}


@functools.lru_cache(maxsize=1)
def _conda_frontend() -> str | None:
    for exe in ("mamba", "micromamba", "conda"):
        if shutil.which(exe):
            return exe
    return None


@functools.lru_cache(maxsize=1)
def _conda_env_names() -> frozenset[str]:
    frontend = _conda_frontend()
    if frontend is None:
        return frozenset()
    try:
        out = subprocess.run(
            [frontend, "env", "list", "--json"], capture_output=True, text=True, timeout=60
        )
        envs = json.loads(out.stdout or "{}").get("envs", [])
    except (OSError, ValueError, subprocess.SubprocessError):
        return frozenset()
    return frozenset(e.replace("\\", "/").rstrip("/").split("/")[-1] for e in envs)


class ToolResolver:
    def __init__(self, tools_config: dict[str, str | list[str]]):
        self._config = dict(tools_config)

    def command(self, name: str) -> list[str] | None:
        """Command prefix for ``name`` or ``None`` when unavailable."""
        spec = self._config.get(name, "auto")
        exe, env_name = TOOL_SPECS.get(name, (name, None))
        if isinstance(spec, list):
            return list(spec)
        if spec != "auto":
            return [spec] if shutil.which(spec) or "/" in spec or "\\" in spec else None
        found = shutil.which(exe)
        if found:
            return [found]
        if env_name and env_name in _conda_env_names():
            frontend = _conda_frontend()
            if frontend:
                return [frontend, "run", "-n", env_name, exe]
        return None

    def available(self, name: str) -> bool:
        return self.command(name) is not None

    def require(self, name: str) -> list[str]:
        cmd = self.command(name)
        if cmd is None:
            exe = TOOL_SPECS.get(name, (name, None))[0]
            raise FileNotFoundError(
                f"Required tool '{exe}' was not found. Run `gwas2m setup` or set tools.{name} "
                f"in your configuration."
            )
        return cmd

    def version(self, name: str) -> str | None:
        cmd = self.command(name)
        if cmd is None:
            return None
        args = VERSION_ARGS.get(name, ["--version"])
        try:
            out = subprocess.run(cmd + args, capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError):
            return None
        text = (out.stdout or "") + (out.stderr or "")
        if name == "vep":
            match = re.search(r"ensembl-vep\s*:\s*([\d.]+)", text)
            return match.group(1) if match else None
        first = next((ln.strip() for ln in text.splitlines() if ln.strip()), "")
        return first[:120] or None
