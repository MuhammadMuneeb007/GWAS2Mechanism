"""Run manifest (manifest.yaml) and resource manifest (manifest.tsv)."""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as md
import platform
import sys
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from gwas2mechanism import __version__
from gwas2mechanism.utils.io import atomic_write_text
from gwas2mechanism.utils.locking import resource_lock

PYTHON_PACKAGES = (
    "polars",
    "pyarrow",
    "duckdb",
    "numpy",
    "scipy",
    "pandas",
    "pydantic",
    "gwaslab",
    "MultiSuSiE",
    "gwaspoker",
    "snakemake",
    "plotly",
)

RESOURCE_COLUMNS = {
    "resource": pl.String,
    "version": pl.String,
    "source": pl.String,
    "local_path": pl.String,
    "size": pl.Int64,
    "checksum": pl.String,
    "download_date": pl.String,
    "status": pl.String,
    "expected_size": pl.Int64,
    "actual_size": pl.Int64,
    "verification_timestamp": pl.String,
    "chromosome": pl.String,
    "ancestry": pl.String,
    "started_at": pl.String,
    "completed_at": pl.String,
    "software_version": pl.String,
    "preparation_command": pl.String,
}


def python_versions() -> dict[str, str | None]:
    out: dict[str, str | None] = {"python": sys.version.split()[0], "gwas2mechanism": __version__}
    for name in PYTHON_PACKAGES:
        try:
            out[name] = md.version(name)
        except md.PackageNotFoundError:
            out[name] = None
    return out


class RunManifest:
    """A YAML document updated incrementally by each stage."""

    def __init__(self, path: Path):
        self.path = path
        self.data: dict[str, Any] = {}
        if path.exists():
            self.data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    def initialise(self, *, phenotype: str, run_id: str, command_line: list[str], parameters: dict[str, Any]) -> None:
        self.data.setdefault("phenotype", phenotype)
        self.data.setdefault("run_id", run_id)
        self.data.setdefault("created", _dt.datetime.now().isoformat(timespec="seconds"))
        self.data["command_line"] = command_line
        self.data["platform"] = platform.platform()
        self.data["software"] = {"python_packages": python_versions()}
        self.data["parameters"] = parameters
        self.data.setdefault("genome_build", "GRCh38")
        self.save()

    def update(self, key: str, value: Any) -> None:
        self.data[key] = value
        self.save()

    def merge(self, key: str, value: dict[str, Any]) -> None:
        current = self.data.get(key) or {}
        current.update(value)
        self.data[key] = current
        self.save()

    def save(self) -> None:
        atomic_write_text(self.path, yaml.safe_dump(self.data, sort_keys=False, default_flow_style=False))


class ResourceManifest:
    """``<resources>/manifest.tsv`` - one row per cached resource (latest wins)."""

    def __init__(self, path: Path):
        self.path = path

    def read(self) -> pl.DataFrame:
        if not self.path.exists():
            return pl.DataFrame(schema=RESOURCE_COLUMNS)
        frame = pl.read_csv(self.path, separator="\t", infer_schema=False, null_values=["NA", ""])
        for name, dtype in RESOURCE_COLUMNS.items():
            if name not in frame.columns:
                frame = frame.with_columns(pl.lit(None, dtype=dtype).alias(name))
            else:
                frame = frame.with_columns(pl.col(name).cast(dtype, strict=False))
        return frame.select(list(RESOURCE_COLUMNS))

    def record(
        self,
        resource: str,
        *,
        version: str | None,
        source: str | None,
        local_path: Path | str,
        size: int | None = None,
        checksum: str | None = None,
        status: str = "OK",
        expected_size: int | None = None,
        actual_size: int | None = None,
        chromosome: str | None = None,
        ancestry: str | None = None,
        started_at: str | None = None,
        completed_at: str | None = None,
        software_version: str | None = None,
        preparation_command: str | None = None,
    ) -> None:
        now = _dt.datetime.now().isoformat(timespec="seconds")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with resource_lock(self.path):
            row = pl.DataFrame(
                [{"resource": resource, "version": version, "source": source, "local_path": str(local_path), "size": size, "checksum": checksum, "download_date": _dt.date.today().isoformat(), "status": status, "expected_size": expected_size, "actual_size": actual_size if actual_size is not None else size, "verification_timestamp": now, "chromosome": chromosome, "ancestry": ancestry, "started_at": started_at, "completed_at": completed_at or now, "software_version": software_version, "preparation_command": preparation_command}],
                schema=RESOURCE_COLUMNS,
            )
            current = self.read().filter(pl.col("resource") != resource)
            out = pl.concat([current, row], how="vertical_relaxed").sort("resource")
            tmp = self.path.with_suffix(".tsv.part")
            out.write_csv(tmp, separator="\t", null_value="NA")
            tmp.replace(self.path)

    def lookup(self, resource: str) -> dict[str, Any] | None:
        df = self.read().filter(pl.col("resource") == resource)
        return df.row(0, named=True) if df.height else None
