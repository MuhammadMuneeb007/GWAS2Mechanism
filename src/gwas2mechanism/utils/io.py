"""Columnar I/O helpers: Parquet+ZSTD intermediates, atomic writes, TSV exports."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import polars as pl

PARQUET_COMPRESSION = "zstd"
PARQUET_LEVEL = 3


def configure(compression: str = "zstd", level: int = 3) -> None:
    global PARQUET_COMPRESSION, PARQUET_LEVEL
    PARQUET_COMPRESSION = compression
    PARQUET_LEVEL = level


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _tmp_path(path: Path) -> Path:
    ensure_dir(path.parent)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".part", dir=path.parent)
    os.close(fd)
    return Path(name)


def write_parquet(frame: pl.DataFrame | pl.LazyFrame, path: Path) -> Path:
    """Write Parquet atomically (tmp file + rename) with ZSTD compression."""
    tmp = _tmp_path(path)
    try:
        if isinstance(frame, pl.LazyFrame):
            frame.sink_parquet(
                tmp, compression=PARQUET_COMPRESSION, compression_level=PARQUET_LEVEL
            )
        else:
            frame.write_parquet(
                tmp, compression=PARQUET_COMPRESSION, compression_level=PARQUET_LEVEL
            )
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return path


def write_tsv(frame: pl.DataFrame, path: Path) -> Path:
    """TSV is an *export* format only; nested columns are JSON-encoded."""
    out = frame
    nested = [c for c, dt in frame.schema.items() if isinstance(dt, (pl.List, pl.Struct))]
    if nested:
        out = frame.with_columns(
            [
                pl.col(c).map_elements(
                    lambda v: json.dumps(v.to_list() if hasattr(v, "to_list") else v),
                    return_dtype=pl.String,
                )
                for c in nested
            ]
        )
    tmp = _tmp_path(path)
    try:
        out.write_csv(tmp, separator="\t", null_value="NA")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return path


def atomic_write_text(path: Path, text: str) -> Path:
    tmp = _tmp_path(path)
    try:
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return path


def atomic_write_json(path: Path, payload: Any) -> Path:
    return atomic_write_text(path, json.dumps(payload, indent=2, default=str))


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def read_parquet_or_empty(path: Path, schema: dict[str, Any] | None = None) -> pl.DataFrame:
    if path.exists():
        return pl.read_parquet(path)
    return pl.DataFrame(schema=schema or {})


def concat_parquets(paths: Iterable[Path]) -> pl.DataFrame:
    existing = [p for p in paths if p.exists()]
    if not existing:
        return pl.DataFrame()
    return pl.concat([pl.read_parquet(p) for p in existing], how="diagonal_relaxed")


def empty_frame(columns: dict[str, pl.DataType]) -> pl.DataFrame:
    return pl.DataFrame(schema=columns)
