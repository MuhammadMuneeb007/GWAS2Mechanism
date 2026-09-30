"""Download only the selected summary-statistics files (cached, verified)."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import polars as pl

from gwas2mechanism.config import Config
from gwas2mechanism.paths import RunLayout
from gwas2mechanism.utils.download import Downloader, fetch_range, is_local
from gwas2mechanism.utils.tools import ToolResolver

log = logging.getLogger(__name__)


def _md5_lookup(files: pl.DataFrame, accession: str, filename: str) -> str | None:
    md5_rows = files.filter((pl.col("study_accession") == accession) & (pl.col("kind") == "md5"))
    for url in md5_rows["url"].to_list() if md5_rows.height else []:
        try:
            text = fetch_range(url, 65536).decode("utf-8", "replace")
        except OSError:
            continue
        for line in text.splitlines():
            parts = line.split()
            if (
                len(parts) >= 2
                and parts[-1].lstrip("*./").endswith(filename)
                and re.fullmatch(r"[0-9a-fA-F]{32}", parts[0])
            ):
                return parts[0].lower()
    return None


def _gwaslab_download(accession: str, out_dir: Path, harmonised: bool) -> Path | None:
    try:
        import gwaslab as gl  # pandas-boundary library
    except ImportError:
        return None
    try:
        path = gl.download_sumstats(accession, output_dir=str(out_dir), harmonised=harmonised, overwrite=False, verbose=False)
    except Exception as exc:  # GWASLab API changes between versions
        log.warning("GWASLab download of %s failed (%s); using built-in downloader", accession, exc)
        return None
    return Path(path) if path and Path(path).exists() else None


def download_analysis(analysis: dict[str, Any], layout: RunLayout, config: Config, tools: ToolResolver) -> dict[str, Any]:
    """Fetch the selected file for one analysis; returns provenance."""
    accession = analysis["study_accession"]
    url = analysis.get("selected_file")
    if not url:
        raise FileNotFoundError(f"{accession}: no summary-statistics file was identified during discovery")
    cache = config.resources_path / "gwas" / accession
    cache.mkdir(parents=True, exist_ok=True)
    filename = url.replace("\\", "/").rstrip("/").split("/")[-1]
    files = pl.read_parquet(layout.candidate_files) if layout.candidate_files.exists() else pl.DataFrame()

    method = "builtin"
    path: Path | None = None
    if not is_local(url) and config.harmonization.engine in ("auto", "gwaslab"):
        path = _gwaslab_download(accession, cache, harmonised=bool(analysis.get("harmonised_available")))
        method = "gwaslab" if path else method
    md5 = None
    if path is None:
        md5 = _md5_lookup(files, accession, filename) if not files.is_empty() else None
        aria2 = tools.command("aria2c")
        result = Downloader(use_aria2=aria2).fetch(url, cache / filename, md5=md5)
        path = result.path
    # GWAS-SSF metadata sidecar (tiny) for build / file-type provenance.
    meta_path = None
    if not files.is_empty():
        meta = files.filter((pl.col("study_accession") == accession) & (pl.col("kind") == "meta_yaml"))
        if meta.height:
            meta_url = meta["url"][0]
            meta_path = cache / meta_url.replace("\\", "/").split("/")[-1]
            try:
                Downloader().fetch(meta_url, meta_path)
            except OSError as exc:
                log.warning("Could not fetch metadata sidecar %s: %s", meta_url, exc)
                meta_path = None
    return {
        "analysis_id": analysis["analysis_id"],
        "study_accession": accession,
        "source_url": url,
        "local_path": str(path),
        "meta_yaml": str(meta_path) if meta_path else None,
        "md5": md5,
        "download_method": method,
    }
