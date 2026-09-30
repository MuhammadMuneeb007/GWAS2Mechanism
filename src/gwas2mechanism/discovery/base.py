"""DiscoveryBackend interface plus shared file listing and header probing."""

from __future__ import annotations

import gzip
import logging
import os
import re
import zlib
from abc import ABC, abstractmethod
from pathlib import Path
from urllib.parse import urljoin

import httpx
import yaml

from gwas2mechanism.discovery.models import (
    CandidateFile,
    CandidateStudy,
    OntologyTerm,
    PhenotypeResolution,
)
from gwas2mechanism.sumstats.schema import capabilities, map_columns
from gwas2mechanism.utils.download import USER_AGENT, fetch_range, is_local, local_path

log = logging.getLogger(__name__)

DATA_SUFFIXES = (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".txt", ".csv", ".gz", ".zip", ".bgz", ".vcf.gz")


class DiscoveryBackend(ABC):
    """Enumerates candidate studies and their files for a resolved phenotype."""

    name: str = "base"

    @abstractmethod
    def ontology_terms(self, phenotype: str) -> list[OntologyTerm]:
        """Terms the backend itself knows for ``phenotype`` (exact label matches are used)."""

    @abstractmethod
    def search_studies(self, resolution: PhenotypeResolution) -> list[CandidateStudy]:
        """Candidate studies (metadata only)."""

    def list_files(self, study: CandidateStudy) -> list[CandidateFile]:
        if not study.summary_stats_url:
            return []
        return list_study_directory(study.study_accession, study.summary_stats_url)

    def probe(self, file: CandidateFile, n_bytes: int) -> CandidateFile:
        return probe_file(file, n_bytes)

    def close(self) -> None:  # pragma: no cover - optional
        """Release backend resources when a backend owns any."""
        return None


# ----------------------------------------------------------------------------
# Directory listing (HTTP index pages or local mirror)
# ----------------------------------------------------------------------------

_HREF = re.compile(r'href="([^"?#]+)"', re.IGNORECASE)


def _classify(filename: str, in_harmonised: bool) -> str:
    lower = filename.lower()
    if lower.endswith("-meta.yaml") or lower.endswith(".yaml"):
        return "meta_yaml"
    if "md5" in lower:
        return "md5"
    if lower.endswith(DATA_SUFFIXES):
        return "harmonised" if (in_harmonised or ".h.tsv" in lower) else "raw"
    return "other"


def _build_from_name(filename: str) -> str | None:
    match = re.search(r"build(GRCh3[78]|hg1[89]|hg38|b3[78])", filename, re.IGNORECASE)
    if not match:
        return None
    token = match.group(1).lower()
    return "GRCh38" if token in ("grch38", "hg38", "b38") else "GRCh37"


def list_study_directory(accession: str, url: str, timeout: float = 60.0) -> list[CandidateFile]:
    url = url.replace("http://", "https://", 1)
    if not url.endswith("/"):
        url += "/"
    files: list[CandidateFile] = []
    for name, full, harmonised in _walk(url, timeout):
        kind = _classify(name, harmonised)
        if kind == "other":
            continue
        is_h = kind == "harmonised"
        files.append(
            CandidateFile(
                study_accession=accession,
                url=full,
                filename=name,
                kind=kind,
                is_harmonised=is_h,
                build="GRCh38" if is_h else _build_from_name(name),
            )
        )
    return files


def _walk(url: str, timeout: float) -> list[tuple[str, str, bool]]:
    out: list[tuple[str, str, bool]] = []
    if is_local(url):
        root = local_path(url.rstrip("/"))
        if not root.exists():
            return out
        for path in sorted(root.rglob("*")):
            if path.is_file():
                out.append((path.name, str(path), "harmonised" in path.parts))
        return out
    try:
        with httpx.Client(timeout=timeout, headers={"User-Agent": USER_AGENT}, follow_redirects=True) as client:
            listing = client.get(url)
            listing.raise_for_status()
            for href in _HREF.findall(listing.text):
                if href.startswith(("/", "..", "http")) and not href.startswith(url):
                    continue
                full = urljoin(url, href)
                name = href.rstrip("/").split("/")[-1]
                if href.endswith("/"):
                    if name.lower() == "harmonised":
                        sub = client.get(full)
                        if sub.status_code == 200:
                            for inner in _HREF.findall(sub.text):
                                if inner.startswith(("/", "..", "http")) or inner.endswith("/"):
                                    continue
                                out.append((inner, urljoin(full, inner), True))
                    continue
                out.append((name, full, False))
    except httpx.HTTPError as exc:
        log.warning("Could not list %s: %s", url, exc)
    return out


# ----------------------------------------------------------------------------
# Bounded header probe
# ----------------------------------------------------------------------------


def decompress_prefix(data: bytes) -> bytes:
    if data[:2] == b"\x1f\x8b":
        decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            return decompressor.decompress(data)
        except zlib.error:
            try:
                return gzip.decompress(data)
            except (OSError, EOFError):
                return b""
    return data


def parse_header(text: str) -> list[str]:
    for line in text.splitlines():
        if not line.strip() or line.startswith("##"):
            continue
        for sep in ("\t", ",", None):
            parts = line.split(sep) if sep else line.split()
            if len(parts) >= 4:
                return [p.strip().strip('"') for p in parts]
        return []
    return []


def probe_file(file: CandidateFile, n_bytes: int) -> CandidateFile:
    """Read at most ``n_bytes`` (HTTP Range) to recover the header and metadata."""
    try:
        if file.kind == "meta_yaml":
            raw = fetch_range(file.url, min(n_bytes, 65536))
            meta = yaml.safe_load(raw.decode("utf-8", "replace")) or {}
            file.build = _normalise_build(meta.get("genome_assembly")) or file.build
            file.file_type = meta.get("file_type")
            file.is_harmonised = bool(meta.get("is_harmonised", file.is_harmonised))
            file.probe_status = "meta_yaml"
            file.probe_bytes = len(raw)
            return file
        raw = fetch_range(file.url, n_bytes)
        text = decompress_prefix(raw).decode("utf-8", "replace")
        header = parse_header(text)
        file.header_columns = header
        file.probe_bytes = len(raw)
        if header:
            caps = capabilities(map_columns(header))
            file.has_beta, file.has_or, file.has_se = caps["has_beta"], caps["has_or"], caps["has_se"]
            file.has_eaf, file.has_n, file.has_alleles = caps["has_eaf"], caps["has_n"], caps["has_alleles"]
            file.probe_status = "header_parsed"
        else:
            file.probe_status = "header_not_found"
    except (httpx.HTTPError, OSError, yaml.YAMLError) as exc:
        file.probe_status = f"probe_failed: {type(exc).__name__}"
    return file


def _normalise_build(value: object) -> str | None:
    if not value:
        return None
    text = str(value).lower()
    if "38" in text:
        return "GRCh38"
    if "37" in text or "19" in text:
        return "GRCh37"
    return str(value)


def gcst_directory(base: str, accession: str) -> str:
    """GWAS Catalog 1000-study bucket directory for an accession (legacy logic)."""
    digits = accession.replace("GCST", "")
    number = int(digits)
    width = max(len(digits), 6)
    start = ((number - 1) // 1000) * 1000 + 1
    block = f"GCST{start:0{width}d}-GCST{start + 999:0{width}d}"
    if is_local(base):
        return str(Path(local_path(base)) / block / accession) + os.sep
    return f"{base.rstrip('/')}/{block}/{accession}/"
