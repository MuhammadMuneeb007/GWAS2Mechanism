"""Versioned, resumable resource setup used by ``gwas2m setup``."""

from __future__ import annotations

import datetime as dt
import sys
import tarfile
from pathlib import Path
from typing import Any

import httpx

from gwas2mechanism.config import Config
from gwas2mechanism.constants import SUPERPOPULATIONS
from gwas2mechanism.provenance import ResourceManifest
from gwas2mechanism.qtl.gtex import discover_release
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.utils.download import Downloader
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.runtime import resolve_threads
from gwas2mechanism.utils.tools import ToolResolver

FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"


def _safe_extract(archive: Path, destination: Path) -> None:
    """Extract a trusted scientific archive while rejecting path traversal."""
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive, "r:*") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if root != target and root not in target.parents:
                raise ValueError(f"Unsafe archive member: {member.name}")
        handle.extractall(destination)


def _gtex_objects(config: Config, release: str) -> list[dict[str, Any]]:
    response = httpx.get(config.qtl.bucket_listing_url, params={"prefix": f"{release}/"}, timeout=60)
    response.raise_for_status()
    objects: list[dict[str, Any]] = []
    token: str | None = None
    while True:
        payload = response.json()
        objects.extend(payload.get("items", []))
        token = payload.get("nextPageToken")
        if not token:
            break
        response = httpx.get(config.qtl.bucket_listing_url, params={"prefix": f"{release}/", "pageToken": token}, timeout=60)
        response.raise_for_status()
    wanted = []
    for item in objects:
        name = str(item.get("name", ""))
        lower = name.lower()
        if "susie" in lower and any(qtl_type in lower for qtl_type in ("eqtl", "sqtl")) and (lower.endswith(".tar") or lower.endswith(".tar.gz") or lower.endswith(".parquet")):
            wanted.append(item)
    return wanted


def setup_resources(
    config: Config,
    *,
    reference: bool,
    gtex: bool,
    vep: bool,
    splice: bool,
    populations: list[str],
    full: bool,
    dry_run: bool,
) -> dict[str, Any]:
    invalid = set(populations) - set(SUPERPOPULATIONS)
    if invalid:
        raise ValueError(f"No 1000 Genomes super-population reference for {sorted(invalid)}")
    root = config.resources_path.resolve()
    plan: dict[str, Any] = {
        "cache": str(root),
        "date": dt.date.today().isoformat(),
        "reference": reference,
        "gtex": gtex,
        "vep": vep,
        "splice": splice,
        "populations": populations,
        "full": full,
        "dry_run": dry_run,
    }
    if dry_run:
        return plan
    root.mkdir(parents=True, exist_ok=True)
    manifest = ResourceManifest(root / "manifest.tsv")
    tools = ToolResolver(config.tools)
    threads = resolve_threads(config.performance.threads)
    # FASTA/GTF are shared by VEP, build validation and both splice predictors.
    if reference or vep or splice:
        genome = root / "genome" / "GRCh38"
        fasta_gz = genome / "GRCh38.primary_assembly.fa.gz"
        gtf = genome / "gencode.v50.annotation.gtf.gz"
        fasta_result = Downloader(use_aria2=tools.command("aria2c")).fetch(FASTA_URL, fasta_gz)
        gtf_result = Downloader(use_aria2=tools.command("aria2c")).fetch(GTF_URL, gtf)
        manifest.record("genome/GRCh38/fasta", version="GENCODE-50", source=FASTA_URL, local_path=fasta_gz, size=fasta_result.size, checksum=fasta_result.md5)
        manifest.record("genome/GRCh38/gtf", version="GENCODE-50", source=GTF_URL, local_path=gtf, size=gtf_result.size, checksum=gtf_result.md5)
        if tools.available("samtools"):
            run([*tools.require("samtools"), "faidx", str(fasta_gz)], log_file=genome / "samtools-faidx.log")
    if reference:
        panel = ReferencePanel(config.reference, root, tools)
        panel.download(config.genome.chromosomes, manifest, workers=min(4, threads))
        for population in populations:
            panel.prepare(population, config.genome.chromosomes, threads, manifest)
    if gtex:
        release = discover_release(config.qtl.bucket_listing_url) if config.qtl.release == "latest" else config.qtl.release.lower()
        objects = _gtex_objects(config, release)
        if not objects:
            raise RuntimeError(f"No compact SuSiE eQTL/sQTL archives found for Adult GTEx {release}")
        destination = root / "gtex" / release
        total = 0
        for item in objects:
            name = str(item["name"])
            path = destination / Path(name).name
            url = f"{config.qtl.bucket_download_url.rstrip('/')}/{name}"
            result = Downloader(use_aria2=tools.command("aria2c")).fetch(url, path, expected_size=int(item["size"]) if item.get("size") else None, md5=item.get("md5Hash"))
            total += result.size
            if path.name.endswith((".tar", ".tar.gz", ".tgz")):
                extracted = destination / "extracted" / path.name.removesuffix(".gz").removesuffix(".tar")
                marker = extracted / ".complete"
                if not marker.exists():
                    _safe_extract(path, extracted)
                    marker.write_text("OK\n", encoding="utf-8")
        manifest.record(f"gtex/{release}/susie", version=release, source=config.qtl.bucket_download_url, local_path=destination, size=total, checksum="per-object")
        plan["gtex_release"] = release
        plan["gtex_objects"] = len(objects)
    if vep:
        cache_dir = root / "vep"
        cache_dir.mkdir(parents=True, exist_ok=True)
        run([*tools.require("vep_install"), "-a", "cf", "-s", config.annotation.species, "-y", "GRCh38", "-c", str(cache_dir), "--NO_UPDATE"], log_file=cache_dir / "install.log")
        manifest.record("vep/cache", version=tools.version("vep"), source="Ensembl VEP installer", local_path=cache_dir)
    if splice:
        for name in ("spliceai", "pangolin"):
            if not tools.available(name):
                raise FileNotFoundError(f"{name} environment is missing; run scripts/install.sh --full")
            manifest.record(f"software/{name}", version=tools.version(name), source="pinned conda environment", local_path="PATH")
        pangolin_dir = root / "splicing" / "pangolin"
        annotation_db = pangolin_dir / "gencode.annotation.db"
        if not annotation_db.exists():
            pangolin_dir.mkdir(parents=True, exist_ok=True)
            pangolin_command = tools.require("pangolin")
            python_command = (
                [*pangolin_command[:-1], "python"]
                if len(pangolin_command) >= 5 and pangolin_command[-1] == "pangolin"
                else [sys.executable]
            )
            script = (
                "import gffutils,sys; "
                "gffutils.create_db(sys.argv[1], dbfn=sys.argv[2], force=True, keep_order=True, "
                "merge_strategy='merge', sort_attribute_values=True)"
            )
            run(
                [*python_command, "-c", script, str(root / "genome" / "GRCh38" / "gencode.v50.annotation.gtf.gz"), str(annotation_db)],
                log_file=pangolin_dir / "create-db.log",
            )
        manifest.record("splicing/pangolin/annotation_db", version="GENCODE-50", source=GTF_URL, local_path=annotation_db, size=annotation_db.stat().st_size)
    return plan
