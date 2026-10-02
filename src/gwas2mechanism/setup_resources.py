"""Parallel, restartable resource setup used by ``gwas2m setup``."""

from __future__ import annotations

import datetime as dt
import gzip
import os
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import httpx

from gwas2mechanism.config import Config
from gwas2mechanism.constants import SUPERPOPULATIONS
from gwas2mechanism.provenance import ResourceManifest
from gwas2mechanism.qtl.gtex import discover_release
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.resource_status import human_bytes, inspect_resources, validate_resources
from gwas2mechanism.utils.download import Downloader, DownloadOptions, DownloadStats
from gwas2mechanism.utils.io import atomic_write_json
from gwas2mechanism.utils.locking import resource_lock
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.progress import progress
from gwas2mechanism.utils.runtime import resolve_threads
from gwas2mechanism.utils.tools import ToolResolver

FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"


def _progress(message: str) -> None:
    print(f"[gwas2m setup] {message}", file=sys.stderr, flush=True)


def _gunzip_cached(source: Path, destination: Path) -> Path:
    """Atomically unpack an ordinary gzip file, reusing a completed output."""
    with resource_lock(destination):
        if destination.is_file() and destination.stat().st_size > 0 and destination.stat().st_mtime >= source.stat().st_mtime:
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + ".part")
        try:
            total = source.stat().st_size
            with gzip.open(source, "rb") as compressed, temporary.open("wb") as uncompressed, progress(total=total, unit="B", unit_scale=True, desc=f"Extract {source.name}") as bar:
                while block := compressed.read(8 * 1024 * 1024):
                    uncompressed.write(block)
                    bar.update(min(len(block), max(0, total - bar.n)))
            if temporary.stat().st_size == 0:
                raise OSError(f"Decompression produced an empty file: {source}")
            os.replace(temporary, destination)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    return destination


def _safe_extract(archive: Path, destination: Path) -> None:
    """Extract an archive while rejecting path traversal and unsafe links."""
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive, "r:*") as handle:
        members = handle.getmembers()
        for member in members:
            target = (destination / member.name).resolve()
            if (root != target and root not in target.parents) or member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
                raise ValueError(f"Unsafe archive member: {member.name}")
        with progress(members, total=len(members), desc=f"Extract {archive.name}", unit="file") as items:
            for member in items:
                handle.extract(member, destination)


def _gtex_objects(config: Config, release: str) -> list[dict[str, Any]]:
    params: dict[str, str] = {"prefix": f"{release}/"}
    objects: list[dict[str, Any]] = []
    while True:
        response = httpx.get(config.qtl.bucket_listing_url, params=params, timeout=60)
        response.raise_for_status()
        payload = response.json()
        objects.extend(payload.get("items", []))
        token = payload.get("nextPageToken")
        if not token:
            break
        params["pageToken"] = token
    wanted = []
    for item in objects:
        lower = str(item.get("name", "")).lower()
        if "susie" in lower and any(kind in lower for kind in ("eqtl", "sqtl")) and lower.endswith((".tar", ".tar.gz", ".tgz", ".parquet")):
            wanted.append(item)
    return wanted


def _components(config: Config, stats: DownloadStats | None = None) -> tuple[Path, ResourceManifest, ToolResolver, DownloadOptions]:
    root = config.resources_path.resolve()
    root.mkdir(parents=True, exist_ok=True)
    return root, ResourceManifest(root / "manifest.tsv"), ToolResolver(config.tools), DownloadOptions.from_config(config.performance.downloads)


def setup_genome(config: Config, stats: DownloadStats | None = None) -> dict[str, Any]:
    root, manifest, tools, options = _components(config, stats)
    genome = root / "genome" / "GRCh38"
    fasta_gz = genome / "GRCh38.primary_assembly.fa.gz"
    fasta = genome / "GRCh38.primary_assembly.fa"
    gtf_gz = genome / "gencode.v50.annotation.gtf.gz"
    gtf = genome / "gencode.v50.annotation.gtf"
    downloader = Downloader(use_aria2=tools.command("aria2c"), options=options, stats=stats)
    fasta_result = downloader.fetch(FASTA_URL, fasta_gz)
    gtf_result = downloader.fetch(GTF_URL, gtf_gz)
    _gunzip_cached(fasta_gz, fasta)
    _gunzip_cached(gtf_gz, gtf)
    manifest.record("genome/GRCh38/fasta", version="GENCODE-50", source=FASTA_URL, local_path=fasta, size=fasta.stat().st_size, checksum=fasta_result.md5, actual_size=fasta.stat().st_size)
    manifest.record("genome/GRCh38/gtf", version="GENCODE-50", source=GTF_URL, local_path=gtf, size=gtf.stat().st_size, checksum=gtf_result.md5, actual_size=gtf.stat().st_size)
    if tools.available("samtools"):
        fai = Path(f"{fasta}.fai")
        with resource_lock(fai):
            if not fai.exists() or fai.stat().st_mtime < fasta.stat().st_mtime:
                run([*tools.require("samtools"), "faidx", str(fasta)], log_file=genome / "samtools-faidx.log")
    return {"stage": "genome", "status": "complete"}


def setup_reference(config: Config, populations: list[str], stats: DownloadStats | None = None) -> dict[str, Any]:
    root, manifest, tools, options = _components(config, stats)
    threads = resolve_threads(config.performance.threads)
    panel = ReferencePanel(config.reference, root, tools, download_options=options, download_stats=stats)
    panel.download(config.genome.chromosomes, manifest, workers=min(config.performance.downloads.concurrent_files, config.performance.setup.chromosome_workers))
    panel.prepare_populations(populations, config.genome.chromosomes, threads, manifest, chromosome_workers=min(config.performance.setup.chromosome_workers, threads), population_workers=min(config.performance.setup.population_workers, threads))
    return {"stage": "reference", "status": "complete", "populations": populations}


def setup_gtex(config: Config, stats: DownloadStats | None = None) -> dict[str, Any]:
    root, manifest, tools, options = _components(config, stats)
    release = discover_release(config.qtl.bucket_listing_url) if config.qtl.release == "latest" else config.qtl.release.lower()
    objects = _gtex_objects(config, release)
    if not objects:
        raise RuntimeError(f"No compact SuSiE eQTL/sQTL archives found for Adult GTEx {release}")
    destination = root / "gtex" / release
    downloader = Downloader(use_aria2=tools.command("aria2c"), options=options, stats=stats)

    def fetch_object(item: dict[str, Any]) -> int:
        name = str(item["name"])
        path = destination / Path(name).name
        url = f"{config.qtl.bucket_download_url.rstrip('/')}/{name}"
        result = downloader.fetch(url, path, expected_size=int(item["size"]) if item.get("size") else None, md5=item.get("md5Hash"))
        if path.name.endswith((".tar", ".tar.gz", ".tgz")):
            stem = path.name.removesuffix(".gz").removesuffix(".tar").removesuffix(".tgz")
            extracted = destination / "extracted" / stem
            marker = extracted / ".complete.json"
            with resource_lock(marker):
                if not marker.exists():
                    _safe_extract(path, extracted)
                    atomic_write_json(marker, {"status": "COMPLETE", "archive": path.name, "completed_at": dt.datetime.now().isoformat(timespec="seconds")})
        manifest.record(f"gtex/{release}/{path.name}", version=release, source=url, local_path=path, size=result.size, expected_size=int(item["size"]) if item.get("size") else None, actual_size=result.size, checksum=item.get("md5Hash"))
        return result.size

    workers = min(config.performance.setup.gtex_workers, config.performance.downloads.concurrent_files, len(objects))
    total = 0
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = [pool.submit(fetch_object, item) for item in objects]
        with progress(total=len(futures), desc=f"GTEx {release.upper()}", unit="file") as bar:
            for future in as_completed(futures):
                total += future.result()
                bar.update(1)
    manifest.record(f"gtex/{release}/susie", version=release, source=config.qtl.bucket_download_url, local_path=destination, size=total, checksum="per-object")
    return {"stage": "gtex", "status": "complete", "gtex_release": release, "gtex_objects": len(objects)}


def setup_vep(config: Config, stats: DownloadStats | None = None) -> dict[str, Any]:
    root, manifest, tools, _ = _components(config, stats)
    cache_dir = root / "vep"
    marker = cache_dir / ".complete.json"
    with resource_lock(marker):
        cache_files = [path for path in cache_dir.rglob("*") if path.is_file() and path.name not in {"install.log", marker.name}] if cache_dir.exists() else []
        if not marker.exists() or not cache_files:
            cache_dir.mkdir(parents=True, exist_ok=True)
            command = [*tools.require("vep_install"), "-a", "cf", "-s", config.annotation.species, "-y", "GRCh38", "-c", str(cache_dir), "--NO_UPDATE"]
            run(command, log_file=cache_dir / "install.log")
            atomic_write_json(marker, {"status": "COMPLETE", "version": tools.version("vep"), "command": command})
    manifest.record("vep/cache", version=tools.version("vep"), source="Ensembl VEP installer", local_path=cache_dir, software_version=tools.version("vep"))
    return {"stage": "vep", "status": "complete"}


def setup_splice(config: Config, stats: DownloadStats | None = None) -> dict[str, Any]:
    root, manifest, tools, _ = _components(config, stats)
    for name in ("spliceai", "pangolin"):
        if not tools.available(name):
            raise FileNotFoundError(f"{name} environment is missing; run bash Install.sh --full")
        manifest.record(f"software/{name}", version=tools.version(name), source="pinned conda environment", local_path="PATH", software_version=tools.version(name))
    pangolin_dir = root / "splicing" / "pangolin"
    annotation_db = pangolin_dir / "gencode.annotation.db"
    with resource_lock(annotation_db):
        if not annotation_db.exists() or annotation_db.stat().st_size == 0:
            pangolin_dir.mkdir(parents=True, exist_ok=True)
            pangolin_command = tools.require("pangolin")
            python_command = [*pangolin_command[:-1], "python"] if len(pangolin_command) >= 5 and pangolin_command[-1] == "pangolin" else [sys.executable]
            script = "import gffutils,sys; gffutils.create_db(sys.argv[1], dbfn=sys.argv[2], force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)"
            run([*python_command, "-c", script, str(root / "genome" / "GRCh38" / "gencode.v50.annotation.gtf"), str(annotation_db)], log_file=pangolin_dir / "create-db.log")
    manifest.record("splicing/pangolin/annotation_db", version="GENCODE-50", source=GTF_URL, local_path=annotation_db, size=annotation_db.stat().st_size)
    return {"stage": "splice", "status": "complete"}


def run_setup_task(config: Config, task: str, *, chromosome: str | None = None, population: str | None = None, populations: list[str] | None = None, reference: bool = False, gtex: bool = False, vep: bool = False, splice: bool = False) -> dict[str, Any]:
    """Execute one idempotent unit used by generated SLURM scripts."""
    root, manifest, tools, options = _components(config)
    threads = resolve_threads(config.performance.threads)
    if task == "environment-check":
        for name in ("plink2", "samtools", "rscript"):
            tools.require(name)
        return {"task": task, "status": "complete"}
    if task == "genome":
        return setup_genome(config)
    if task == "reference-download":
        if not chromosome:
            raise ValueError("reference-download requires --chromosome")
        panel = ReferencePanel(config.reference, root, tools, download_options=options)
        panel.download_common()
        panel.download_chromosome(chromosome, panel._manifest_checks(), manifest)
        return {"task": task, "chromosome": chromosome, "status": "complete"}
    if task == "reference-convert":
        if not chromosome:
            raise ValueError("reference-convert requires --chromosome")
        ReferencePanel(config.reference, root, tools, download_options=options).convert_all_chromosome(chromosome, threads, manifest)
        return {"task": task, "chromosome": chromosome, "status": "complete"}
    if task == "reference-population":
        if not chromosome or not population:
            raise ValueError("reference-population requires --chromosome and --population")
        ReferencePanel(config.reference, root, tools, download_options=options).subset_population_chromosome(population, chromosome, threads, manifest)
        return {"task": task, "chromosome": chromosome, "population": population, "status": "complete"}
    if task == "gtex":
        return setup_gtex(config)
    if task == "vep":
        return setup_vep(config)
    if task == "splice":
        return setup_splice(config)
    if task == "verify":
        errors = validate_resources(config, reference=reference, gtex=gtex, vep=vep, splice=splice, populations=populations)
        if errors:
            raise RuntimeError("Resource validation failed:\n- " + "\n- ".join(errors))
        return {"task": task, "status": "complete"}
    raise ValueError(f"Unknown setup task: {task}")


def setup_resources(config: Config, *, reference: bool, gtex: bool, vep: bool, splice: bool, populations: list[str], full: bool, dry_run: bool) -> dict[str, Any]:
    invalid = set(populations) - set(SUPERPOPULATIONS)
    if invalid:
        raise ValueError(f"No 1000 Genomes super-population reference for {sorted(invalid)}")
    root = config.resources_path.resolve()
    selected = {"reference": reference, "gtex": gtex, "vep": vep, "splice": splice}
    plan: dict[str, Any] = {"cache": str(root), "date": dt.date.today().isoformat(), **selected, "populations": populations, "full": full, "dry_run": dry_run, "executor": "local"}
    if dry_run:
        return plan
    _progress(f"Resource cache: {root}")
    started = time.monotonic()
    stats = DownloadStats()
    needs_genome = reference or splice
    independent: list[tuple[str, Any]] = []
    if needs_genome:
        independent.append(("genome", lambda: setup_genome(config, stats)))
    if reference:
        independent.append(("reference", lambda: setup_reference(config, populations, stats)))
    if gtex:
        independent.append(("gtex", lambda: setup_gtex(config, stats)))
    if vep:
        independent.append(("vep", lambda: setup_vep(config, stats)))
    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, min(4, len(independent)))) as pool:
        futures = {pool.submit(call): name for name, call in independent}
        total_stages = len(futures) + int(splice)
        with progress(total=total_stages, desc="Resource setup", unit="stage") as bar:
            for future in as_completed(futures):
                results.append(future.result())
                bar.update(1)
            if splice:
                results.append(setup_splice(config, stats))
                bar.update(1)
    status = inspect_resources(config)
    elapsed = time.monotonic() - started
    performance = {"total_elapsed_seconds": round(elapsed, 2), **stats.as_dict(), "bytes_downloaded_human": human_bytes(stats.bytes_downloaded), "resources_complete": status["complete"], "resources_partial": status["partial"], "resources_missing": status["missing"]}
    plan.update({"stages": results, "performance": performance})
    _progress(f"Resource cache: {root}")
    _progress(f"Completed in {elapsed:.1f}s; downloaded {stats.files_downloaded} files ({human_bytes(stats.bytes_downloaded)}), cache hits {stats.cache_hits}, resumed {stats.resumed}")
    return plan
