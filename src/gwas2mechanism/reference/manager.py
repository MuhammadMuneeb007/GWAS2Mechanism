"""Restartable ancestry-specific 1000 Genomes GRCh38 reference preparation.

Each raw chromosome is converted once to an ALL-sample PGEN intermediate.
EUR/AFR/EAS/SAS/AMR are then independently subset from that binary dataset;
ALL is never exposed as a pooled cross-ancestry LD panel.
"""

from __future__ import annotations

import datetime as dt
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from gwas2mechanism.config import ReferenceConfig
from gwas2mechanism.constants import SUPERPOPULATIONS
from gwas2mechanism.provenance import ResourceManifest
from gwas2mechanism.utils.download import Downloader, DownloadOptions, DownloadStats, fetch_range
from gwas2mechanism.utils.io import atomic_write_json, read_json, write_parquet
from gwas2mechanism.utils.locking import resource_lock
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.progress import progress
from gwas2mechanism.utils.runtime import thread_env
from gwas2mechanism.utils.tools import ToolResolver

log = logging.getLogger(__name__)

PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"
MANIFEST_FILE = "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"
README_FILE = "20190312_biallelic_SNV_and_INDEL_README.txt"


class ReferencePanel:
    def __init__(self, config: ReferenceConfig, resources_dir: Path, tools: ToolResolver, *, download_options: DownloadOptions | object | None = None, download_stats: DownloadStats | None = None):
        self.config = config
        self.root = resources_dir / "reference" / config.panel
        self.raw = self.root / "raw"
        self.all_dir = self.root / "ALL"
        self.tools = tools
        self.download_options = download_options
        self.download_stats = download_stats

    def pop_dir(self, population: str) -> Path:
        if population not in SUPERPOPULATIONS:
            raise ValueError(f"No single-ancestry reference for {population!r}; pooled panels are not supported")
        return self.root / population

    def prefix(self, population: str, chrom: str) -> Path:
        return self.pop_dir(population) / f"chr{chrom}"

    def all_prefix(self, chrom: str) -> Path:
        return self.all_dir / f"chr{chrom}"

    def pvar_parquet(self, population: str, chrom: str) -> Path:
        return self.pop_dir(population) / f"chr{chrom}.pvar.parquet"

    def all_pvar_parquet(self, chrom: str) -> Path:
        return self.all_dir / f"chr{chrom}.pvar.parquet"

    def dosage_file(self, population: str, chrom: str) -> Path:
        return self.pop_dir(population) / f"chr{chrom}.dosage.npz"

    def chromosome_marker(self, directory: Path, chrom: str) -> Path:
        return directory / f"chr{chrom}.complete.json"

    def sample_size(self, population: str) -> int:
        samples = self.pop_dir(population) / "samples.txt"
        if samples.exists():
            return sum(1 for line in samples.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#"))
        return 0

    def is_ready(self, population: str, chromosomes: list[str]) -> bool:
        return all(self.population_chromosome_complete(population, chrom) for chrom in chromosomes)

    def vcf_name(self, chrom: str) -> str:
        return self.config.vcf_template.format(chrom=chrom)

    def raw_chromosome_complete(self, chrom: str) -> bool:
        vcf = self.raw / self.vcf_name(chrom)
        index = self.raw / f"{self.vcf_name(chrom)}.tbi"
        return all(path.is_file() and path.stat().st_size > 0 and path.with_name(path.name + ".verified.json").is_file() for path in (vcf, index))

    def _dataset_complete(self, prefix: Path, parquet: Path, marker: Path) -> bool:
        outputs = [Path(f"{prefix}.{suffix}") for suffix in ("pgen", "pvar", "psam")] + [parquet]
        if not marker.is_file() or not all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            return False
        try:
            data = read_json(marker)
            sizes = data.get("output_sizes", {})
            return data.get("status") == "COMPLETE" and int(data.get("n_samples", 0)) > 0 and int(data.get("n_variants", 0)) > 0 and all(int(sizes.get(path.name, -1)) == path.stat().st_size for path in outputs)
        except (OSError, ValueError, TypeError):
            return False

    def all_chromosome_complete(self, chrom: str) -> bool:
        return self._dataset_complete(self.all_prefix(chrom), self.all_pvar_parquet(chrom), self.chromosome_marker(self.all_dir, chrom))

    def population_chromosome_complete(self, population: str, chrom: str) -> bool:
        return self._dataset_complete(self.prefix(population, chrom), self.pvar_parquet(population, chrom), self.chromosome_marker(self.pop_dir(population), chrom))

    def _manifest_checks(self) -> dict[str, tuple[int, str]]:
        url = f"{self.config.base_url}/{MANIFEST_FILE}"
        path = self.raw / MANIFEST_FILE
        try:
            Downloader(options=self.download_options, stats=self.download_stats).fetch(url, path)
        except OSError as exc:
            log.warning("Could not fetch release MANIFEST (%s); size-only verification", exc)
            return {}
        checks: dict[str, tuple[int, str]] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and parts[1].isdigit():
                checks[parts[0].lstrip("./")] = (int(parts[1]), parts[2].strip())
        return checks

    def _downloader(self, large: bool = False) -> Downloader:
        aria2 = self.tools.command("aria2c") if large else None
        return Downloader(use_aria2=aria2, options=self.download_options, stats=self.download_stats)

    def download_common(self) -> None:
        self.raw.mkdir(parents=True, exist_ok=True)
        self._downloader().fetch(self.config.panel_url, self.raw / PANEL_FILE)
        self._downloader().fetch(f"{self.config.base_url}/{README_FILE}", self.raw / README_FILE)

    def download_chromosome(self, chrom: str, checks: dict[str, tuple[int, str]] | None = None, manifest: ResourceManifest | None = None) -> None:
        checks = checks or {}
        started = dt.datetime.now().isoformat(timespec="seconds")
        total = 0
        for name in (self.vcf_name(chrom), f"{self.vcf_name(chrom)}.tbi"):
            size, md5 = checks.get(name, (None, None))
            result = self._downloader(large=(size or 0) > 50_000_000).fetch(f"{self.config.base_url}/{name}", self.raw / name, expected_size=size, md5=md5)
            total += result.size
        if manifest:
            manifest.record(f"reference/{self.config.panel}/raw/chr{chrom}", version="20190312_biallelic_SNV_and_INDEL", source=self.config.base_url, local_path=self.raw / self.vcf_name(chrom), size=total, checksum="md5:MANIFEST" if checks else "verified-size", chromosome=chrom, started_at=started)

    def download(self, chromosomes: list[str], manifest: ResourceManifest, workers: int = 4) -> None:
        self.download_common()
        checks = self._manifest_checks()
        workers = max(1, min(workers, len(chromosomes)))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(self.download_chromosome, chrom, checks, manifest) for chrom in chromosomes]
            with progress(total=len(futures), desc="1000G downloads", unit="chr") as bar:
                for future in as_completed(futures):
                    future.result()
                    bar.update(1)
        manifest.record(f"reference/{self.config.panel}/raw", version="20190312_biallelic_SNV_and_INDEL", source=self.config.base_url, local_path=self.raw, status="OK")

    def write_sample_lists(self, populations: list[str]) -> dict[str, Path]:
        invalid = set(populations) - set(SUPERPOPULATIONS)
        if invalid:
            raise ValueError(f"No single-ancestry reference for {sorted(invalid)}")
        result: dict[str, Path] = {}
        with resource_lock(self.root / "sample-lists"):
            panel_path = self.raw / PANEL_FILE
            if not panel_path.exists():
                raise FileNotFoundError(f"{panel_path} missing - download the 1000 Genomes sample panel first")
            outputs = {population: self.pop_dir(population) / "samples.txt" for population in populations}
            missing = [population for population, path in outputs.items() if not path.is_file() or self._sample_list_count(path) <= 0]
            panel = None
            if missing:
                panel = pl.read_csv(panel_path, separator="\t", infer_schema=False, truncate_ragged_lines=True)
                panel = panel.rename({column: column.strip() for column in panel.columns})
            for population, out in outputs.items():
                if population not in missing:
                    result[population] = out
                    continue
                assert panel is not None
                samples = panel.filter(pl.col("super_pop") == population)["sample"].to_list()
                if not samples:
                    raise RuntimeError(f"No 1000 Genomes samples for super-population {population}")
                expected = "#IID\n" + "\n".join(samples) + "\n"
                if not out.exists() or out.read_text(encoding="utf-8") != expected:
                    out.parent.mkdir(parents=True, exist_ok=True)
                    temporary = out.with_name(out.name + ".part")
                    temporary.write_text(expected, encoding="utf-8")
                    temporary.replace(out)
                result[population] = out
        return result

    def _sample_list_count(self, path: Path) -> int:
        return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#"))

    def write_sample_list(self, population: str) -> Path:
        return self.write_sample_lists([population])[population]

    def _cache_pvar_prefix(self, prefix: Path, destination: Path) -> int:
        frame = pl.scan_csv(Path(f"{prefix}.pvar"), separator="\t", comment_prefix="##", infer_schema=False)
        cols = frame.collect_schema().names()
        selected = frame.select(pl.col(cols[0]).str.replace(r"(?i)^chr", "").alias("CHR"), pl.col("POS").cast(pl.Int64), pl.col("ID"), pl.col("REF").str.to_uppercase(), pl.col("ALT").str.to_uppercase())
        n_variants = int(selected.select(pl.len()).collect().item())
        write_parquet(selected, destination)
        return n_variants

    def _sample_count(self, psam: Path) -> int:
        return sum(1 for line in psam.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#"))

    def _write_dataset_marker(self, prefix: Path, parquet: Path, marker: Path, *, chrom: str, ancestry: str, command: list[str], n_variants: int) -> None:
        outputs = [Path(f"{prefix}.{suffix}") for suffix in ("pgen", "pvar", "psam")] + [parquet]
        if not all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            raise RuntimeError(f"PLINK2 did not create every expected output for {prefix}")
        n_samples = self._sample_count(Path(f"{prefix}.psam"))
        if n_samples <= 0 or n_variants <= 0:
            raise RuntimeError(f"Invalid empty PGEN dataset: {prefix}")
        atomic_write_json(marker, {"status": "COMPLETE", "chromosome": chrom, "ancestry": ancestry, "n_samples": n_samples, "n_variants": n_variants, "tool_version": self.tools.version("plink2"), "command": command, "completed_at": dt.datetime.now().isoformat(timespec="seconds"), "output_sizes": {path.name: path.stat().st_size for path in outputs}})

    def convert_all_chromosome(self, chrom: str, threads: int, manifest: ResourceManifest | None = None) -> None:
        prefix = self.all_prefix(chrom)
        parquet = self.all_pvar_parquet(chrom)
        marker = self.chromosome_marker(self.all_dir, chrom)
        with resource_lock(marker):
            if self.all_chromosome_complete(chrom):
                return
            vcf = self.raw / self.vcf_name(chrom)
            if not self.raw_chromosome_complete(chrom):
                raise FileNotFoundError(f"Verified VCF/index missing for chromosome {chrom}")
            prefix.parent.mkdir(parents=True, exist_ok=True)
            command = [*self.tools.require("plink2"), "--vcf", str(vcf), "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "1000", "missing", "--rm-dup", "force-first", "--max-alleles", "2", "--make-pgen", "--threads", str(threads), "--out", str(prefix)]
            started = dt.datetime.now().isoformat(timespec="seconds")
            run(command, log_file=Path(f"{prefix}.convert.log"), env=thread_env(threads))
            n_variants = self._cache_pvar_prefix(prefix, parquet)
            self._write_dataset_marker(prefix, parquet, marker, chrom=chrom, ancestry="ALL_INTERMEDIATE", command=command, n_variants=n_variants)
            if manifest:
                manifest.record(f"reference/{self.config.panel}/ALL/chr{chrom}", version=self.config.panel, source=f"derived:{vcf.name}", local_path=prefix.parent, chromosome=chrom, ancestry="ALL_INTERMEDIATE", started_at=started, software_version=self.tools.version("plink2"), preparation_command=" ".join(command))

    def subset_population_chromosome(self, population: str, chrom: str, threads: int, manifest: ResourceManifest | None = None) -> None:
        if population not in SUPERPOPULATIONS:
            raise ValueError("A pooled cross-ancestry LD panel is not supported")
        if not self.all_chromosome_complete(chrom):
            raise FileNotFoundError(f"ALL PGEN intermediate for chromosome {chrom} is incomplete")
        keep = self.write_sample_list(population)
        prefix = self.prefix(population, chrom)
        parquet = self.pvar_parquet(population, chrom)
        marker = self.chromosome_marker(self.pop_dir(population), chrom)
        with resource_lock(marker):
            if self.population_chromosome_complete(population, chrom):
                return
            prefix.parent.mkdir(parents=True, exist_ok=True)
            command = [*self.tools.require("plink2"), "--pfile", str(self.all_prefix(chrom)), "--keep", str(keep), "--make-pgen", "--threads", str(threads), "--out", str(prefix)]
            started = dt.datetime.now().isoformat(timespec="seconds")
            run(command, log_file=Path(f"{prefix}.prepare.log"), env=thread_env(threads))
            n_variants = self._cache_pvar_prefix(prefix, parquet)
            self._write_dataset_marker(prefix, parquet, marker, chrom=chrom, ancestry=population, command=command, n_variants=n_variants)
            if manifest:
                manifest.record(f"reference/{self.config.panel}/{population}/chr{chrom}", version=self.config.panel, source="derived: ALL PGEN --keep super_pop samples", local_path=prefix.parent, chromosome=chrom, ancestry=population, started_at=started, software_version=self.tools.version("plink2"), preparation_command=" ".join(command))

    def prepare_populations(self, populations: list[str], chromosomes: list[str], threads: int, manifest: ResourceManifest, *, chromosome_workers: int = 4, population_workers: int = 8) -> None:
        self.write_sample_lists(populations)
        conversion_workers = max(1, min(chromosome_workers, len(chromosomes)))
        per_conversion = max(1, threads // conversion_workers)
        with ThreadPoolExecutor(max_workers=conversion_workers) as pool:
            futures = [pool.submit(self.convert_all_chromosome, chrom, per_conversion, manifest) for chrom in chromosomes]
            with progress(total=len(futures), desc="1000G PGEN conversion", unit="chr") as bar:
                for future in as_completed(futures):
                    future.result()
                    bar.update(1)
        tasks = [(population, chrom) for population in populations for chrom in chromosomes]
        subset_workers = max(1, min(population_workers, len(tasks), threads))
        per_subset = max(1, threads // subset_workers)
        bars = {population: progress(total=len(chromosomes), desc=f"1000G {population}", unit="chr") for population in populations}
        try:
            with ThreadPoolExecutor(max_workers=subset_workers) as pool:
                future_map = {pool.submit(self.subset_population_chromosome, population, chrom, per_subset, manifest): population for population, chrom in tasks}
                for future in as_completed(future_map):
                    future.result()
                    bars[future_map[future]].update(1)
        finally:
            for bar in bars.values():
                bar.close()
        for population in populations:
            manifest.record(f"reference/{self.config.panel}/{population}", version=self.config.panel, source="derived: ALL PGEN --keep super_pop samples", local_path=self.pop_dir(population), ancestry=population)

    def prepare(self, population: str, chromosomes: list[str], threads: int, manifest: ResourceManifest) -> None:
        self.prepare_populations([population], chromosomes, threads, manifest, chromosome_workers=min(4, threads), population_workers=min(8, threads))

    def cache_pvar(self, population: str, chrom: str) -> Path:
        self._cache_pvar_prefix(self.prefix(population, chrom), self.pvar_parquet(population, chrom))
        return self.pvar_parquet(population, chrom)

    def region_variants(self, population: str, chrom: str, start: int, end: int) -> pl.DataFrame:
        path = self.pvar_parquet(population, chrom)
        if not path.exists():
            if Path(f"{self.prefix(population, chrom)}.pvar").exists():
                self.cache_pvar(population, chrom)
            else:
                raise FileNotFoundError(f"No {population} reference variants for chr{chrom} ({path})")
        return pl.scan_parquet(path).filter(pl.col("POS").is_between(start, end)).collect()


def _chrom_key(c: str) -> tuple[int, str]:
    return (int(c), "") if c.isdigit() else (99, c)


__all__ = ["ReferencePanel", "fetch_range"]
