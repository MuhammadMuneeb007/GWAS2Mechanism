"""Ancestry-specific LD reference panels (1000 Genomes GRCh38, 20190312 release).

Layout (cached once under ``<resources>/reference/<panel>/``)::

    raw/                          chromosome VCF + TBI, sample panel, MANIFEST
    EUR/ AFR/ EAS/ SAS/ AMR/
        samples.txt               super-population sample list
        chr{c}.pgen/.pvar/.psam   PLINK2 per population (one pass: plink2 --keep)
        chr{c}.pvar.parquet       variant table used for vectorised matching
        manifest.tsv, .complete

Variant IDs are set to ``CHR:POS:REF:ALT`` for every variant. There is **no**
pooled multi-population panel: each population keeps its own LD.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

from gwas2mechanism.config import ReferenceConfig
from gwas2mechanism.constants import SUPERPOPULATIONS
from gwas2mechanism.provenance import ResourceManifest
from gwas2mechanism.utils.download import Downloader, fetch_range
from gwas2mechanism.utils.io import atomic_write_json, read_json, write_parquet
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.runtime import thread_env
from gwas2mechanism.utils.tools import ToolResolver

log = logging.getLogger(__name__)

PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"
MANIFEST_FILE = "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"
README_FILE = "20190312_biallelic_SNV_and_INDEL_README.txt"


class ReferencePanel:
    def __init__(self, config: ReferenceConfig, resources_dir: Path, tools: ToolResolver):
        self.config = config
        self.root = resources_dir / "reference" / config.panel
        self.raw = self.root / "raw"
        self.tools = tools

    # ---------------------------------------------------------------- paths
    def pop_dir(self, population: str) -> Path:
        if population not in SUPERPOPULATIONS:
            raise ValueError(f"No single-ancestry reference for {population!r}; pooled panels are not supported")
        return self.root / population

    def prefix(self, population: str, chrom: str) -> Path:
        return self.pop_dir(population) / f"chr{chrom}"

    def pvar_parquet(self, population: str, chrom: str) -> Path:
        return self.pop_dir(population) / f"chr{chrom}.pvar.parquet"

    def dosage_file(self, population: str, chrom: str) -> Path:
        """Small numpy genotype panels (``ld_engine: dosage``; synthetic/test use)."""
        return self.pop_dir(population) / f"chr{chrom}.dosage.npz"

    def sample_size(self, population: str) -> int:
        samples = self.pop_dir(population) / "samples.txt"
        if samples.exists():
            return sum(1 for line in samples.read_text().splitlines() if line.strip() and not line.startswith("#"))
        marker = self.pop_dir(population) / ".complete"
        if marker.exists():
            return int(read_json(marker).get("n_samples", 0))
        return 0

    def is_ready(self, population: str, chromosomes: list[str]) -> bool:
        marker = self.pop_dir(population) / ".complete"
        if not marker.exists():
            return False
        done = set(read_json(marker).get("chromosomes", []))
        return set(chromosomes) <= done

    # ------------------------------------------------------------- download
    def vcf_name(self, chrom: str) -> str:
        return self.config.vcf_template.format(chrom=chrom)

    def _manifest(self) -> dict[str, tuple[int, str]]:
        url = f"{self.config.base_url}/{MANIFEST_FILE}"
        path = self.raw / MANIFEST_FILE
        try:
            Downloader().fetch(url, path)
        except OSError as exc:
            log.warning("Could not fetch release MANIFEST (%s); size-only verification", exc)
            return {}
        out = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and parts[1].isdigit():
                out[parts[0].lstrip("./")] = (int(parts[1]), parts[2].strip())
        return out

    def download(self, chromosomes: list[str], manifest: ResourceManifest, workers: int = 4) -> None:
        self.raw.mkdir(parents=True, exist_ok=True)
        checks = self._manifest()
        aria2 = self.tools.command("aria2c")
        jobs = [(self.config.panel_url, self.raw / PANEL_FILE, None, None)]
        jobs.append((f"{self.config.base_url}/{README_FILE}", self.raw / README_FILE, None, None))
        for chrom in chromosomes:
            for name in (self.vcf_name(chrom), self.vcf_name(chrom) + ".tbi"):
                size, md5 = checks.get(name, (None, None))
                jobs.append((f"{self.config.base_url}/{name}", self.raw / name, size, md5))

        def fetch(job: tuple[str, Path, int | None, str | None]):
            url, dest, size, md5 = job
            return Downloader(use_aria2=aria2 if (size or 0) > 50_000_000 else None).fetch(
                url, dest, md5=md5, expected_size=size
            )

        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(fetch, jobs))
        total = sum(r.size for r in results)
        manifest.record(
            f"reference/{self.config.panel}/raw",
            version="20190312_biallelic_SNV_and_INDEL",
            source=self.config.base_url,
            local_path=self.raw,
            size=total,
            checksum="md5:MANIFEST" if checks else "size-only",
        )

    # ---------------------------------------------------------- preparation
    def write_sample_list(self, population: str) -> Path:
        panel = pl.read_csv(self.raw / PANEL_FILE, separator="\t", infer_schema=False, truncate_ragged_lines=True)
        panel = panel.rename({c: c.strip() for c in panel.columns})
        samples = panel.filter(pl.col("super_pop") == population)
        if samples.is_empty():
            raise RuntimeError(f"No 1000 Genomes samples for super-population {population}")
        out = self.pop_dir(population) / "samples.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("#IID\n" + "\n".join(samples["sample"].to_list()) + "\n", encoding="utf-8")
        return out

    def prepare(self, population: str, chromosomes: list[str], threads: int, manifest: ResourceManifest) -> None:
        if self.is_ready(population, chromosomes):
            log.info("[cache] %s reference ready", population)
            return
        plink2 = self.tools.require("plink2")
        keep = self.write_sample_list(population)
        done: list[str] = []
        marker = self.pop_dir(population) / ".complete"
        if marker.exists():
            done = list(read_json(marker).get("chromosomes", []))
        for chrom in chromosomes:
            prefix = self.prefix(population, chrom)
            if chrom in done and Path(f"{prefix}.pgen").exists():
                continue
            vcf = self.raw / self.vcf_name(chrom)
            if not vcf.exists():
                raise FileNotFoundError(f"{vcf} missing - run `gwas2m setup --reference` first")
            run(
                [
                    *plink2,
                    "--vcf", str(vcf),
                    "--keep", str(keep),
                    "--set-all-var-ids", "@:#:$r:$a",
                    "--new-id-max-allele-len", "1000", "missing",
                    "--rm-dup", "force-first",
                    "--max-alleles", "2",
                    "--make-pgen",
                    "--threads", str(threads),
                    "--out", str(prefix),
                ],
                log_file=Path(f"{prefix}.prepare.log"),
                env=thread_env(threads),
            )
            self.cache_pvar(population, chrom)
            done.append(chrom)
            atomic_write_json(marker, {"chromosomes": sorted(set(done), key=_chrom_key), "n_samples": self.sample_size(population)})
        manifest.record(
            f"reference/{self.config.panel}/{population}",
            version=self.config.panel,
            source="derived: plink2 --keep super_pop samples",
            local_path=self.pop_dir(population),
            status="OK",
        )

    def cache_pvar(self, population: str, chrom: str) -> Path:
        pvar = Path(f"{self.prefix(population, chrom)}.pvar")
        frame = pl.scan_csv(pvar, separator="\t", comment_prefix="##", infer_schema=False)
        cols = frame.collect_schema().names()
        chrom_col = cols[0]
        out = frame.select(
            pl.col(chrom_col).str.replace(r"(?i)^chr", "").alias("CHR"),
            pl.col("POS").cast(pl.Int64),
            pl.col("ID"),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
        )
        return write_parquet(out, self.pvar_parquet(population, chrom))

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


def remote_file_size(url: str) -> int | None:  # pragma: no cover - network helper
    try:
        return len(fetch_range(url, 1))
    except OSError:
        return None
