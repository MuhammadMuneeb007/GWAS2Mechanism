"""Signed LD matrices per population.

* ``plink2`` engine: ``--r-unphased square bin4 ref-based`` (legacy Step04
  logic) read with a memory map; reused when size-verified.
* ``dosage`` engine: numpy correlation of genotype dosages stored in
  ``chr{c}.dosage.npz`` (small synthetic panels / tests).

Both return correlations for one consistent allele coding per variant. Because
``r(-x, -y) = r(x, y)``, a REF-based matrix equals an ALT-based matrix, so it is
valid with ALT-oriented z-scores. LD matrices are **never** pooled across
populations.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.runtime import thread_env

log = logging.getLogger(__name__)


@dataclass
class LDResult:
    ids: list[str]
    R: np.ndarray  # float64, symmetric, unit diagonal after stabilisation
    diagnostics: dict[str, float]


def stabilize(R: np.ndarray, ridge: float) -> tuple[np.ndarray, dict[str, float]]:
    """Symmetrise, unit diagonal, optional shrinkage ``R <- (1-ridge) R + ridge I``."""
    R = np.asarray(R, dtype=np.float64)
    if not np.all(np.isfinite(R)):
        raise ValueError("LD matrix contains non-finite values")
    sym_err = float(np.max(np.abs(R - R.T))) if R.size else 0.0
    R = (R + R.T) / 2.0
    np.fill_diagonal(R, 1.0)
    if ridge and ridge > 0:
        R *= 1.0 - ridge
        np.fill_diagonal(R, 1.0)
    try:
        min_eig = float(np.linalg.eigvalsh(R).min()) if R.shape[0] <= 5000 else float("nan")
    except np.linalg.LinAlgError:
        min_eig = float("nan")
    return R, {"symmetry_error": sym_err, "ld_ridge": ridge, "min_eigenvalue": min_eig}


class LDProvider:
    def __init__(self, panel: ReferencePanel, engine: str, maf: float, geno: float, plink2: list[str] | None, threads: int = 1):
        self.panel = panel
        self.engine = engine
        self.maf = maf
        self.geno = geno
        self.plink2 = plink2
        self.threads = threads

    def compute(self, population: str, chrom: str, ids: list[str], workdir: Path, ridge: float, max_variants: int) -> LDResult:
        if self.engine == "dosage":
            ids_kept, R = self._dosage(population, chrom, ids)
        else:
            ids_kept, R = self._plink2(population, chrom, ids, workdir)
        if len(ids_kept) > max_variants:
            raise ValueError(f"{len(ids_kept)} LD variants exceeds max_ld_variants={max_variants}")
        R, diag = stabilize(R, ridge)
        diag["n_ld_variants"] = float(len(ids_kept))
        diag["n_requested"] = float(len(ids))
        return LDResult(ids_kept, R, diag)

    # -- plink2 ------------------------------------------------------------
    def _plink2(self, population: str, chrom: str, ids: list[str], workdir: Path) -> tuple[list[str], np.ndarray]:
        if self.plink2 is None:
            raise FileNotFoundError("plink2 is required for LD computation (reference.ld_engine: plink2)")
        workdir.mkdir(parents=True, exist_ok=True)
        prefix = workdir / f"{population}_LD"
        ld_file = Path(f"{prefix}.unphased.vcor1.bin")
        vars_file = Path(f"{prefix}.unphased.vcor1.bin.vars")
        extract = workdir / f"{population}_ld_ids.txt"
        wanted = "\n".join(ids) + "\n"
        reuse = (
            ld_file.exists()
            and vars_file.exists()
            and extract.exists()
            and extract.read_text(encoding="utf-8") == wanted
        )
        if reuse:
            kept = [x for x in vars_file.read_text(encoding="utf-8").split() if x]
            if ld_file.stat().st_size == len(kept) ** 2 * 4:
                return kept, self._read_bin4(ld_file, len(kept))
        extract.write_text(wanted, encoding="utf-8")
        run(
            [
                *self.plink2,
                "--pfile", str(self.panel.prefix(population, chrom)),
                "--extract", str(extract),
                "--maf", str(self.maf),
                "--geno", str(self.geno),
                "--min-alleles", "2",
                "--max-alleles", "2",
                "--r-unphased", "square", "bin4", "ref-based", "yes-really",
                "--threads", str(self.threads),
                "--out", str(prefix),
            ],
            log_file=workdir / f"{population}_LD.log",
            env=thread_env(self.threads),
        )
        kept = [x for x in vars_file.read_text(encoding="utf-8").split() if x]
        expected = len(kept) ** 2 * 4
        if ld_file.stat().st_size != expected:
            raise ValueError(f"LD matrix size mismatch: {ld_file.stat().st_size} != {expected}")
        return kept, self._read_bin4(ld_file, len(kept))

    @staticmethod
    def _read_bin4(path: Path, p: int) -> np.ndarray:
        mm = np.memmap(path, dtype="<f4", mode="r", shape=(p, p))
        return np.array(mm, dtype=np.float64)

    # -- dosage (numpy) -------------------------------------------------------
    def _dosage(self, population: str, chrom: str, ids: list[str]) -> tuple[list[str], np.ndarray]:
        data = np.load(self.panel.dosage_file(population, chrom), allow_pickle=False)
        all_ids = data["ids"].astype(str)
        index = {v: i for i, v in enumerate(all_ids)}
        rows = [index[v] for v in ids if v in index]
        G = data["dosage"][rows].astype(np.float64)  # (variants, samples) ALT counts
        missing = np.isnan(G).mean(axis=1) if np.isnan(G).any() else np.zeros(len(rows))
        freq = np.nanmean(G, axis=1) / 2.0
        maf = np.minimum(freq, 1.0 - freq)
        keep = (maf >= self.maf) & (missing <= self.geno)
        G = G[keep]
        kept_ids = [all_ids[r] for r, k in zip(rows, keep, strict=True) if k]
        if G.shape[0] < 2:
            return kept_ids, np.eye(G.shape[0])
        G = np.where(np.isnan(G), np.nanmean(G, axis=1, keepdims=True), G)
        R = np.corrcoef(G)
        return kept_ids, R
