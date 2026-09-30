"""SuSiE-RSS execution and result normalisation.

The R ``susieR`` implementation is the production default.  A compact
single-effect approximation is included only for deterministic smoke tests and
small demonstrations; diagnostics identify the engine so it cannot be confused
with a production susieR result.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import logsumexp

from gwas2mechanism.utils.proc import run
from gwas2mechanism.utils.tools import ToolResolver


@dataclass
class FineMapResult:
    variants: pl.DataFrame
    credible_sets: pl.DataFrame
    summary: pl.DataFrame
    diagnostics: pl.DataFrame


def _credible_set(alpha: np.ndarray, coverage: float) -> tuple[np.ndarray, float]:
    order = np.argsort(-alpha, kind="stable")
    stop = int(np.searchsorted(np.cumsum(alpha[order]), coverage, side="left")) + 1
    keep = order[: min(stop, order.size)]
    return keep, float(alpha[keep].sum())


def _validate_inputs(z: np.ndarray, ld: np.ndarray) -> None:
    if z.ndim != 1 or ld.shape != (z.size, z.size):
        raise ValueError("z must be length M and LD must be M x M")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(ld)):
        raise ValueError("Fine-mapping inputs must be finite")
    if not np.allclose(ld, ld.T, atol=1e-7):
        raise ValueError("LD matrix is not symmetric")


def run_numpy_smoke_susie(
    variants: pl.DataFrame,
    ld: np.ndarray,
    *,
    locus_id: str,
    population: str,
    coverage: float = 0.95,
    max_effects: int = 10,
    prior_variance: float = 0.2,
) -> FineMapResult:
    """Deterministic ABF approximation for tests; never labelled as susieR."""
    if "Z" not in variants.columns:
        raise ValueError("Variant table requires Z")
    z = variants["Z"].to_numpy().astype(np.float64)
    _validate_inputs(z, ld)
    if z.size == 0:
        raise ValueError("Cannot fine-map an empty locus")
    # Wakefield-style log ABF, vectorised.  Repeated effects are selected after
    # LD-residualisation to provide realistic smoke-test credible sets.
    residual = z.copy()
    alphas: list[np.ndarray] = []
    for _ in range(min(max_effects, max(1, z.size))):
        log_abf = 0.5 * ((prior_variance / (1.0 + prior_variance)) * residual**2 - np.log1p(prior_variance))
        alpha = np.exp(log_abf - logsumexp(log_abf))
        alphas.append(alpha)
        lead = int(np.argmax(alpha))
        if abs(residual[lead]) < 1.0:
            break
        residual = residual - ld[:, lead] * residual[lead]
    alpha_matrix = np.vstack(alphas)
    pip = 1.0 - np.prod(1.0 - alpha_matrix, axis=0)
    cs_rows: list[dict[str, object]] = []
    cs_membership: list[list[str]] = [[] for _ in range(z.size)]
    ids = variants["VARIANT_ID"].cast(pl.String).to_list()
    for effect, alpha in enumerate(alphas, start=1):
        keep, achieved = _credible_set(alpha, coverage)
        cs_id = f"{locus_id}:CS{effect}"
        purity = float(np.min(np.abs(ld[np.ix_(keep, keep)]))) if keep.size > 1 else 1.0
        cs_rows.append(
            {
                "LOCUS_ID": locus_id,
                "POPULATION": population,
                "CS_ID": cs_id,
                "EFFECT": effect,
                "COVERAGE": achieved,
                "PURITY_MIN_ABS_COR": purity,
                "N_VARIANTS": int(keep.size),
                "VARIANTS": [ids[i] for i in keep],
            }
        )
        for i in keep:
            cs_membership[int(i)].append(cs_id)
    out = variants.with_columns(
        pl.Series("PIP", pip),
        pl.Series("CREDIBLE_SETS", cs_membership, dtype=pl.List(pl.String)),
        pl.lit(locus_id).alias("LOCUS_ID"),
        pl.lit(population).alias("POPULATION"),
    )
    cs = pl.DataFrame(cs_rows)
    summary = pl.DataFrame(
        [
            {
                "LOCUS_ID": locus_id,
                "POPULATION": population,
                "N_VARIANTS": variants.height,
                "N_EFFECTS": len(alphas),
                "N_CREDIBLE_SETS": cs.height,
                "MAX_PIP": float(np.max(pip)),
                "STATUS": "OK",
            }
        ]
    )
    diagnostics = pl.DataFrame(
        [
            {
                "LOCUS_ID": locus_id,
                "POPULATION": population,
                "ENGINE": "numpy_smoke_abf",
                "PRODUCTION_SUSIE": False,
                "LD_MIN_EIGENVALUE": float(np.linalg.eigvalsh(ld).min()),
                "LD_CONDITION_NUMBER": float(np.linalg.cond(ld)),
                "NOTE": "Smoke-test approximation; use susieR for scientific analysis",
            }
        ]
    )
    return FineMapResult(out, cs, summary, diagnostics)


def run_susier(
    variants: pl.DataFrame,
    ld: np.ndarray,
    *,
    sample_size: float,
    locus_id: str,
    population: str,
    workdir: Path,
    tools: ToolResolver,
    coverage: float = 0.95,
    max_effects: int = 10,
) -> FineMapResult:
    """Run the reference ``susieR::susie_rss`` implementation."""
    _validate_inputs(variants["Z"].to_numpy().astype(np.float64), ld)
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "variants.tsv"
    ld_path = workdir / "ld.tsv"
    output_path = workdir / "susie.json"
    variants.select("VARIANT_ID", "Z").write_csv(input_path, separator="\t")
    np.savetxt(ld_path, ld, delimiter="\t")
    script = workdir / "run_susie.R"
    script.write_text(
        """args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(jsonlite))
v <- read.delim(args[1], check.names=FALSE)
R <- as.matrix(read.delim(args[2], header=FALSE))
fit <- susie_rss(
  z=v$Z, R=R, n=as.numeric(args[4]), L=as.integer(args[5]),
  coverage=as.numeric(args[6]), estimate_residual_variance=FALSE
)
sets <- if (is.null(fit$sets$cs)) list() else lapply(fit$sets$cs, function(x) as.integer(x))
write_json(
  list(pip=as.numeric(fit$pip), sets=sets,
       converged=isTRUE(fit$converged), niter=fit$niter),
  args[3], auto_unbox=TRUE
)
""",
        encoding="utf-8",
    )
    run(
        [*tools.require("rscript"), str(script), str(input_path), str(ld_path), str(output_path), str(sample_size), str(max_effects), str(coverage)],
        log_file=workdir / "susie.log",
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    pip = np.asarray(payload["pip"], dtype=float)
    ids = variants["VARIANT_ID"].cast(pl.String).to_list()
    membership: list[list[str]] = [[] for _ in ids]
    cs_rows: list[dict[str, object]] = []
    for effect, members_one_based in enumerate(payload.get("sets", []), start=1):
        members = np.asarray(members_one_based, dtype=int) - 1
        cs_id = f"{locus_id}:CS{effect}"
        for idx in members:
            membership[int(idx)].append(cs_id)
        cs_rows.append(
            {
                "LOCUS_ID": locus_id,
                "POPULATION": population,
                "CS_ID": cs_id,
                "EFFECT": effect,
                "COVERAGE": None,
                "PURITY_MIN_ABS_COR": float(np.min(np.abs(ld[np.ix_(members, members)]))) if members.size > 1 else 1.0,
                "N_VARIANTS": int(members.size),
                "VARIANTS": [ids[i] for i in members],
            }
        )
    out = variants.with_columns(
        pl.Series("PIP", pip),
        pl.Series("CREDIBLE_SETS", membership, dtype=pl.List(pl.String)),
        pl.lit(locus_id).alias("LOCUS_ID"),
        pl.lit(population).alias("POPULATION"),
    )
    cs = pl.DataFrame(cs_rows) if cs_rows else pl.DataFrame()
    summary = pl.DataFrame([{"LOCUS_ID": locus_id, "POPULATION": population, "N_VARIANTS": variants.height, "N_EFFECTS": len(cs_rows), "N_CREDIBLE_SETS": len(cs_rows), "MAX_PIP": float(np.max(pip)), "STATUS": "OK" if payload.get("converged") else "NOT_CONVERGED"}])
    diagnostics = pl.DataFrame([{"LOCUS_ID": locus_id, "POPULATION": population, "ENGINE": "susieR", "PRODUCTION_SUSIE": True, "CONVERGED": bool(payload.get("converged")), "N_ITER": payload.get("niter"), "LD_MIN_EIGENVALUE": float(np.linalg.eigvalsh(ld).min()), "LD_CONDITION_NUMBER": float(np.linalg.cond(ld))}])
    return FineMapResult(out, cs, summary, diagnostics)
