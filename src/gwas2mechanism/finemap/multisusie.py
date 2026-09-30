"""MultiSuSiE adapter retaining one summary-statistic vector and LD matrix per ancestry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl


@dataclass
class MultiSuSiEInputs:
    populations: list[str]
    z: list[np.ndarray]
    ld: list[np.ndarray]
    sample_sizes: list[float]
    variant_ids: list[str]

    def validate(self) -> None:
        if len(self.populations) < 2:
            raise ValueError("Multi-ancestry fine-mapping requires at least two populations")
        if not (len(self.z) == len(self.ld) == len(self.sample_sizes) == len(self.populations)):
            raise ValueError("One z vector, LD matrix and sample size are required per population")
        m = len(self.variant_ids)
        for population, z_one, ld_one in zip(self.populations, self.z, self.ld, strict=True):
            if z_one.shape != (m,) or ld_one.shape != (m, m):
                raise ValueError(f"Dimension mismatch for {population}")
        if len({id(matrix) for matrix in self.ld}) != len(self.ld):
            raise ValueError("Each population must supply its own LD matrix object")


def run_multisusie(inputs: MultiSuSiEInputs, *, coverage: float = 0.95, rho: float = 0.75) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Call the installed MultiSuSiE package and normalise its public outputs."""
    inputs.validate()
    try:
        from multisusie.multisusie import multisusie_rss  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("MultiSuSiE is not installed; install the pinned optional dependency") from exc
    fit = multisusie_rss(
        z_list=inputs.z,
        R_list=inputs.ld,
        n_list=inputs.sample_sizes,
        rho=rho,
        coverage=coverage,
    )
    pip = np.asarray(getattr(fit, "pip", fit["pip"]), dtype=float)
    variants = pl.DataFrame({"VARIANT_ID": inputs.variant_ids, "CROSS_ANCESTRY_PIP": pip})
    summary = pl.DataFrame([{"METHOD": "MultiSuSiE", "POPULATIONS": inputs.populations, "N_POPULATIONS": len(inputs.populations), "N_VARIANTS": len(inputs.variant_ids), "MAX_PIP": float(pip.max()), "STATUS": "OK"}])
    return variants, summary
