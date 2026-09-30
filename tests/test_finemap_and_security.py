from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from gwas2mechanism.finemap.multisusie import MultiSuSiEInputs
from gwas2mechanism.finemap.susie import run_numpy_smoke_susie


def test_population_specific_finemap_smoke() -> None:
    variants = pl.DataFrame({"VARIANT_ID": ["a", "b", "c"], "Z": [1.0, 8.0, 0.5]})
    result = run_numpy_smoke_susie(variants, np.eye(3), locus_id="L1", population="AFR", max_effects=1)
    assert result.variants.sort("PIP", descending=True)["VARIANT_ID"][0] == "b"
    assert result.diagnostics["PRODUCTION_SUSIE"][0] is False


def test_multisusie_requires_distinct_population_ld() -> None:
    ld = np.eye(2)
    inputs = MultiSuSiEInputs(["AFR", "EUR"], [np.ones(2), np.ones(2)], [ld, ld], [1000, 2000], ["a", "b"])
    with pytest.raises(ValueError, match="own LD matrix"):
        inputs.validate()


def test_no_inappropriate_core_defaults() -> None:
    root = Path(__file__).parents[1] / "src" / "gwas2mechanism"
    forbidden = [r"coronary artery disease", r"Heart_Atrial_Appendage", r"68 expected loci", r"55,251 cases", r"402,261 controls"]
    # Avoid matching ordinary words that merely contain the three letters.
    forbidden.append(r"\bCAD\b")
    offenders = []
    for path in root.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        for pattern in forbidden:
            if re.search(pattern, text, flags=re.IGNORECASE):
                offenders.append((path.relative_to(root), pattern))
    assert not offenders


def test_no_hardcoded_secrets() -> None:
    root = Path(__file__).parents[1]
    pattern = re.compile(r"(?:sk-|ghp_|gho_)[A-Za-z0-9_]{12,}")
    offenders = [path for path in root.rglob("*") if path.is_file() and ".git" not in path.parts and pattern.search(path.read_text(encoding="utf-8", errors="ignore"))]
    assert not offenders
