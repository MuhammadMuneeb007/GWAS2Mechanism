"""Adult GTEx compact SuSiE eQTL/sQTL resource integration across all tissues."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

import httpx
import polars as pl


def requested_tissues(available: Iterable[str], configured: str | list[str] = "all") -> list[str]:
    found = sorted(set(available))
    if configured == "all":
        return found
    requested = set(configured)
    missing = requested - set(found)
    if missing:
        raise ValueError(f"Requested GTEx tissues are unavailable: {sorted(missing)}")
    return [tissue for tissue in found if tissue in requested]


def discover_release(listing_url: str, timeout: float = 30.0) -> str:
    """Discover the highest stable Adult GTEx release from its bucket listing."""
    response = httpx.get(listing_url, params={"prefix": "v"}, timeout=timeout)
    response.raise_for_status()
    names = [item.get("name", "") for item in response.json().get("items", [])]
    versions = sorted({int(match.group(1)) for name in names if (match := re.search(r"(?:^|/)v(\d+)(?:/|\.)", name))})
    if not versions:
        raise RuntimeError("No stable Adult GTEx release was found")
    return f"v{versions[-1]}"


def parse_gtex_variant(value: str) -> str | None:
    parts = str(value).replace("chr", "", 1).split("_")
    if len(parts) < 4 or not parts[1].isdigit():
        return None
    chrom = {"23": "X", "24": "Y", "25": "MT", "M": "MT"}.get(parts[0].upper(), parts[0].upper())
    return f"{chrom}:{int(parts[1])}:{parts[2].upper()}:{parts[3].upper()}"


def _normalise_qtl(frame: pl.DataFrame, tissue: str, qtl_type: str) -> pl.DataFrame:
    aliases = {
        "variant_id": "GTEX_VARIANT_ID", "variant": "GTEX_VARIANT_ID", "gene_id": "GENE_ID",
        "gene_name": "GENE", "phenotype_id": "QTL_PHENOTYPE", "pip": "QTL_PIP",
        "susie_pip": "QTL_PIP", "molecular_trait_id": "QTL_PHENOTYPE",
    }
    frame = frame.rename({key: value for key, value in aliases.items() if key in frame.columns and value not in frame.columns})
    if "GTEX_VARIANT_ID" not in frame.columns:
        raise ValueError("GTEx table lacks a variant_id column")
    if "QTL_PIP" not in frame.columns:
        frame = frame.with_columns(pl.lit(None, dtype=pl.Float64).alias("QTL_PIP"))
    if "GENE_ID" not in frame.columns:
        frame = frame.with_columns(pl.lit(None, dtype=pl.String).alias("GENE_ID"))
    if "GENE" not in frame.columns:
        frame = frame.with_columns(pl.col("GENE_ID").alias("GENE"))
    if "QTL_PHENOTYPE" not in frame.columns:
        frame = frame.with_columns(pl.col("GENE_ID").alias("QTL_PHENOTYPE"))
    return frame.with_columns(
        pl.col("GTEX_VARIANT_ID").cast(pl.String).map_elements(parse_gtex_variant, return_dtype=pl.String).alias("VARIANT_ID"),
        pl.lit(tissue).alias("TISSUE"),
        pl.lit(qtl_type).alias("QTL_TYPE"),
        pl.col("QTL_PIP").cast(pl.Float64, strict=False),
    ).filter(pl.col("VARIANT_ID").is_not_null()).select("VARIANT_ID", "TISSUE", "QTL_TYPE", "GENE_ID", "GENE", "QTL_PHENOTYPE", "QTL_PIP")


def integrate_local_archive(
    prioritised: pl.DataFrame,
    archive_dir: Path,
    *,
    tissues: str | list[str] = "all",
    qtl_types: tuple[str, ...] = ("eQTL", "sQTL"),
) -> pl.DataFrame:
    """Scan compact local Parquet resources with predicate pushdown."""
    files: list[tuple[Path, str, str]] = []
    for path in archive_dir.rglob("*.parquet"):
        lower = path.name.lower()
        qtl_type = "eQTL" if "eqtl" in lower else "sQTL" if "sqtl" in lower else "apaQTL" if "apa" in lower else ""
        if qtl_type not in qtl_types:
            continue
        tissue = re.sub(r"(?i)[._-]?(?:e|s|apa)qtl.*$", "", path.stem).strip("._-")
        files.append((path, tissue, qtl_type))
    selected_tissues = set(requested_tissues((item[1] for item in files), tissues))
    wanted = prioritised.select("VARIANT_ID").unique()
    matches = []
    for path, tissue, qtl_type in files:
        if tissue not in selected_tissues:
            continue
        normal = _normalise_qtl(pl.read_parquet(path), tissue, qtl_type)
        matches.append(normal.join(wanted, on="VARIANT_ID", how="semi"))
    return pl.concat(matches, how="diagonal_relaxed") if matches else pl.DataFrame(
        schema={"VARIANT_ID": pl.String, "TISSUE": pl.String, "QTL_TYPE": pl.String, "GENE_ID": pl.String, "GENE": pl.String, "QTL_PHENOTYPE": pl.String, "QTL_PIP": pl.Float64}
    )
