"""Harmonisation to the canonical GRCh38 schema (Polars lazy).

Engines
-------
* ``gwaslab``  - GWASLab format detection, ``basic_check`` QC/normalisation and
  (if needed) liftover. GWASLab is pandas-based: we convert to Polars
  immediately after it returns.
* ``builtin``  - Polars lazy scan with column-alias recognition.
* ``auto``     - GWASLab when importable, otherwise builtin.

Statistical conventions (never violated)
----------------------------------------
* ``BETA = ln(OR)`` when only OR is available; a raw OR is never used as BETA.
* ``Z = BETA / SE``; P recomputed from Z only when missing; ``NEG_LOG10_P`` is
  computed from Z via ``scipy.stats.norm.logsf`` so P-values that underflow to 0
  keep their magnitude.
* ``N_EFF = 4 / (1/N_cases + 1/N_controls)`` for binary traits (per-SNP counts
  win over study-level metadata); quantitative traits keep ``N``.
* Sample sizes are never fabricated: unknown stays null.
* An aggregated multi-ancestry file is **never** relabelled as a single ancestry.
"""

from __future__ import annotations

import gzip
import logging
import math
import shutil
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy import stats

from gwas2mechanism.config import Config
from gwas2mechanism.constants import COMBINED_PUBLISHED, MULTI, SUMSTATS_COLUMNS, SUPERPOPULATIONS
from gwas2mechanism.discovery.base import parse_header
from gwas2mechanism.sumstats.schema import CANONICAL_DTYPES, map_columns
from gwas2mechanism.utils.fasta import FastaReader
from gwas2mechanism.utils.io import atomic_write_json, write_parquet
from gwas2mechanism.utils.variants import allele_pair_key, is_valid_dna, norm_allele, norm_chr

log = logging.getLogger(__name__)

LN10 = math.log(10.0)
MEDIAN_CHI2_1DF = 0.454936423119572


class AggregatedSplitError(RuntimeError):
    """Raised if anything attempts to treat aggregated multi-ancestry statistics
    as ancestry-specific. Summary statistics cannot be split retrospectively."""


class BuildMismatchError(RuntimeError):
    pass


def guard_population_label(analysis: dict[str, Any]) -> str:
    """Return the population label for an analysis, refusing any split of MULTI data."""
    population = analysis.get("population")
    label = analysis.get("analysis_label")
    if population == MULTI:
        if label != COMBINED_PUBLISHED:
            raise AggregatedSplitError(
                f"{analysis.get('study_accession')}: aggregated multi-ancestry summary statistics "
                f"cannot be analysed as {label!r}; they may only be used as {COMBINED_PUBLISHED}."
            )
        return COMBINED_PUBLISHED
    if label in SUPERPOPULATIONS and population != label:
        raise AggregatedSplitError(f"analysis label {label} does not match study population {population}")
    return label


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


def _decompress(path: Path, workdir: Path) -> Path:
    """Decompress once to a scratch file so Polars can scan lazily."""
    name = path.name.lower()
    if name.endswith((".gz", ".bgz")):
        out = workdir / (path.name.rsplit(".", 1)[0] + ".uncompressed")
        if not out.exists():
            tmp = out.with_suffix(".part")
            pigz = shutil.which("pigz")
            if pigz:
                import subprocess

                with open(tmp, "wb") as handle:
                    subprocess.run([pigz, "-dc", str(path)], stdout=handle, check=True)
            else:
                with gzip.open(path, "rb") as src, open(tmp, "wb") as dst:
                    shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
            tmp.replace(out)
        return out
    if name.endswith(".zip"):
        out = workdir / (path.stem + ".uncompressed")
        if not out.exists():
            with zipfile.ZipFile(path) as archive:
                member = next(m for m in archive.namelist() if not m.endswith("/"))
                with archive.open(member) as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        return out
    return path


def _detect_separator(path: Path) -> tuple[str, list[str], int]:
    with open(path, encoding="utf-8", errors="replace") as handle:
        skip = 0
        for line in handle:
            if line.startswith("##"):
                skip += 1
                continue
            header = parse_header(line)
            if "\t" in line:
                return "\t", header, skip
            if "," in line and len(line.split(",")) >= 4:
                return ",", header, skip
            return " ", header, skip
    raise ValueError(f"{path}: no header line found")


def scan_builtin(path: Path, workdir: Path) -> tuple[pl.LazyFrame, dict[str, str]]:
    if path.name.lower().endswith((".vcf", ".vcf.gz")):
        raise ValueError("GWAS-VCF input requires the GWASLab engine")
    plain = _decompress(path, workdir)
    separator, header, skip = _detect_separator(plain)
    mapping = map_columns(header)
    missing = [c for c in ("CHR", "POS", "EA", "NEA") if c not in mapping]
    if missing:
        raise ValueError(f"{path.name}: could not identify required columns {missing}; header={header[:30]}")
    if separator == " ":
        lf = pl.scan_csv(plain, separator=" ", skip_rows=skip, infer_schema=False, truncate_ragged_lines=True)
    else:
        lf = pl.scan_csv(
            plain, separator=separator, skip_rows=skip, infer_schema=False, quote_char=None, truncate_ragged_lines=True
        )
    lf = lf.select([pl.col(src).alias(dst) for dst, src in mapping.items()])
    return lf, mapping


GWASLAB_RENAME = {
    "rsID": "RSID",
    "CHR": "CHR",
    "POS": "POS",
    "EA": "EA",
    "NEA": "NEA",
    "BETA": "BETA",
    "SE": "SE",
    "P": "P",
    "MLOG10P": "NEG_LOG10_P",
    "EAF": "EAF",
    "N": "N",
    "N_CASE": "N_CASES",
    "N_CONTROL": "N_CONTROLS",
    "OR": "OR",
    "OR_95L": "OR_L95",
    "OR_95U": "OR_U95",
    "Z": "Z",
    "INFO": "INFO",
}


def scan_gwaslab(path: Path, build: str | None, threads: int, liftover: bool) -> tuple[pl.LazyFrame, dict[str, str]]:
    import gwaslab as gl  # pandas boundary

    source_build = "38" if (build or "GRCh38") == "GRCh38" else "19"
    ss = gl.Sumstats(str(path), fmt="auto", build=source_build, verbose=False)
    ss.basic_check(remove=True, remove_dup=True, normalize=True, threads=threads, verbose=False)
    if source_build != "38":
        if not liftover:
            raise BuildMismatchError(f"{path.name} is on build {build}; liftover disabled")
        ss.liftover(n_cores=threads, from_build="19", to_build="38", remove=True)
    frame = pl.from_pandas(ss.data.reset_index(drop=True))  # back to Polars immediately
    del ss
    keep = {src: dst for src, dst in GWASLAB_RENAME.items() if src in frame.columns}
    frame = frame.select([pl.col(src).cast(pl.String).alias(dst) for src, dst in keep.items()])
    return frame.lazy(), {dst: src for src, dst in keep.items()}


# ----------------------------------------------------------------------------
# Statistics derivation (pure, unit-tested)
# ----------------------------------------------------------------------------


def _p_from_z(z: pl.Series) -> pl.Series:
    arr = np.abs(z.to_numpy().astype(np.float64))
    return pl.Series(2.0 * stats.norm.sf(arr))


def _neglog10p_from_z(z: pl.Series) -> pl.Series:
    arr = np.abs(z.to_numpy().astype(np.float64))
    return pl.Series(-(math.log(2.0) + stats.norm.logsf(arr)) / LN10)


def derive_statistics(lf: pl.LazyFrame, meta: dict[str, Any]) -> pl.LazyFrame:
    """Cast, normalise and derive BETA/SE/Z/P/N fields. Vectorised end-to-end."""
    schema = lf.collect_schema()
    num = lambda c: pl.col(c).cast(pl.Float64, strict=False) if c in schema else pl.lit(None, pl.Float64)  # noqa: E731
    string = lambda c: pl.col(c).cast(pl.String) if c in schema else pl.lit(None, pl.String)  # noqa: E731

    odds = num("OR")
    beta = pl.when(num("BETA").is_not_null()).then(num("BETA")).when(odds > 0).then(odds.log()).otherwise(None)
    se_from_ci = (num("OR_U95").log() - num("OR_L95").log()) / (2.0 * stats.norm.isf(0.025))
    se = pl.coalesce(num("SE"), pl.when((num("OR_L95") > 0) & (num("OR_U95") > 0)).then(se_from_ci))

    n_cases_meta = meta.get("n_cases")
    n_controls_meta = meta.get("n_controls")
    n_total_meta = meta.get("n_total")
    n_cases = pl.coalesce(num("N_CASES"), pl.lit(n_cases_meta, pl.Float64))
    n_controls = pl.coalesce(num("N_CONTROLS"), pl.lit(n_controls_meta, pl.Float64))
    n = pl.coalesce(num("N"), n_cases + n_controls, pl.lit(n_total_meta, pl.Float64))
    n_eff = pl.when((n_cases > 0) & (n_controls > 0)).then(4.0 / (1.0 / n_cases + 1.0 / n_controls))

    out = lf.with_columns(
        CHR=norm_chr(string("CHR")),
        POS=num("POS").cast(pl.Int64, strict=False),
        RSID=pl.when(string("RSID").str.contains(r"^rs\d+$")).then(string("RSID")),
        EA=norm_allele(string("EA")),
        NEA=norm_allele(string("NEA")),
        BETA=beta,
        SE=se,
        OR=odds,
        EAF=num("EAF"),
        N=n,
        N_CASES=n_cases,
        N_CONTROLS=n_controls,
        N_EFF=n_eff,
        P_IN=num("P"),
        NLP_IN=num("NEG_LOG10_P"),
        Z_IN=num("Z"),
        INFO=num("INFO"),
    ).with_columns(
        Z=pl.coalesce(pl.when(pl.col("SE") > 0).then(pl.col("BETA") / pl.col("SE")), pl.col("Z_IN")),
    )
    out = out.with_columns(
        P=pl.coalesce(
            pl.col("P_IN"),
            pl.when(pl.col("NLP_IN").is_not_null()).then(pl.lit(10.0).pow(-pl.col("NLP_IN"))),
            pl.col("Z").map_batches(_p_from_z, return_dtype=pl.Float64),
        ),
    ).with_columns(
        NEG_LOG10_P=pl.when(pl.col("NLP_IN").is_not_null())
        .then(pl.col("NLP_IN"))
        .when(pl.col("P") > 1e-300)
        .then(-pl.col("P").log10())
        .otherwise(pl.col("Z").map_batches(_neglog10p_from_z, return_dtype=pl.Float64)),
    )
    return out


def apply_qc(lf: pl.LazyFrame, config: Config) -> pl.LazyFrame:
    chroms = config.genome.chromosomes
    filters = (
        pl.col("CHR").is_in(chroms)
        & pl.col("POS").is_not_null()
        & (pl.col("POS") > 0)
        & is_valid_dna(pl.col("EA"))
        & is_valid_dna(pl.col("NEA"))
        & (pl.col("EA") != pl.col("NEA"))
        & pl.col("BETA").is_finite()
        & pl.col("SE").is_finite()
        & (pl.col("SE") > 0)
        & pl.col("P").is_not_null()
        & (pl.col("P") >= 0)
        & (pl.col("P") <= 1)
        & (pl.col("EAF").is_null() | pl.col("EAF").is_between(0.0, 1.0))
    )
    if config.harmonization.min_info is not None:
        filters = filters & (pl.col("INFO").is_null() | (pl.col("INFO") >= config.harmonization.min_info))
    return lf.filter(filters)


def finalize(lf: pl.LazyFrame, *, study_id: str, trait: str, population: str, source_file: str) -> pl.LazyFrame:
    out = lf.with_columns(
        STUDY_ID=pl.lit(study_id),
        TRAIT=pl.lit(trait),
        POPULATION=pl.lit(population),
        REF=pl.lit(None, pl.String),
        ALT=pl.lit(None, pl.String),
        BUILD=pl.lit("GRCh38"),
        SOURCE_FILE=pl.lit(source_file),
        VARIANT_ID=pl.concat_str(
            [pl.col("CHR"), pl.col("POS").cast(pl.String), allele_pair_key("EA", "NEA").str.replace("/", ":")],
            separator=":",
        ),
    )
    # One row per variant: keep the strongest association for duplicates.
    out = out.sort("P").unique(subset=["VARIANT_ID"], keep="first", maintain_order=True)
    return out.select([pl.col(c).cast(CANONICAL_DTYPES[c]) for c in SUMSTATS_COLUMNS]).sort(["CHR", "POS"])


# ----------------------------------------------------------------------------
# Build validation against the GRCh38 FASTA
# ----------------------------------------------------------------------------


def validate_build(frame: pl.DataFrame, fasta_path: Path, sample: int, min_concordance: float) -> dict[str, Any]:
    snvs = frame.filter((pl.col("EA").str.len_chars() == 1) & (pl.col("NEA").str.len_chars() == 1))
    if snvs.is_empty():
        return {"status": "skipped", "reason": "no SNVs"}
    snvs = snvs.sample(n=min(sample, snvs.height), seed=1)
    hits = total = 0
    with FastaReader(fasta_path) as fasta:
        for chrom, pos, ea, nea in snvs.select("CHR", "POS", "EA", "NEA").iter_rows():
            base = fasta.fetch(chrom, int(pos), 1)
            if base is None:
                continue
            total += 1
            hits += base in (ea, nea)
    concordance = hits / total if total else float("nan")
    result = {"status": "checked", "n_checked": total, "ref_concordance": concordance}
    if total and concordance < min_concordance:
        raise BuildMismatchError(
            f"Only {concordance:.1%} of sampled SNVs match the GRCh38 reference; the file is probably "
            "not on GRCh38 (or positions are 0-based)."
        )
    return result


# ----------------------------------------------------------------------------
# Stage entry point
# ----------------------------------------------------------------------------


def harmonize_analysis(
    analysis: dict[str, Any],
    download: dict[str, Any],
    out_sumstats: Path,
    out_significant: Path,
    out_qc: Path,
    config: Config,
    trait: str,
    threads: int,
    fasta: Path | None = None,
) -> dict[str, Any]:
    population = guard_population_label(analysis)
    source = Path(download["local_path"])
    workdir = out_sumstats.parent / "tmp"
    workdir.mkdir(parents=True, exist_ok=True)

    engine = config.harmonization.engine
    if engine == "auto":
        try:
            import gwaslab  # noqa: F401

            engine = "gwaslab"
        except ImportError:
            engine = "builtin"
    build = analysis.get("build") or ("GRCh38" if analysis.get("harmonised_available") else None)
    if engine == "gwaslab":
        lf, mapping = scan_gwaslab(source, build, threads, config.harmonization.liftover == "gwaslab")
    else:
        if config.harmonization.require_grch38 and build not in (None, "GRCh38"):
            raise BuildMismatchError(f"{source.name} is {build}; the builtin engine cannot liftover (use GWASLab)")
        lf, mapping = scan_builtin(source, workdir)

    meta = {k: analysis.get(k) for k in ("n_cases", "n_controls", "n_total")}
    derived = derive_statistics(lf, meta)
    n_input = derived.select(pl.len()).collect().item()

    beta_check = derived.select(
        pl.col("BETA").min().alias("min"), pl.col("BETA").median().alias("median")
    ).collect().row(0, named=True)
    if "BETA" in mapping and beta_check["min"] is not None and beta_check["min"] > 0 and 0.5 < beta_check["median"] < 2:
        raise ValueError(
            f"{source.name}: the column mapped to BETA is strictly positive with median "
            f"{beta_check['median']:.3f}; it looks like an odds ratio. Refusing to use OR as BETA."
        )

    qc = apply_qc(derived, config)
    final = finalize(
        qc,
        study_id=analysis["study_accession"],
        trait=trait,
        population=population,
        source_file=download.get("source_url") or str(source),
    )
    write_parquet(final, out_sumstats)
    written = pl.scan_parquet(out_sumstats)
    significant = written.filter(pl.col("P") < config.gwas.p_threshold).sort("P")
    write_parquet(significant, out_significant)

    summary = written.select(
        pl.len().alias("n_variants"),
        (pl.col("Z") ** 2).median().alias("median_chi2"),
        pl.col("EAF").is_not_null().mean().alias("frac_eaf"),
        pl.col("N").is_not_null().mean().alias("frac_n"),
        pl.col("RSID").is_not_null().mean().alias("frac_rsid"),
        (pl.col("P") < config.gwas.p_threshold).sum().alias("n_significant"),
    ).collect().row(0, named=True)
    build_check: dict[str, Any] = {"status": "skipped", "reason": "no GRCh38 FASTA available"}
    if fasta is not None and fasta.exists() and config.harmonization.validate_build_with_fasta:
        build_check = validate_build(
            pl.read_parquet(out_sumstats, columns=["CHR", "POS", "EA", "NEA"]),
            fasta,
            config.harmonization.build_validation_sample,
            config.harmonization.build_validation_min_concordance,
        )
    qc_report = {
        "analysis_id": analysis["analysis_id"],
        "study_accession": analysis["study_accession"],
        "population": population,
        "engine": engine,
        "source_file": str(source),
        "column_mapping": mapping,
        "derived": {
            "beta_from_or": "BETA" not in mapping and "OR" in mapping,
            "se_from_ci": "SE" not in mapping and "OR_L95" in mapping,
            "n_from_metadata": "N" not in mapping and meta.get("n_total") is not None,
        },
        "n_input": n_input,
        "n_after_qc": summary["n_variants"],
        "n_removed": n_input - summary["n_variants"],
        "n_significant": summary["n_significant"],
        "lambda_gc": (summary["median_chi2"] / MEDIAN_CHI2_1DF) if summary["median_chi2"] else None,
        "fraction_with_eaf": summary["frac_eaf"],
        "fraction_with_n": summary["frac_n"],
        "fraction_with_rsid": summary["frac_rsid"],
        "n_cases": meta.get("n_cases"),
        "n_controls": meta.get("n_controls"),
        "n_total": meta.get("n_total"),
        "build": "GRCh38",
        "build_validation": build_check,
    }
    atomic_write_json(out_qc, qc_report)
    shutil.rmtree(workdir, ignore_errors=True)
    return qc_report
