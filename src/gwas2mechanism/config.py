"""Typed configuration.

The default configuration lives in ``config/default.yaml`` (bundled into the
wheel). User configuration files and ``--set key=value`` overrides are deep
merged on top and validated with Pydantic, with unknown keys rejected so a typo
never silently falls back to a default.
"""

from __future__ import annotations

import copy
import os
from importlib import resources
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class PhenotypeResolutionConfig(_Strict):
    synonyms: list[str] = Field(default_factory=list)
    efo_ids: list[str] = Field(default_factory=list)
    include_descendants: bool = False
    allow_fuzzy: bool = False
    fuzzy_min_ratio: float = 0.92
    ols_url: str = "https://www.ebi.ac.uk/ols4/api"


class GwasCatalogConfig(_Strict):
    mode: Literal["rest", "bulk"] = "rest"
    rest_url: str
    studies_source: str
    ancestries_source: str
    sumstats_base: str
    page_size: int = 200
    max_studies: int = 2000


class DiscoveryConfig(_Strict):
    backend: Literal["auto", "gwaspoker", "gwas_catalog"] = "auto"
    gwas_catalog: GwasCatalogConfig
    list_files: bool = True
    probe_headers: bool = True
    probe_top_n: int = 5
    probe_bytes: int = 262144
    timeout_s: float = 60
    concurrency: int = 8
    all_studies: bool = False
    prefer_harmonised: bool = True
    exclude_gxe: bool = True
    exclude_background_traits: bool = True


class RankingConfig(_Strict):
    weights: dict[str, float]
    penalties: dict[str, float]
    min_phenotype_match: float = 0.8


class PopulationMapEntry(_Strict):
    population: str
    approximate: bool = False


class PopulationsConfig(_Strict):
    mode: str = "auto"
    allowed: list[str] = Field(default_factory=lambda: ["EUR", "AFR", "EAS", "SAS", "AMR"])
    mapping: dict[str, PopulationMapEntry]


class CombinedConfig(_Strict):
    published_combined: bool = True
    fixed_effect_meta: bool = True
    multiancestry_finemap: bool = True
    min_populations_for_meta: int = 2
    finemap_published_with_population: str | None = None


class GenomeConfig(_Strict):
    build: Literal["GRCh38"] = "GRCh38"
    chromosomes: list[str]


class HarmonizationConfig(_Strict):
    engine: Literal["auto", "gwaslab", "builtin"] = "auto"
    require_grch38: bool = True
    liftover: Literal["gwaslab", "none"] = "gwaslab"
    min_info: float | None = None
    drop_palindromic_ambiguous_maf: float = 0.4
    validate_build_with_fasta: bool = True
    build_validation_sample: int = 20000
    build_validation_min_concordance: float = 0.9


class GwasConfig(_Strict):
    p_threshold: float = 5e-8


class ClumpingConfig(_Strict):
    engine: Literal["plink2", "builtin"] = "plink2"
    p1: float = 5e-8
    p2: float = 1.0
    r2: float = 0.1
    kb: int = 1000


class LociConfig(_Strict):
    window_bp: int = 1_000_000
    combined_merge_bp: int = 0


class ReferenceConfig(_Strict):
    panel: str
    base_url: str
    panel_url: str
    vcf_template: str
    ld_engine: Literal["plink2", "dosage"] = "plink2"
    maf: float = 0.01
    geno: float = 0.05
    keep_population_vcf: bool = False


class FinemappingConfig(_Strict):
    method: Literal["susie_rss"] = "susie_rss"
    engine: Literal["susier", "numpy", "auto"] = "susier"
    credible_set_coverage: float = 0.95
    min_purity: float = 0.5
    min_variants: int = 10
    max_ld_variants: int = 15000
    ld_ridge: float = 1e-4
    default_L: int = 10
    max_L: int = 20
    estimate_s_max_variants: int = 5000
    refine: bool = False
    max_iter: int = 100
    tol: float = 1e-3
    cross_ancestry_method: Literal["multisusie", "susiex", "none"] = "multisusie"
    multisusie_rho: float = 0.75
    n_jobs: int | Literal["auto"] = "auto"


class MetaConfig(_Strict):
    engine: Literal["builtin", "metal"] = "builtin"
    mrmega: bool = False


class AnnotationConfig(_Strict):
    vep: bool = True
    min_pip: float = 0.01
    species: str = "homo_sapiens"
    vep_cache_version: int | None = None
    vep_extra_args: list[str] = Field(default_factory=list)


class QtlConfig(_Strict):
    source: Literal["gtex"] = "gtex"
    release: str = "latest"
    tissues: Literal["all"] | list[str] = "all"
    types: list[Literal["eqtl", "sqtl", "apaqtl"]] = Field(default_factory=lambda: ["eqtl", "sqtl"])
    bucket_listing_url: str
    bucket_download_url: str
    local_archive_dir: str | None = None


class SplicingThresholds(_Strict):
    high_recall: float = 0.2
    recommended: float = 0.5
    high_precision: float = 0.8


class SplicingConfig(_Strict):
    spliceai: bool = True
    pangolin: bool = True
    spliceai_distance: int = 50
    spliceai_mask: int = 0
    spliceai_annotation: str = "grch38"
    pangolin_distance: int = 500
    pangolin_mask: bool = False
    thresholds: SplicingThresholds = Field(default_factory=SplicingThresholds)
    pangolin_ensembl_release: int = 116


class FormalColocConfig(_Strict):
    enabled: bool = False
    top_pairs: int = 100
    max_download_mb: int = 500


class ColocalizationConfig(_Strict):
    pip_screen: bool = True
    screen_classes: dict[str, float]
    formal: FormalColocConfig = Field(default_factory=FormalColocConfig)


class PriorityScoreConfig(_Strict):
    enabled: bool = True
    weights: dict[str, float]


class EvidenceConfig(_Strict):
    priority_score: PriorityScoreConfig
    support_pip: float = 0.1


class OptionalModulesConfig(_Strict):
    lncrna: bool = False


class LncrnaConfig(_Strict):
    gencode_release: int = 50
    pip_thresholds: list[float] = Field(default_factory=lambda: [0.01, 0.1, 0.5])


class PerformanceConfig(_Strict):
    parquet_compression: str = "zstd"
    parquet_compression_level: int = 3
    threads: int | Literal["auto"] = "auto"
    memory_safe: bool = True
    blas_threads_per_job: int = 1


class ReportConfig(_Strict):
    html: bool = True
    markdown: bool = True
    tsv: bool = True
    top_n: int = 50


class Config(_Strict):
    phenotype: str | None = None
    mode: Literal["fast", "full"] = "fast"
    run_root: str = "runs"
    resources_dir: str | None = None
    phenotype_resolution: PhenotypeResolutionConfig
    discovery: DiscoveryConfig
    ranking: RankingConfig
    populations: PopulationsConfig
    combined: CombinedConfig
    genome: GenomeConfig
    harmonization: HarmonizationConfig
    gwas: GwasConfig
    clumping: ClumpingConfig
    loci: LociConfig
    reference: ReferenceConfig
    finemapping: FinemappingConfig
    meta: MetaConfig
    annotation: AnnotationConfig
    qtl: QtlConfig
    splicing: SplicingConfig
    colocalization: ColocalizationConfig
    evidence: EvidenceConfig
    optional_modules: OptionalModulesConfig
    lncrna: LncrnaConfig
    performance: PerformanceConfig
    report: ReportConfig
    tools: dict[str, str | list[str]]

    @field_validator("tools")
    @classmethod
    def _tools_non_empty(cls, value: dict[str, Any]) -> dict[str, Any]:
        for name, cmd in value.items():
            if cmd in ("", []):
                raise ValueError(f"tools.{name} must not be empty")
        return value

    # ------------------------------------------------------------------
    @property
    def resources_path(self) -> Path:
        if self.resources_dir:
            return Path(self.resources_dir).expanduser()
        env = os.environ.get("GWAS2M_CACHE")
        if env:
            return Path(env).expanduser()
        return Path.home() / ".cache" / "gwas2mechanism"

    def dump(self) -> dict[str, Any]:
        return self.model_dump(mode="json")


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------


def package_root() -> Path:
    """Repository root when running from a source checkout."""
    return Path(__file__).resolve().parents[2]


def bundled_path(*parts: str) -> Path:
    """Locate a file shipped via hatch force-include, else in the source tree."""
    try:
        candidate = Path(str(resources.files("gwas2mechanism").joinpath("_bundled", *parts)))
        if candidate.exists():
            return candidate
    except (ModuleNotFoundError, FileNotFoundError):  # pragma: no cover - defensive
        pass
    return package_root().joinpath(*parts)


def default_config_path() -> Path:
    return bundled_path("config", "default.yaml")


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _coerce_scalar(text: str) -> Any:
    parsed = yaml.safe_load(text)
    return parsed


def parse_set_overrides(items: list[str] | None) -> dict[str, Any]:
    """Turn ``["a.b=1", "c=[x, y]"]`` into a nested dict (YAML-typed values)."""
    result: dict[str, Any] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--set expects key=value, got {item!r}")
        key, value = item.split("=", 1)
        node = result
        parts = key.strip().split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = _coerce_scalar(value)
    return result


def load_config(
    paths: list[Path] | None = None,
    overrides: dict[str, Any] | None = None,
) -> Config:
    with open(default_config_path(), encoding="utf-8") as handle:
        data: dict[str, Any] = yaml.safe_load(handle)
    for path in paths or []:
        with open(path, encoding="utf-8") as handle:
            data = deep_merge(data, yaml.safe_load(handle) or {})
    if overrides:
        data = deep_merge(data, overrides)
    return Config.model_validate(data)


def save_config(config: Config, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(config.dump(), handle, sort_keys=False)
