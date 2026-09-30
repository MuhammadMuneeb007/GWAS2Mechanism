"""Long-format evidence and cross-population comparison tables."""

from __future__ import annotations

import polars as pl


def integrate_evidence(
    finemap: pl.DataFrame,
    *,
    phenotype: str,
    annotation: pl.DataFrame | None = None,
    qtl: pl.DataFrame | None = None,
    coloc: pl.DataFrame | None = None,
    splicing: pl.DataFrame | None = None,
    cross_ancestry: pl.DataFrame | None = None,
) -> pl.DataFrame:
    frame = finemap.with_columns(pl.lit(phenotype).alias("PHENOTYPE"))
    for other in (cross_ancestry, annotation, qtl, coloc, splicing):
        if other is not None and not other.is_empty() and "VARIANT_ID" in other.columns:
            join_keys = ["VARIANT_ID"]
            for candidate in ("LOCUS_ID", "POPULATION", "GENE", "TISSUE", "QTL_TYPE"):
                if candidate in frame.columns and candidate in other.columns:
                    join_keys.append(candidate)
            overlap = (set(frame.columns) & set(other.columns)) - set(join_keys)
            if overlap:
                other = other.rename({column: f"{column}_EVIDENCE" for column in overlap})
            frame = frame.join(other, on=join_keys, how="left")
    return frame


def add_priority_components(frame: pl.DataFrame, weights: dict[str, float]) -> pl.DataFrame:
    """Optional transparent score: every normalized component is retained."""
    aliases = {
        "finemap_pip": "PIP",
        "cross_ancestry_pip": "CROSS_ANCESTRY_PIP",
        "lclpp": "LCLPP",
        "spliceai": "SPLICEAI_MAX",
        "pangolin": "PANGOLIN_SCORE",
        "n_supporting_populations": "N_SUPPORTING_POPULATIONS",
    }
    components = []
    out = frame
    for name, weight in weights.items():
        source = aliases.get(name, name.upper())
        target = f"PRIORITY_COMPONENT_{name.upper()}"
        if name == "coding_or_splice_consequence":
            expression = pl.col("FUNCTIONAL_CATEGORY").is_in(["CODING", "SPLICING"]).cast(pl.Float64) if "FUNCTIONAL_CATEGORY" in out.columns else pl.lit(0.0)
        elif source in out.columns:
            expression = pl.col(source).cast(pl.Float64, strict=False).fill_null(0.0).clip(0.0, 1.0)
        else:
            expression = pl.lit(0.0)
        out = out.with_columns((expression * float(weight)).alias(target))
        components.append(pl.col(target))
    return out.with_columns(pl.sum_horizontal(components).alias("PRIORITY_SCORE"))


def cross_population_comparison(finemap: pl.DataFrame, support_pip: float = 0.1) -> pl.DataFrame:
    if finemap.is_empty():
        return pl.DataFrame()
    populations = sorted(finemap["POPULATION"].drop_nulls().unique().to_list())
    pivot = finemap.select("VARIANT_ID", "POPULATION", "PIP").pivot(on="POPULATION", index="VARIANT_ID", values="PIP", aggregate_function="max")
    pivot = pivot.rename({population: f"{population}_PIP" for population in populations})
    pip_columns = [f"{population}_PIP" for population in populations]
    return pivot.with_columns(
        pl.sum_horizontal([pl.col(column).fill_null(0.0).ge(support_pip).cast(pl.Int64) for column in pip_columns]).alias("N_SUPPORTING_POPULATIONS"),
        pl.max_horizontal([pl.col(column) for column in pip_columns]).alias("MAX_POPULATION_PIP"),
        pl.mean_horizontal([pl.col(column) for column in pip_columns]).alias("MEAN_POPULATION_PIP"),
    )


def effect_direction_concordance(sumstats: pl.DataFrame) -> pl.DataFrame:
    return sumstats.group_by("VARIANT_ID").agg(
        pl.col("BETA").drop_nulls().sign().n_unique().eq(1).alias("EFFECT_DIRECTION_CONCORDANT"),
        pl.len().alias("N_POPULATIONS_WITH_EFFECT"),
        pl.col("N").alias("SAMPLE_SIZES"),
    )
