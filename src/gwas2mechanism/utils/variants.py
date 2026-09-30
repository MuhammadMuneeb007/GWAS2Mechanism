"""Vectorised variant normalisation expressions (Polars).

Variant identity throughout the framework is ``CHR:POS:REF:ALT`` on GRCh38 with
``CHR`` lacking a ``chr`` prefix and alleles upper-cased.
"""

from __future__ import annotations

import polars as pl

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
PALINDROMIC_PAIRS = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


def norm_chr(expr: pl.Expr) -> pl.Expr:
    """Strip ``chr``; map 23/24/25/M to X/Y/MT."""
    base = expr.cast(pl.String).str.strip_chars().str.replace(r"(?i)^chr", "")
    return (
        pl.when(base == "23")
        .then(pl.lit("X"))
        .when(base == "24")
        .then(pl.lit("Y"))
        .when(base.is_in(["25", "M", "m", "MT"]))
        .then(pl.lit("MT"))
        .otherwise(base.str.to_uppercase())
    )


def norm_allele(expr: pl.Expr) -> pl.Expr:
    return expr.cast(pl.String).str.strip_chars().str.to_uppercase()


def variant_id(chr_: str = "CHR", pos: str = "POS", ref: str = "REF", alt: str = "ALT") -> pl.Expr:
    return pl.concat_str(
        [pl.col(chr_), pl.col(pos).cast(pl.String), pl.col(ref), pl.col(alt)], separator=":"
    )


def gtex_variant_id(
    chr_: str = "CHR", pos: str = "POS", ref: str = "REF", alt: str = "ALT"
) -> pl.Expr:
    """GTEx GRCh38 variant id: ``chr1_12345_A_G_b38``."""
    return pl.concat_str(
        [
            pl.lit("chr") + pl.col(chr_),
            pl.col(pos).cast(pl.String),
            pl.col(ref),
            pl.col(alt),
            pl.lit("b38"),
        ],
        separator="_",
    )


def parse_gtex_variant(expr: pl.Expr) -> pl.Expr:
    """``chr1_12345_A_G_b38`` -> ``1:12345:A:G`` (null when unparsable)."""
    parts = expr.str.extract_groups(r"^(?:chr)?([^_]+)_(\d+)_([ACGTN]+)_([ACGTN]+)_(?:b38|GRCh38)$")
    return pl.concat_str(
        [
            norm_chr(parts.struct.field("1")),
            parts.struct.field("2"),
            parts.struct.field("3"),
            parts.struct.field("4"),
        ],
        separator=":",
    )


def allele_pair_key(a: str, b: str) -> pl.Expr:
    """Order-independent allele pair ``A1/A2`` with A1 <= A2 lexicographically."""
    lo = pl.min_horizontal(pl.col(a), pl.col(b))
    hi = pl.max_horizontal(pl.col(a), pl.col(b))
    return pl.concat_str([lo, hi], separator="/")


def is_palindromic(a: str, b: str) -> pl.Expr:
    pair = pl.concat_str([pl.col(a), pl.col(b)])
    return pair.is_in(["AT", "TA", "CG", "GC"])


def complement(expr: pl.Expr) -> pl.Expr:
    return expr.str.replace_all("A", "t").str.replace_all("T", "a").str.replace_all(
        "C", "g"
    ).str.replace_all("G", "c").str.to_uppercase()


def complement_py(allele: str) -> str:
    return allele.upper().translate(_COMPLEMENT)


def is_valid_dna(expr: pl.Expr) -> pl.Expr:
    return expr.str.contains(r"^[ACGT]+$")


def split_variant_id(frame: pl.DataFrame | pl.LazyFrame, column: str = "VARIANT_ID"):
    parts = pl.col(column).str.split_exact(":", 3)
    return frame.with_columns(
        parts.struct.field("field_0").alias("CHR"),
        parts.struct.field("field_1").cast(pl.Int64).alias("POS"),
        parts.struct.field("field_2").alias("REF"),
        parts.struct.field("field_3").alias("ALT"),
    )
