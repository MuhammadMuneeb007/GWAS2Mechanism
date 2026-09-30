from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from gwas2mechanism.coloc.screen import locus_clpp, variant_clpp
from gwas2mechanism.config import load_config
from gwas2mechanism.discovery.ancestry import classify_study, parse_ancestry_string
from gwas2mechanism.discovery.models import AncestryRecord, CandidateStudy, OntologyTerm
from gwas2mechanism.discovery.phenotype import PhenotypeResolver, match_study
from gwas2mechanism.discovery.ranking import effective_n
from gwas2mechanism.finemap.region import match_to_reference
from gwas2mechanism.meta.fixed_effect import inverse_variance_meta
from gwas2mechanism.qtl.gtex import parse_gtex_variant, requested_tissues
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.sumstats.harmonize import (
    AggregatedSplitError,
    derive_statistics,
    guard_population_label,
    validate_build,
)
from gwas2mechanism.utils.tools import ToolResolver
from gwas2mechanism.utils.variants import split_variant_id, variant_id


def test_phenotype_exact_mapped_trait_does_not_broaden() -> None:
    cfg = load_config().phenotype_resolution
    migraine = OntologyTerm("EFO_0003821", "migraine", ("migraine disorder",), source="fixture")
    headache = OntologyTerm("EFO_0003843", "headache", source="fixture")
    resolution = PhenotypeResolver(cfg, extra_terms=[migraine, headache]).resolve("migraine")
    exact = CandidateStudy("GCST1", reported_trait="Migraine", mapped_trait_ids=[migraine.term_id], mapped_trait_labels=[migraine.label])
    broad = CandidateStudy("GCST2", reported_trait="Headache", mapped_trait_ids=[headache.term_id], mapped_trait_labels=[headache.label])
    assert match_study(exact, resolution).match_type == "exact_mapped_trait"
    assert match_study(broad, resolution).match_score == 0.0


def test_population_classification_and_combined_published() -> None:
    cfg = load_config()
    eur = CandidateStudy("A", ancestries=[AncestryRecord("initial", "European", 100)])
    assert classify_study(eur, cfg.populations.mapping).population == "EUR"
    mixed = CandidateStudy("B", ancestries=[AncestryRecord("initial", "European", 100), AncestryRecord("initial", "East Asian", 100)])
    mixed = classify_study(mixed, cfg.populations.mapping)
    assert mixed.population == "MULTI"
    assert mixed.analysis_label == "COMBINED_PUBLISHED"


def test_embedded_ancestry_count_and_country_are_preserved() -> None:
    record = parse_ancestry_string("513266 European (U.S.)")[0]
    assert record.category == "European"
    assert record.n == 513266
    assert record.country == "U.S."


def test_aggregated_multiancestry_is_never_split() -> None:
    with pytest.raises(AggregatedSplitError):
        guard_population_label({"population": "MULTI", "analysis_label": "EUR"})
    assert guard_population_label({"population": "MULTI", "analysis_label": "COMBINED_PUBLISHED"}) == "COMBINED_PUBLISHED"


def test_or_to_log_or_z_and_effective_n() -> None:
    frame = pl.DataFrame({"CHR": ["1"], "POS": [2], "EA": ["A"], "NEA": ["G"], "OR": [2.0], "SE": [0.1], "P": [None], "N_CASES": [100.0], "N_CONTROLS": [300.0]}).lazy()
    row = derive_statistics(frame, {}).collect().row(0, named=True)
    assert row["BETA"] == pytest.approx(math.log(2.0))
    assert row["Z"] == pytest.approx(math.log(2.0) / 0.1)
    assert row["N_EFF"] == pytest.approx(300.0)
    assert effective_n(100, 300) == pytest.approx(300.0)


def test_quantitative_trait_uses_n_without_neff() -> None:
    frame = pl.DataFrame({"CHR": ["1"], "POS": [2], "EA": ["A"], "NEA": ["G"], "BETA": [0.2], "SE": [0.1], "N": [1234.0]}).lazy()
    row = derive_statistics(frame, {}).collect().row(0, named=True)
    assert row["N"] == 1234
    assert row["N_EFF"] is None


def test_meta_analysis_q_i2() -> None:
    long = pl.DataFrame(
        {
            "VARIANT_ID": ["1:1:A:G", "1:1:A:G"], "CHR": ["1", "1"], "POS": [1, 1],
            "A1": ["A", "A"], "A2": ["G", "G"], "RSID": ["rs1", "rs1"],
            "POP": ["AFR", "EUR"], "BETA": [0.1, 0.2], "SE": [0.1, 0.1],
            "EAF": [0.2, 0.3], "N": [1000.0, 2000.0], "N_CASES": [None, None], "N_CONTROLS": [None, None],
        }
    )
    row = inverse_variance_meta(long, ["AFR", "EUR"]).row(0, named=True)
    assert row["BETA"] == pytest.approx(0.15)
    assert row["SE"] == pytest.approx(math.sqrt(1 / 200))
    assert row["Q"] == pytest.approx(0.5)
    assert row["I2"] == pytest.approx(0.0)
    assert row["N_POPULATIONS"] == 2


def test_variant_normalization_and_gtex_key() -> None:
    frame = pl.DataFrame({"CHR": ["1"], "POS": [12], "REF": ["A"], "ALT": ["G"]}).with_columns(variant_id().alias("VARIANT_ID"))
    assert split_variant_id(frame)["POS"][0] == 12
    assert parse_gtex_variant("chr1_12_A_G_b38") == "1:12:A:G"


def test_allele_harmonization_orients_effect_to_reference_alt() -> None:
    gwas = pl.DataFrame(
        {
            "CHR": ["1"], "POS": [10], "RSID": ["rs1"], "EA": ["A"], "NEA": ["G"],
            "BETA": [0.2], "SE": [0.1], "P": [0.05], "EAF": [0.3], "N": [1000.0],
            "N_CASES": [None], "N_CONTROLS": [None], "N_EFF": [None],
        }
    )
    reference = pl.DataFrame({"CHR": ["1"], "POS": [10], "ID": ["1:10:A:G"], "REF": ["A"], "ALT": ["G"]})
    row = match_to_reference(gwas, reference).row(0, named=True)
    assert row["BETA_ALT"] == pytest.approx(-0.2)
    assert row["Z_ALT"] == pytest.approx(-2.0)
    assert row["EAF_ALT"] == pytest.approx(0.7)


def test_build_validation(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    variants = pl.DataFrame({"CHR": ["1", "1"], "POS": [1, 2], "EA": ["G", "G"], "NEA": ["A", "C"]})
    result = validate_build(variants, fasta, sample=100, min_concordance=0.9)
    assert result["status"] == "checked"
    assert result["ref_concordance"] == 1.0


def test_reference_selection_rejects_combined(tmp_path: Path) -> None:
    cfg = load_config()
    panel = ReferencePanel(cfg.reference, tmp_path, ToolResolver(cfg.tools))
    with pytest.raises(ValueError, match="pooled panels are not supported"):
        panel.pop_dir("MULTI")
    assert panel.pop_dir("AFR").name == "AFR"


def test_all_tissues_default_and_clpp() -> None:
    tissues = ["Liver", "Brain_Cortex", "Artery_Aorta"]
    assert requested_tissues(tissues) == sorted(tissues)
    values = variant_clpp(np.array([0.5, 0.2]), np.array([0.4, 0.5]))
    assert values.tolist() == pytest.approx([0.2, 0.1])
    assert locus_clpp(values) == pytest.approx(1 - (1 - 0.2) * (1 - 0.1))
