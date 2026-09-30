from __future__ import annotations

import numpy as np
import polars as pl

from gwas2mechanism.clumping.loci import construct_loci, greedy_ld_clump
from gwas2mechanism.config import load_config
from gwas2mechanism.discovery.models import CandidateFile, CandidateStudy
from gwas2mechanism.discovery.ranking import build_candidate_table, score_candidates, select_studies


def test_ranking_selects_strongest_per_population() -> None:
    cfg = load_config()
    studies = [
        CandidateStudy("A", reported_trait="trait", population="EUR", analysis_label="EUR", n_total=1000, match_type="exact_mapped_trait", match_score=1.0, summary_stats_available=True),
        CandidateStudy("B", reported_trait="trait", population="EUR", analysis_label="EUR", n_total=10000, match_type="exact_mapped_trait", match_score=1.0, summary_stats_available=True),
    ]
    files = {
        study.study_accession: [CandidateFile(study.study_accession, f"file:///{study.study_accession}.tsv", f"{study.study_accession}.tsv", "harmonised", build="GRCh38", is_harmonised=True, has_beta=True, has_se=True, has_alleles=True)]
        for study in studies
    }
    table = select_studies(score_candidates(build_candidate_table(studies, files), cfg), cfg)
    assert table.filter(pl.col("selected"))["study_accession"].to_list() == ["B"]
    assert "selection_score" in table.columns


def test_ld_clump_and_locus_merge() -> None:
    frame = pl.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 150, 5000], "P": [1e-10, 1e-9, 1e-12], "VARIANT_ID": ["a", "b", "c"]})
    ld = np.eye(3)
    ld[0, 1] = ld[1, 0] = 0.9
    leads = greedy_ld_clump(frame, ld, kb=1, r2=0.1)
    assert set(leads["LEAD_VARIANT"]) == {"a", "c"}
    loci = construct_loci(leads, population="EUR", window_bp=100)
    assert loci.height == 2
