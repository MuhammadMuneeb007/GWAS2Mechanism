# Contributing

Use a feature branch, add tests for scientific behavior, and run `ruff check .` plus `pytest -q` before opening a pull request. Core code must remain phenotype-agnostic. Any ancestry mapping must preserve the source metadata and state whether the mapping is approximate. Never add downloaded individual- or summary-level genetic data, credentials, model weights, or reference caches to Git.

Changes to harmonisation, LD selection, meta-analysis, fine-mapping, QTL matching, or colocalisation require a short scientific rationale and a regression test. Outputs must distinguish unavailable evidence from evidence of absence.
