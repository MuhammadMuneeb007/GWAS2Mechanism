# Scientific limitations and interpretation

GWAS2Mechanism is an evidence-integration framework, not a causal-proof
machine. Its outputs are ranked hypotheses that require study-specific review
and, where appropriate, experimental validation.

## Population and study design

- An aggregated multi-ancestry summary-statistics file cannot be separated
  into ancestry-specific datasets. Such inputs remain
  `COMBINED_PUBLISHED`.
- `COMBINED_META` is created only from independently reported,
  ancestry-stratified effects after allele alignment. Cochran's Q and I² are
  retained because an average effect can hide real heterogeneity.
- Reference-panel labels are broad approximations of cohort ancestry. The
  1000 Genomes panel may be a poor LD match for admixed, founder, or otherwise
  underrepresented cohorts.
- Smaller studies have less power. The systematic imbalance in available
  sample sizes means missing evidence must not be interpreted as biological
  absence.
- Cross-ancestry fine-mapping receives a separate Z vector, sample size, and LD
  matrix for each eligible population. It does not use pooled LD.

## Statistical assumptions

- SuSiE-RSS assumes well-aligned summary statistics and LD from an appropriate
  population. Allele mismatches, duplicated variants, imputation differences,
  or LD/sample mismatch can distort credible sets.
- The deterministic NumPy fine-mapper is a smoke-test fixture only. Scientific
  runs use `susieR`; reports and diagnostics state which engine produced each
  result.
- The PIP-product shared-variant score is a scalable screen, not formal
  colocalisation. Formal `coloc.susie` is attempted only when compatible dense
  regional inputs and fitted SuSiE objects are available.
- Multiple-testing thresholds, locus windows, priors, and credible-set coverage
  are configurable analytical choices rather than universally correct values.

## Functional evidence

- VEP consequence labels, QTL fine-mapping, SpliceAI, and Pangolin are
  predictive evidence. They do not establish molecular mediation or clinical
  pathogenicity.
- Gene-body overlap is reported separately from coding-sequence overlap. A
  variant within a lncRNA body is not thereby shown to regulate that lncRNA.
- GTEx evidence is limited by tissues, donors, ancestry composition, assay,
  expression level, and statistical power. A top tissue is not necessarily the
  causal tissue.
- Variant normalization and genome build must agree across every source.
  Liftover and allele harmonisation can introduce ambiguity; unresolved records
  are excluded and counted.
- SpliceAI weights and software are subject to upstream licensing constraints.
  Users must confirm that their intended use is permitted.

Use calibrated language: **candidate causal variant**, **credible-set member**,
**colocalisation evidence**, **predicted splice effect**, and **candidate
regulatory mechanism**.
