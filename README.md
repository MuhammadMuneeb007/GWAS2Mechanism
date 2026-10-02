# GWAS2Mechanism

GWAS2Mechanism is a phenotype-agnostic, population-aware framework that
discovers GWAS summary statistics and integrates ancestry-specific and
cross-ancestry fine-mapping with VEP annotation, Adult GTEx molecular-QTL
evidence, SpliceAI, and Pangolin.

The `gwas2m` CLI uses Polars/Parquet, PLINK2, SuSiE-RSS, MultiSuSiE or SuSiEx,
compact GTEx SuSiE eQTL/sQTL resources, and a restartable Snakemake workflow.

## Installation

Mamba is the supported solver. Environments, package cache, and scientific
resources remain inside the checkout under `.gwas2m/`.

From an empty Linux/HPC directory:

```bash
curl -fsSL \
  https://raw.githubusercontent.com/MuhammadMuneeb007/GWAS2Mechanism/main/Install.sh \
  -o Install.sh
bash Install.sh
cd GWAS2Mechanism
source scripts/activate.sh
```

`bash Install.sh` clones the public repository when necessary, creates the
project-local software environments, installs the package, and runs software
diagnostics. It deliberately does not start the large scientific downloads.
Use `bash Install.sh --full` only when software plus full local resources are
wanted in one command.

## Resource setup

Read-only status:

```bash
gwas2m setup --status
```

Bounded local setup:

```bash
gwas2m setup --reference --populations EUR --executor local
gwas2m setup --all --executor local
```

SLURM setup (generates and submits the dependency DAG):

```bash
gwas2m setup --all --executor slurm --partition ascher
```

Generate the scripts without submitting:

```bash
gwas2m setup --generate-slurm
```

Jobs are written to `.gwas2m/setup_jobs/` and stdout/stderr to
`.gwas2m/setup_logs/`. No partition is hard-coded. GTEx and VEP are independent
of the reference pipeline. 1000 Genomes download, ALL-sample conversion, and
population subsetting use bounded chromosome arrays and `afterok` dependencies.

Resources are cached under `.gwas2m/resources/` unless `GWAS2M_CACHE` is set.
Major storage categories are compressed 1000 Genomes VCFs, one-time ALL-sample
PGEN intermediates, separate ancestry PGENs, GENCODE, compact GTEx SuSiE data,
and VEP. Sizes vary by release and requested populations.

Downloads retain `.part` data for resume, atomically rename completed files,
reuse `.verified.json` stamps, and use writer locks. Derived chromosome markers
validate PGEN/PVAR/PSAM/Parquet together. Rerunning local or SLURM setup skips
valid work and resumes partial work.

Raw VCFs are converted once per chromosome to an ALL-sample PGEN intermediate.
EUR, AFR, EAS, SAS, and AMR are then independently subset from that binary
intermediate. ALL is never used as pooled cross-ancestry LD.

See [resource setup documentation](docs/RESOURCE_SETUP.md) and the
[local/HPC installation guide](docs/LOCAL_HPC_INSTALL.md).

## Quick start

```bash
gwas2m discover --phenotype "migraine"
gwas2m run --phenotype "migraine" --populations auto --tissues all --mode fast --threads 32
gwas2m resume --run runs/migraine/latest
gwas2m report --run runs/migraine/latest
```

A download-free validation run is available:

```bash
gwas2m run --phenotype "synthetic phenotype" --synthetic
```

## Population strategy

The framework maps supported metadata to EUR, AFR, EAS, SAS, or AMR while
retaining reported ancestry and flagging approximate mappings. An aggregated
multi-ancestry file is never split into artificial ancestry datasets; it is
analysed as `COMBINED_PUBLISHED`. `COMBINED_META` is created only from compatible
ancestry-stratified effects. Cross-population fine-mapping supplies a separate
summary-statistic vector, sample size, and LD matrix for every ancestry.

## Adult GTEx strategy

`qtl.tissues: all` is the default. Every tissue in the selected stable Adult
GTEx compact SuSiE release is eligible for eQTL and sQTL matching. Tissues may
be ranked in reports but are not pre-filtered by phenotype.

## Outputs and interpretation

Each run is isolated under `runs/<phenotype>/<run_id>/` with its resolved
configuration, provenance manifest, discovery tables, population-specific
GWAS/loci/fine-mapping directories, multi-ancestry results, annotations, QTL,
splicing, colocalisation, evidence, and reports.

Summary statistics cannot be retrospectively split by ancestry. 1000 Genomes
may not perfectly match a study cohort. Missing evidence is not biological
absence. VEP, SpliceAI, and Pangolin provide predictions rather than
experimental validation. See [scientific limitations](docs/SCIENTIFIC_LIMITATIONS.md).

SpliceAI source/model licensing restricts model use to qualifying
non-commercial purposes unless a commercial license is obtained; review the
upstream license before using that optional stage.

## Development and citation

Run `ruff check .` and `pytest -q`; see [CONTRIBUTING.md](CONTRIBUTING.md).
Cite this software using [CITATION.cff](CITATION.cff), plus the upstream GWAS,
reference, GTEx, VEP, fine-mapping, and splicing sources recorded in manifests.
