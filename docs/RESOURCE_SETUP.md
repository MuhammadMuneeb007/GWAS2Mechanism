# Resource setup architecture

## Commands

```bash
gwas2m setup --status
gwas2m setup --reference --populations EUR --executor local
gwas2m setup --gtex
gwas2m setup --vep
gwas2m setup --splice
gwas2m setup --all --executor local
gwas2m setup --all --executor slurm --partition ascher
gwas2m setup --generate-slurm
```

`--status` is read-only. Local setup prints the absolute cache path at startup
and completion. At completion it reports elapsed time, files and bytes
downloaded, cache hits, resumed transfers, and complete/partial/missing counts.

## Execution graph

```text
environment check
├── GENCODE genome ────────────────┐
│                                  └── splicing
├── 1000G raw chromosome array
│   └── ALL-PGEN chromosome array
│       └── ancestry/chromosome array
├── GTEx (bounded object workers)
└── VEP cache

all selected terminal jobs ─────────── validation
```

GTEx does not wait for 1000 Genomes, and VEP does not wait for GTEx. SLURM
submission uses `sbatch --parsable` and `afterok` dependencies.

## Reference preparation

For every requested chromosome, the raw GRCh38 VCF and TBI are downloaded once
and converted once using canonical `CHR:POS:REF:ALT` IDs:

```text
raw VCF -> ALL/chrN.pgen,pvar,psam,pvar.parquet
        -> EUR/chrN.*
        -> AFR/chrN.*
        -> EAS/chrN.*
        -> SAS/chrN.*
        -> AMR/chrN.*
```

Only requested populations are derived. The ALL dataset is an implementation
intermediate and is never selected as an LD reference.

## Restart and concurrency safety

- Downloads use configurable aria2 retries and exponential backoff, then resume
  the same `.part` file with HTTPX if aria2 remains unsuccessful.
- HTTPX shows byte-level tqdm progress and resumes with HTTP Range requests.
- `.verified.json` reuses size, mtime, and checksum results without rehashing
  giant unchanged files.
- Atomic lock directories prevent two processes from writing the same `.part`,
  PGEN dataset, extraction directory, or manifest simultaneously.
- A chromosome is complete only when non-empty PGEN, PVAR, PSAM, and Parquet
  files match its atomic completion marker and have positive sample/variant
  counts.
- GTEx extraction rejects path traversal and archive links.

## Configuration

All worker and network caps live under `performance` in `config/default.yaml`:

```yaml
performance:
  downloads:
    concurrent_files: 4
    connections_per_file: 4
    min_split_size: 16M
    retry_attempts: 8
    retry_initial_seconds: 2
    retry_max_seconds: 60
  setup:
    chromosome_workers: 4
    population_workers: 8
    gtex_workers: 4
    slurm_download_array_limit: 6
    slurm_convert_array_limit: 8
    slurm_prepare_array_limit: 12
```

Override settings in a normal GWAS2Mechanism YAML config and pass it with
`--config`. Local CPU use also honours `SLURM_CPUS_PER_TASK`.

## Status example

```text
RESOURCE                         STATUS    DETAIL
GRCh38 FASTA                     COMPLETE
GENCODE 50 annotation            COMPLETE
GRCh38 FASTA index               COMPLETE
1000G raw                        COMPLETE  22/22
1000G ALL PGEN (intermediate)    COMPLETE  22/22
EUR                              COMPLETE  22/22
AFR                              PARTIAL   18/22
GTEx SuSiE eQTL                  COMPLETE
GTEx SuSiE sQTL                  COMPLETE
VEP cache                        COMPLETE
SpliceAI                         COMPLETE
Pangolin annotation DB           COMPLETE
Total: 15  Complete: 13  Partial: 1  Missing: 1  Percentage complete: 86.7%
```

## Storage

The main categories are raw 1000 Genomes VCFs, ALL PGEN intermediates,
ancestry-specific PGENs, GENCODE, compact GTEx SuSiE data, and VEP. Actual disk
use and duration depend on upstream releases, selected populations, network,
filesystem, and scheduler load; the package does not promise fixed values.
