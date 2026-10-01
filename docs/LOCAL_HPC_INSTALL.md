# Project-local HPC installation

Run every command from the `GWAS2Mechanism` checkout. The installer intentionally
avoids `$HOME/.cache`, named Conda environments, and shared `/data` paths.

```bash
git pull --ff-only origin main
bash scripts/install.sh --full
source scripts/activate.sh
pytest -q
gwas2m run --phenotype "HPC smoke test" --populations auto --tissues all \
  --mode fast --threads "${SLURM_CPUS_PER_TASK:-4}" --synthetic
```

The complete local layout is:

```text
GWAS2Mechanism/
├── .gwas2m/
│   ├── envs/
│   │   ├── gwas2mechanism/
│   │   ├── gwas2mechanism-vep/
│   │   ├── gwas2mechanism-spliceai/
│   │   ├── gwas2mechanism-pangolin/
│   │   └── gwas2mechanism-r-finemap/
│   ├── pkgs/
│   └── resources/
└── runs/
```

Do not pass a site-wide path such as
`--set run_root=/data/gwas2mechanism/runs` unless that exact directory is
writable. The default `runs` directory is already inside the checkout.

The installation and downloads are restartable. If a network transfer or
login session stops, rerun:

```bash
bash scripts/install.sh --full
```

For later shells, reactivate the project-local paths with:

```bash
cd /path/to/GWAS2Mechanism
source scripts/activate.sh
```
