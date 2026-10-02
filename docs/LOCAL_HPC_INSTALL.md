# Project-local Linux/HPC installation

Run these commands in the directory where `GWAS2Mechanism/` should be created:

```bash
git clone https://github.com/MuhammadMuneeb007/GWAS2Mechanism.git
cd GWAS2Mechanism
bash Install.sh
source scripts/activate.sh
gwas2m doctor
gwas2m setup --status
```

The default install creates software only. It does not download the large
scientific datasets. Everything stays in the current checkout:

```text
GWAS2Mechanism/
├── .gwas2m/
│   ├── envs/
│   ├── pkgs/
│   ├── resources/
│   ├── setup_jobs/
│   └── setup_logs/
└── runs/
```

For a workstation or an interactive allocation:

```bash
gwas2m setup --all --executor local
```

For SLURM:

```bash
gwas2m setup --all --executor slurm --partition ascher
```

The partition is an example supplied explicitly by the user. Generate scripts
without submitting them with `gwas2m setup --generate-slurm`.

To prepare only one ancestry:

```bash
gwas2m setup --reference --populations EUR --executor slurm
```

The same commands are safe to rerun after interruption. Completed downloads,
chromosome conversions, and population subsets are validated and skipped.
Use `GWAS2M_CACHE=/writable/shared/path` before activation only when a shared
cache is intentional; otherwise resources stay in `.gwas2m/resources/`.
