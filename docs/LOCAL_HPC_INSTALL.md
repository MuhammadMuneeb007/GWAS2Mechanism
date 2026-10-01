# Project-local HPC installation

The top-level installer can be downloaded into an empty working directory. It
clones the public repository and then completes the entire installation without
requiring a separate `git clone` command:

```bash
curl -fsSL \
  https://raw.githubusercontent.com/MuhammadMuneeb007/GWAS2Mechanism/main/Install.sh \
  -o Install.sh
bash Install.sh
cd GWAS2Mechanism
source scripts/activate.sh
```

The installer intentionally avoids `$HOME/.cache`, named Conda environments,
and shared `/data` paths.

`Install.sh` installs every environment, downloads resources in four labelled
stages, runs the test suite, and creates a synthetic report. The equivalent
lower-level command is `bash scripts/install.sh --full --example`.

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
bash Install.sh
```

For later shells, reactivate the project-local paths with:

```bash
cd /path/to/GWAS2Mechanism
source scripts/activate.sh
```
