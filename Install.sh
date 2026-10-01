#!/usr/bin/env bash
set -Eeuo pipefail

# Complete, standalone Linux/HPC bootstrap and installer.
#
# When downloaded into an empty directory, this file first clones the public
# GitHub repository into ./GWAS2Mechanism and then continues with the copy of
# itself inside that checkout. When run inside a checkout, it proceeds directly
# to installation. It never deletes or overwrites a non-Git directory.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_URL="https://github.com/MuhammadMuneeb007/GWAS2Mechanism.git"

if [[ ! -f "$SCRIPT_DIR/pyproject.toml" || ! -d "$SCRIPT_DIR/src/gwas2mechanism" ]]; then
  if ! command -v git >/dev/null 2>&1; then
    echo "ERROR: git is required to download GWAS2Mechanism." >&2
    exit 2
  fi

  INSTALL_PARENT="$PWD"
  CHECKOUT="$INSTALL_PARENT/GWAS2Mechanism"

  echo "GWAS2Mechanism bootstrap installer"
  echo "Repository: $REPOSITORY_URL"
  echo "Destination: $CHECKOUT"

  if [[ -d "$CHECKOUT/.git" ]]; then
    echo "An existing Git checkout was found; updating it safely."
    git -C "$CHECKOUT" pull --ff-only origin main
  elif [[ -e "$CHECKOUT" ]]; then
    echo "ERROR: $CHECKOUT exists but is not a Git checkout." >&2
    echo "Move it aside or run this installer from another empty directory." >&2
    exit 2
  else
    git clone --branch main --single-branch "$REPOSITORY_URL" "$CHECKOUT"
  fi

  echo "Repository downloaded. Continuing with $CHECKOUT/Install.sh"
  exec bash "$CHECKOUT/Install.sh"
fi

ROOT_DIR="$SCRIPT_DIR"
LOCAL_ROOT="$ROOT_DIR/.gwas2m"
ENV_ROOT="$LOCAL_ROOT/envs"
CORE_ENV="$ENV_ROOT/gwas2mechanism"
RESOURCE_ROOT="$LOCAL_ROOT/resources"
RUN_ROOT="$ROOT_DIR/runs"

on_error() {
  local status=$?
  echo >&2
  echo "Installation failed at line $1 with exit status $status." >&2
  echo "Fix the reported problem and rerun: bash Install.sh" >&2
  echo "Completed environments and downloads will be reused." >&2
  exit "$status"
}
trap 'on_error "$LINENO"' ERR

heading() {
  echo
  echo "======================================================================"
  echo "$1"
  echo "======================================================================"
}

if ! command -v mamba >/dev/null 2>&1; then
  echo "ERROR: mamba is required and was not found on PATH." >&2
  echo "Load your HPC Miniforge/Conda module, then rerun bash Install.sh." >&2
  exit 2
fi

if ! command -v sha256sum >/dev/null 2>&1; then
  echo "ERROR: sha256sum is required and was not found on PATH." >&2
  exit 2
fi

cd "$ROOT_DIR"
mkdir -p "$ENV_ROOT" "$LOCAL_ROOT/pkgs" "$RESOURCE_ROOT" "$RUN_ROOT"

# Keep packages, environments, resources, and outputs inside this checkout.
export CONDA_PKGS_DIRS="$LOCAL_ROOT/pkgs"
export CONDA_CHANNEL_PRIORITY=strict
export GWAS2M_PROJECT_ROOT="$ROOT_DIR"
export GWAS2M_ENV_ROOT="$ENV_ROOT"
export GWAS2M_CACHE="$RESOURCE_ROOT"

echo "GWAS2Mechanism complete project-local installation"
echo "Project:   $ROOT_DIR"
echo "Envs:      $ENV_ROOT"
echo "Packages:  $CONDA_PKGS_DIRS"
echo "Resources: $GWAS2M_CACHE"
echo "Runs:      $RUN_ROOT"
echo "Mamba:     $(command -v mamba)"
echo
df -h "$ROOT_DIR" || true
echo
echo "The full genomic resource phase is large and may take hours."
echo "The installer is restartable; rerunning it reuses completed work."

install_environment() {
  local label="$1"
  local prefix="$2"
  local specification="$3"
  local marker="$prefix/.gwas2m-spec.sha256"
  local digest
  digest="$(sha256sum "$specification" | awk '{print $1}')"

  heading "Installing environment: $label"
  echo "Prefix: $prefix"
  echo "Specification: $specification"

  if [[ -d "$prefix/conda-meta" && -f "$marker" && "$(<"$marker")" == "$digest" ]]; then
    echo "Status: already installed with the current specification; reusing it."
    return
  fi

  if [[ -d "$prefix/conda-meta" ]]; then
    echo "Status: updating existing environment."
    mamba env update \
      --yes \
      --prefix "$prefix" \
      --file "$specification" \
      --prune
  else
    echo "Status: creating new environment."
    mamba env create \
      --yes \
      --prefix "$prefix" \
      --file "$specification" \
      --override-channels
  fi

  printf '%s\n' "$digest" > "$marker"
  echo "Completed: $label"
}

# 1. Install every environment.
install_environment \
  "Core GWAS2Mechanism, Snakemake, PLINK2, bcftools, R, SuSiE and coloc" \
  "$CORE_ENV" \
  "$ROOT_DIR/environment.yml"

install_environment \
  "SpliceAI" \
  "$ENV_ROOT/gwas2mechanism-spliceai" \
  "$ROOT_DIR/workflow/envs/spliceai.yml"

install_environment \
  "Pangolin" \
  "$ENV_ROOT/gwas2mechanism-pangolin" \
  "$ROOT_DIR/workflow/envs/pangolin.yml"

install_environment \
  "Ensembl VEP" \
  "$ENV_ROOT/gwas2mechanism-vep" \
  "$ROOT_DIR/workflow/envs/vep.yml"

install_environment \
  "R fine-mapping" \
  "$ENV_ROOT/gwas2mechanism-r-finemap" \
  "$ROOT_DIR/workflow/envs/r-finemap.yml"

# 2. Download and prepare every configured scientific resource. Each command
# is separate so the terminal identifies the active long-running phase.
heading "Resource stage 1/4: GRCh38 and 1000 Genomes reference panels"
echo "Downloading chromosomes and preparing EUR, AFR, EAS, SAS and AMR panels."
echo "This is normally the longest stage."
mamba run --prefix "$CORE_ENV" gwas2m setup \
  --reference \
  --populations EUR \
  --populations AFR \
  --populations EAS \
  --populations SAS \
  --populations AMR

heading "Resource stage 2/4: Adult GTEx all-tissue SuSiE resources"
mamba run --prefix "$CORE_ENV" gwas2m setup --gtex

heading "Resource stage 3/4: Ensembl VEP cache"
mamba run --prefix "$CORE_ENV" gwas2m setup --vep

heading "Resource stage 4/4: SpliceAI and Pangolin annotations"
mamba run --prefix "$CORE_ENV" gwas2m setup --splice

# 3. Validate the resulting installation.
heading "Validation stage 1/3: environment and executable diagnostics"
mamba run --prefix "$CORE_ENV" gwas2m doctor

heading "Validation stage 2/3: automated test suite"
mamba run --prefix "$CORE_ENV" pytest -q

# 4. Run a complete, deterministic example. It intentionally uses synthetic
# fixtures so validation does not depend on a particular live GWAS study.
heading "Validation stage 3/3: synthetic end-to-end example"
mamba run --prefix "$CORE_ENV" gwas2m run \
  --phenotype "HPC installation example" \
  --populations auto \
  --tissues all \
  --mode fast \
  --threads "${SLURM_CPUS_PER_TASK:-4}" \
  --synthetic

EXAMPLE_REPORT="$RUN_ROOT/hpc_installation_example/synthetic-smoke/report/report.html"
if [[ ! -s "$EXAMPLE_REPORT" ]]; then
  echo "ERROR: the expected example report was not created: $EXAMPLE_REPORT" >&2
  exit 1
fi

heading "Installation completed successfully"
echo "Core environment: $CORE_ENV"
echo "Resource cache:   $RESOURCE_ROOT"
echo "Example report:   $EXAMPLE_REPORT"
echo
echo "Activate the project in this and future shells with:"
echo "  cd '$ROOT_DIR'"
echo "  source scripts/activate.sh"
