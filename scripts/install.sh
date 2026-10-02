#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_ROOT="$ROOT_DIR/.gwas2m"
ENV_ROOT="$LOCAL_ROOT/envs"
CORE_ENV="$ENV_ROOT/gwas2mechanism"
FULL=0
EXAMPLE=0

usage() {
  cat <<'EOF'
Usage: bash scripts/install.sh [--full] [--example]

  Default    Install every software environment; do not download large resources.
  --full     Also download/prepare all resources with the bounded local executor.
  --example  Run tests and a synthetic end-to-end example after installation.
EOF
}

for argument in "$@"; do
  case "$argument" in
    --full) FULL=1 ;;
    --example) EXAMPLE=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $argument" >&2; usage >&2; exit 2 ;;
  esac
done

if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba is required. Install Miniforge/Mambaforge first." >&2
  exit 2
fi

cd "$ROOT_DIR"
mkdir -p "$ENV_ROOT" "$LOCAL_ROOT/pkgs" "$LOCAL_ROOT/resources" "$ROOT_DIR/runs"

# Keep environments, downloaded conda packages, and scientific resources in
# this checkout.  This is important on HPC systems where $HOME is quota-limited
# and shared /data paths may not be writable.
export CONDA_PKGS_DIRS="$LOCAL_ROOT/pkgs"
export CONDA_CHANNEL_PRIORITY=strict
export GWAS2M_ENV_ROOT="$ENV_ROOT"
export GWAS2M_CACHE="$LOCAL_ROOT/resources"

install_or_update() {
  local prefix="$1"
  local spec="$2"
  local marker="$prefix/.gwas2m-spec.sha256"
  local digest
  digest="$(sha256sum "$spec" | awk '{print $1}')"
  echo
  echo "======================================================================"
  echo "Environment: $prefix"
  echo "Specification: $spec"
  echo "======================================================================"
  if [[ -d "$prefix/conda-meta" && -f "$marker" && "$(<"$marker")" == "$digest" ]]; then
    echo "Already installed with the current specification; reusing it."
    return
  fi
  if [[ -d "$prefix/conda-meta" ]]; then
    mamba env update --yes --prefix "$prefix" --file "$spec" --prune
  else
    mamba env create --yes --prefix "$prefix" --file "$spec" --override-channels
  fi
  printf '%s\n' "$digest" > "$marker"
}

install_or_update "$CORE_ENV" environment.yml
install_or_update "$ENV_ROOT/gwas2mechanism-spliceai" workflow/envs/spliceai.yml
install_or_update "$ENV_ROOT/gwas2mechanism-pangolin" workflow/envs/pangolin.yml
install_or_update "$ENV_ROOT/gwas2mechanism-vep" workflow/envs/vep.yml
install_or_update "$ENV_ROOT/gwas2mechanism-r-finemap" workflow/envs/r-finemap.yml

if [[ "$FULL" -eq 1 ]]; then
  echo
  echo "======================================================================"
  echo "Full bounded local resource setup"
  echo "This can take hours and substantial disk space on an HPC filesystem."
  echo "======================================================================"
  mamba run --prefix "$CORE_ENV" gwas2m setup --all --executor local
fi
mamba run --prefix "$CORE_ENV" gwas2m doctor

if [[ "$EXAMPLE" -eq 1 ]]; then
  echo
  echo "======================================================================"
  echo "Validation: test suite"
  echo "======================================================================"
  mamba run --prefix "$CORE_ENV" pytest -q

  echo
  echo "======================================================================"
  echo "Example: complete download-free synthetic workflow"
  echo "======================================================================"
  mamba run --prefix "$CORE_ENV" gwas2m run \
    --phenotype "HPC installation example" \
    --populations auto --tissues all --mode fast \
    --threads "${SLURM_CPUS_PER_TASK:-4}" --synthetic
  test -s "$ROOT_DIR/runs/hpc_installation_example/synthetic-smoke/report/report.html"
  echo "Example report: $ROOT_DIR/runs/hpc_installation_example/synthetic-smoke/report/report.html"
fi

echo
echo "Local installation complete:"
echo "  core environment: $CORE_ENV"
echo "  package cache:    $LOCAL_ROOT/pkgs"
echo "  resources:        $GWAS2M_CACHE"
echo "  runs:             $ROOT_DIR/runs"
echo "Run: source scripts/activate.sh"
echo "Inspect resources: gwas2m setup --status"
echo "Local resource setup: gwas2m setup --all --executor local"
echo "SLURM resource setup: gwas2m setup --all --executor slurm --partition PARTITION"
