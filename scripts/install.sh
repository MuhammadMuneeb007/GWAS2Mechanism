#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_ROOT="$ROOT_DIR/.gwas2m"
ENV_ROOT="$LOCAL_ROOT/envs"
CORE_ENV="$ENV_ROOT/gwas2mechanism"
FULL=0
if [[ "${1:-}" == "--full" ]]; then FULL=1; fi

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
  if [[ -d "$prefix/conda-meta" ]]; then
    mamba env update --yes --prefix "$prefix" --file "$spec" --prune
  else
    mamba env create --yes --prefix "$prefix" --file "$spec" --override-channels
  fi
}

install_or_update "$CORE_ENV" environment.yml
if [[ "$FULL" -eq 1 ]]; then
  install_or_update "$ENV_ROOT/gwas2mechanism-spliceai" workflow/envs/spliceai.yml
  install_or_update "$ENV_ROOT/gwas2mechanism-pangolin" workflow/envs/pangolin.yml
  install_or_update "$ENV_ROOT/gwas2mechanism-vep" workflow/envs/vep.yml
  install_or_update "$ENV_ROOT/gwas2mechanism-r-finemap" workflow/envs/r-finemap.yml
  mamba run --prefix "$CORE_ENV" gwas2m setup --full
fi
mamba run --prefix "$CORE_ENV" gwas2m doctor

echo
echo "Local installation complete:"
echo "  core environment: $CORE_ENV"
echo "  package cache:    $LOCAL_ROOT/pkgs"
echo "  resources:        $GWAS2M_CACHE"
echo "  runs:             $ROOT_DIR/runs"
echo "Run: source scripts/activate.sh"
