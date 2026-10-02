#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba is required." >&2
  exit 2
fi
LOCAL_ROOT="$ROOT_DIR/.gwas2m"
CORE_ENV="$LOCAL_ROOT/envs/gwas2mechanism"
export CONDA_PKGS_DIRS="$LOCAL_ROOT/pkgs"
export GWAS2M_ENV_ROOT="$LOCAL_ROOT/envs"
export GWAS2M_CACHE="$LOCAL_ROOT/resources"

if [[ ! -d "$CORE_ENV/conda-meta" ]]; then
  echo "Local environment not found. Run: bash Install.sh" >&2
  exit 2
fi

mamba run --prefix "$CORE_ENV" gwas2m setup --all --executor local --populations EUR --populations AFR --populations EAS --populations SAS --populations AMR "$@"
