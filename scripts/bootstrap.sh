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
export CONDA_CHANNEL_PRIORITY=strict
export GWAS2M_ENV_ROOT="$LOCAL_ROOT/envs"
export GWAS2M_CACHE="$LOCAL_ROOT/resources"
mkdir -p "$LOCAL_ROOT/pkgs" "$LOCAL_ROOT/resources" "$ROOT_DIR/runs"

if [[ -d "$CORE_ENV/conda-meta" ]]; then
  mamba env update --yes --prefix "$CORE_ENV" --file environment.yml --prune
else
  mamba env create --yes --prefix "$CORE_ENV" --file environment.yml --override-channels
fi
mamba run --prefix "$CORE_ENV" python -m pip install --no-deps --editable .
