#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="gwas2mechanism"
FULL=0
if [[ "${1:-}" == "--full" ]]; then FULL=1; fi

if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba is required. Install Miniforge/Mambaforge first." >&2
  exit 2
fi

cd "$ROOT_DIR"
mamba config set channel_priority strict
mamba env create --yes --file environment.yml
if [[ "$FULL" -eq 1 ]]; then
  mamba env create --yes --file workflow/envs/spliceai.yml
  mamba env create --yes --file workflow/envs/pangolin.yml
  mamba env create --yes --file workflow/envs/vep.yml
  mamba env create --yes --file workflow/envs/r-finemap.yml
  mamba run -n "$ENV_NAME" gwas2m setup --full
fi
mamba run -n "$ENV_NAME" gwas2m doctor
