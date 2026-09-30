#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba is required." >&2
  exit 2
fi
mamba config set channel_priority strict
mamba env update --name gwas2mechanism --file environment.yml --prune
mamba run -n gwas2mechanism python -m pip install --no-deps --editable .
