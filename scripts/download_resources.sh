#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba is required." >&2
  exit 2
fi
mamba run -n gwas2mechanism gwas2m setup --all --populations EUR --populations AFR --populations EAS --populations SAS --populations AMR "$@"
