#!/usr/bin/env bash

# Source this file from the repository root:
#   source scripts/activate.sh

GWAS2M_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export GWAS2M_PROJECT_ROOT
export GWAS2M_ENV_ROOT="$GWAS2M_PROJECT_ROOT/.gwas2m/envs"
export GWAS2M_CACHE="$GWAS2M_PROJECT_ROOT/.gwas2m/resources"
export CONDA_PKGS_DIRS="$GWAS2M_PROJECT_ROOT/.gwas2m/pkgs"
export PATH="$GWAS2M_ENV_ROOT/gwas2mechanism/bin:$PATH"

cd "$GWAS2M_PROJECT_ROOT"
echo "GWAS2Mechanism local environment active"
echo "  project:   $GWAS2M_PROJECT_ROOT"
echo "  resources: $GWAS2M_CACHE"
echo "  runs:      $GWAS2M_PROJECT_ROOT/runs"
