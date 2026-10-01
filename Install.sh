#!/usr/bin/env bash
set -euo pipefail

# One-command project-local installation for Linux/HPC systems.
#
# Run from anywhere after cloning:
#   bash Install.sh
#
# It installs all environments under ./.gwas2m/envs, downloads all configured
# scientific resources under ./.gwas2m/resources, validates the installation,
# and creates a synthetic example report under ./runs.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

echo "GWAS2Mechanism complete local installation"
echo "Project:   $ROOT_DIR"
echo "Envs:      $ROOT_DIR/.gwas2m/envs"
echo "Packages:  $ROOT_DIR/.gwas2m/pkgs"
echo "Resources: $ROOT_DIR/.gwas2m/resources"
echo "Runs:      $ROOT_DIR/runs"
echo
echo "The full resource phase is large and may take hours. It is restartable."

bash "$ROOT_DIR/scripts/install.sh" --full --example

echo
echo "Installation and example completed successfully."
echo "For this and future shells run:"
echo "  cd '$ROOT_DIR'"
echo "  source scripts/activate.sh"
