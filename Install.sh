#!/usr/bin/env bash
set -Eeuo pipefail

# Standalone bootstrap. By default this installs software only; --full also
# downloads the large scientific resources through the modular local executor.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_URL="https://github.com/MuhammadMuneeb007/GWAS2Mechanism.git"

if [[ ! -f "$SCRIPT_DIR/pyproject.toml" || ! -d "$SCRIPT_DIR/src/gwas2mechanism" ]]; then
  command -v git >/dev/null 2>&1 || { echo "ERROR: git is required." >&2; exit 2; }
  CHECKOUT="$PWD/GWAS2Mechanism"
  echo "Repository:  $REPOSITORY_URL"
  echo "Destination: $CHECKOUT"
  if [[ -d "$CHECKOUT/.git" ]]; then
    git -C "$CHECKOUT" pull --ff-only origin main
  elif [[ -e "$CHECKOUT" ]]; then
    echo "ERROR: $CHECKOUT exists but is not a Git checkout." >&2
    exit 2
  else
    git clone --branch main --single-branch "$REPOSITORY_URL" "$CHECKOUT"
  fi
  exec bash "$CHECKOUT/Install.sh" "$@"
fi

cd "$SCRIPT_DIR"
exec bash scripts/install.sh "$@"
