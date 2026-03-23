#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python -m PyInstaller \
  --noconfirm \
  --clean \
  --name MTGDeckDash \
  --windowed \
  --add-data "mtgdeckdash/rules/default_rules.yaml:mtgdeckdash/rules" \
  --add-data "mtgdeckdash/assets:mtgdeckdash/assets" \
  mtgdeckdash/gui.py

python -m PyInstaller \
  --noconfirm \
  --clean \
  --name mtgdeckdash-cli \
  --console \
  --add-data "mtgdeckdash/rules/default_rules.yaml:mtgdeckdash/rules" \
  --add-data "mtgdeckdash/assets:mtgdeckdash/assets" \
  mtgdeckdash/cli.py

echo "Build complete. Outputs are in dist/."
