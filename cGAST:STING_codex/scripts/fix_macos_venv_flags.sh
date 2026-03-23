#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="$ROOT_DIR/mtgdeckdash/.venv"
SITE_PACKAGES="$VENV_DIR/lib/python3.11/site-packages"
PYSIDE_DIR="$SITE_PACKAGES/PySide6"
SHIBOKEN_DIR="$SITE_PACKAGES/shiboken6"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This script is macOS-only."
  exit 0
fi

if [[ ! -d "$VENV_DIR" ]]; then
  echo "Missing venv: $VENV_DIR"
  echo "Create/install first, then rerun."
  exit 1
fi

if [[ -d "$SITE_PACKAGES" ]]; then
  if ls "$SITE_PACKAGES"/*.pth >/dev/null 2>&1; then
    chflags nohidden "$SITE_PACKAGES"/*.pth || true
  fi
fi

if [[ -d "$PYSIDE_DIR" ]]; then
  chflags -R nohidden "$PYSIDE_DIR" || true
fi

if [[ -d "$SHIBOKEN_DIR" ]]; then
  chflags -R nohidden "$SHIBOKEN_DIR" || true
fi

# Remove potentially problematic extended attributes and ad-hoc sign Qt binaries
# so Python can load them under stricter macOS policies.
xattr -cr "$PYSIDE_DIR" "$SHIBOKEN_DIR" 2>/dev/null || true

for base in "$PYSIDE_DIR" "$SHIBOKEN_DIR"; do
  if [[ ! -d "$base" ]]; then
    continue
  fi
  while IFS= read -r -d '' f; do
    if file -b "$f" 2>/dev/null | grep -q 'Mach-O'; then
      codesign --force --sign - "$f" >/dev/null 2>&1 || true
    fi
  done < <(find "$base" -type f -print0)
done

echo "macOS flags/xattrs fixed and PySide6 binaries ad-hoc signed."
