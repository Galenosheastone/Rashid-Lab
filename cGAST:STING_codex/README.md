# MTGDeckDash

Standalone MTG Commander deck analyzer that reads plain-text decklists, resolves card data from Scryfall, and outputs:

- `report.html` (interactive Plotly dashboard)
- `summary.json`
- `tagged_cards.json`
- `deck_features.json`
- `cards.csv`
- `deck_summary.csv`

Both a GUI (`PySide6`) and CLI are included.

## Features

- Robust decklist parser for MTGO/Moxfield-style lines like:
  - `1 Karn, Legacy Reforged (MAT) 49`
  - `19 Wastes (FIN) 309`
- Ignores comment/header lines starting with `//`
- Commander detection from `// COMMANDER` header, with override support
- Scryfall-only card resolution with local disk cache and retry/backoff handling
- YAML-configurable regex tagging engine (`mtgdeckdash/rules/default_rules.yaml`)
- Deterministic deck metrics + archetype inference (transparent heuristic scoring)
- Offline-friendly re-runs via cache

## Install and Run From Source

```bash
pip install -e .
```

CLI usage:

```bash
mtgdeckdash analyze /path/to/deck.txt --out /path/to/output --format commander
```

With commander override, custom rules, and auto-open report:

```bash
mtgdeckdash analyze /path/to/deck.txt --out /path/to/output --format commander --commander "Karn, Legacy Reforged" --rules /path/to/rules.yaml --open
```

Simulation usage (Path A, Forge via Docker):

```bash
mtgdeckdash simulate --deck-a /path/to/deck_a.txt --deck-b /path/to/deck_b.txt --games 200 --seed 1 --out /path/to/sim_output --open
```

Simulation usage (local heuristic backend, no Docker/GHCR required):

```bash
mtgdeckdash simulate --deck-a /path/to/deck_a.txt --deck-b /path/to/deck_b.txt --games 200 --sim-backend local
```

Simulation with explicit Forge Docker image override:

```bash
mtgdeckdash simulate --deck-a /path/to/deck_a.txt --deck-b /path/to/deck_b.txt --games 200 --forge-image <your/forge-image:tag>
```

Simulation using existing analyzed folders:

```bash
mtgdeckdash simulate --deck-a-out /path/to/deck_a_analysis --deck-b-out /path/to/deck_b_analysis --games 200
```

GUI usage:

```bash
mtgdeckdash-gui
```

GUI flow:

1. Choose a deck file or paste deck text.
2. Optionally set commander override and custom rules YAML.
3. Set output folder and click **Analyze**.
4. Open `report.html` and exports from the generated output folder.

## Build Standalone Executables (PyInstaller)

Install build dependency:

```bash
pip install pyinstaller
```

Build scripts:

- macOS: `scripts/build_mac.sh`
- Windows (PowerShell): `scripts/build_win.ps1`
- Linux: `scripts/build_linux.sh`

Each script builds:

- GUI app: `dist/MTGDeckDash` (`.app` bundle on macOS)
- CLI binary: `dist/mtgdeckdash-cli`

## Output Files

The output directory always contains:

- `report.html`
- `summary.json`
- `tagged_cards.json`
- `deck_features.json`
- `cards.csv`
- `deck_summary.csv`

Simulation output directory contains:

- `matchup_report.html`
- `matchup_summary.json`
- `raw_logs/stdout.log`
- `raw_logs/stderr.log`

## Cache Location

Scryfall cache root:

- macOS: `~/Library/Application Support/mtgdeckdash/cache/scryfall/`
- Windows: `%APPDATA%\mtgdeckdash\cache\scryfall\`
- Linux: `~/.local/share/mtgdeckdash/cache/scryfall/`

Analysis/simulation cache roots:

- analyses: `<user_data_dir>/analyses/<deck_id>/`
- simulations: `<user_data_dir>/simulations/<timestamp_hash>/`

## Troubleshooting

- **Scryfall 429/rate limits**: client retries with exponential backoff.
- **Network unavailable**: previously cached cards are reused; unresolved cards appear in error lists.
- **Some cards unresolved**: check spelling, set code, collector number, or remove printing suffix to allow exact-name lookup.
- **Commander warnings**: deck size, duplicate checks, and legality checks are best-effort based on resolved Scryfall data.
- **Docker unavailable for simulate**: install/start Docker Desktop and verify `docker info` works.
- **Need simulation without Docker/GHCR**: use `--sim-backend local` in CLI, or select `Local heuristic (no Docker)` in the GUI backend dropdown.
- **Forge image pull denied (GHCR/private image)**:
  - Authenticate to GHCR if needed:
    - `echo <GH_PAT> | docker login ghcr.io -u <github_user> --password-stdin`
  - Or point MTGDeckDash to an accessible image:
    - `export FORGE_DOCKER_IMAGE=<your/forge-image:tag>`
    - or use CLI `--forge-image <your/forge-image:tag>`
    - or set **Forge Image** in the GUI Simulation tab
- **Forge unresolved cards**: simulation surfaces unresolved names; re-run analyze to normalize names and retry.
- **`ModuleNotFoundError: No module named 'mtgdeckdash'` after editable install on macOS**: clear hidden flags on `.pth` files in the venv:
  - `chflags nohidden mtgdeckdash/.venv/lib/python3.11/site-packages/*.pth`
- **Qt plugin error (`Could not find the Qt platform plugin "cocoa"`) on macOS**: clear hidden flags on PySide6 plugin directories/files:
  - `./scripts/fix_macos_venv_flags.sh`

## Tests

```bash
pytest
```

Current unit tests cover parser behavior and default tagging rules.
