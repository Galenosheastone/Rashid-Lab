$ErrorActionPreference = "Stop"
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $Root

python -m PyInstaller `
  --noconfirm `
  --clean `
  --name MTGDeckDash `
  --windowed `
  --add-data "mtgdeckdash/rules/default_rules.yaml;mtgdeckdash/rules" `
  --add-data "mtgdeckdash/assets;mtgdeckdash/assets" `
  mtgdeckdash/gui.py

python -m PyInstaller `
  --noconfirm `
  --clean `
  --name mtgdeckdash-cli `
  --console `
  --add-data "mtgdeckdash/rules/default_rules.yaml;mtgdeckdash/rules" `
  --add-data "mtgdeckdash/assets;mtgdeckdash/assets" `
  mtgdeckdash/cli.py

Write-Host "Build complete. Outputs are in dist/."
