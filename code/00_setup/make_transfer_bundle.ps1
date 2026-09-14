# =============================================================================
# Assemble a self-contained bundle for moving the project to another machine.
#
#   powershell -File code/00_setup/make_transfer_bundle.ps1 [-Dest D:\transfer] [-Full]
#
# The layout produced is the one the scripts expect: code/, docs/, results/ and
# work/ at the same level, so that the relative paths inside the scripts resolve
# without editing.
#
#   -Full  also includes the large regenerable items (CellOracle results 282 MB,
#          downloaded published networks 108 MB, raw h5ad 128 MB).
# =============================================================================
param(
  [string]$Dest = "transfer",
  [switch]$Full
)

$ErrorActionPreference = "Stop"
$root = Join-Path $Dest "vko-project"
New-Item -ItemType Directory -Force -Path $root | Out-Null

function Copy-Set([string]$src, [string]$rel) {
  if (-not (Test-Path $src)) { Write-Host "skip (missing): $src"; return }
  $dst = Join-Path $root $rel
  New-Item -ItemType Directory -Force -Path (Split-Path $dst -Parent) | Out-Null
  Copy-Item -Recurse -Force -Path $src -Destination $dst
  $sz = (Get-ChildItem $dst -Recurse -File | Measure-Object Length -Sum).Sum / 1MB
  Write-Host ("copied {0,-45} {1,8:N1} MB" -f $rel, $sz)
}

Write-Host "=== code, documentation and results (always) ==="
Copy-Set "outputs/repo/code"    "code"
Copy-Set "outputs/repo/docs"    "docs"
Copy-Set "outputs/repo/env"     "env"
Copy-Set "outputs/repo/results" "results"
Copy-Set "outputs/repo/data"    "data"
foreach ($f in "README.md", "LICENSE", "CITATION.cff", ".gitignore") {
  if (Test-Path "outputs/repo/$f") { Copy-Item "outputs/repo/$f" -Destination $root -Force }
}

Write-Host "=== derived data and completed analyses (always) ==="
Copy-Set "work/data"                            "work/data"
Copy-Set "work/stability"                       "work/stability"
Copy-Set "outputs/results/vko_ms_microglia"     "outputs/results/vko_ms_microglia"
Copy-Set "outputs/results/technical_null"       "outputs/results/technical_null"
Copy-Set "outputs/results/param_sensitivity"    "outputs/results/param_sensitivity"
Copy-Set "outputs/results/replication"          "outputs/results/replication"
Copy-Set "outputs/results/vko_ms_depthmatched"  "outputs/results/vko_ms_depthmatched"

if ($Full) {
  Write-Host "=== large regenerable items (-Full) ==="
  Copy-Set "work/geo"                "work/geo"
  Copy-Set "work/author_networks"    "work/author_networks"
  Copy-Set "outputs/results/celloracle" "outputs/results/celloracle"
}

$total = (Get-ChildItem $root -Recurse -File | Measure-Object Length -Sum).Sum / 1MB
Write-Host ""
Write-Host ("bundle: {0}  ({1:N0} MB, {2} files)" -f (Resolve-Path $root),
            $total, (Get-ChildItem $root -Recurse -File).Count)
Write-Host "compress with: Compress-Archive -Path '$root' -DestinationPath '$Dest\vko-project.zip'"
