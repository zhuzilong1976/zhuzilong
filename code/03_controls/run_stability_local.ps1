# =============================================================================
# Run the stability replicates on a local Windows machine (no cluster).
#
# Jobs = Replicates x 2 conditions. Each job is single-threaded, so we run
# $Parallel of them at a time: wall-clock ~ (jobs / Parallel) x per-job time.
#
# Memory note: one 800-gene job peaks at 1-2 GB. Choose $Parallel so that
#   $Parallel x 2 GB  +  other running analyses  <  available RAM.
#
# Usage:
#   powershell -File code/03_controls/run_stability_local.ps1 `
#       [-Parallel 3] [-Replicates 5] [-Genes 800] [-Cells 800] [-Rscript "Rscript"]
# =============================================================================
param(
  [int]$Parallel = 3,
  [int]$Replicates = 5,
  [int]$Genes = 800,
  [int]$Cells = 800,
  [string]$Rscript = "Rscript"
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path "logs" | Out-Null

$jobs = @()
foreach ($rep in 1..$Replicates) { $jobs += [pscustomobject]@{ Tag = "MSvC_Control"; Rep = $rep } }
foreach ($rep in 1..$Replicates) { $jobs += [pscustomobject]@{ Tag = "MSvC_MS";      Rep = $rep } }

$running = @()
foreach ($j in $jobs) {
  while (($running | Where-Object { -not $_.HasExited }).Count -ge $Parallel) {
    Start-Sleep -Seconds 20
    # @() is required: without it PowerShell unwraps a single surviving process
    # to a scalar and the subsequent "+=" fails (System.Diagnostics.Process has
    # no op_Addition).
    $running = @($running | Where-Object { -not $_.HasExited })
  }
  $log = "logs\stability_$($j.Tag)_rep$($j.Rep).log"
  $p = Start-Process -FilePath $Rscript `
       -ArgumentList "code/03_controls/28_stability_replicates.R", $j.Tag, $j.Rep, $Cells, $Genes `
       -RedirectStandardOutput $log -RedirectStandardError "$log.err" `
       -WindowStyle Hidden -PassThru
  Write-Host "[$(Get-Date -Format HH:mm:ss)] started $($j.Tag) rep $($j.Rep) (pid $($p.Id))"
  $running = @($running) + $p
}
$running | ForEach-Object { $_.WaitForExit() }
Write-Host "[$(Get-Date -Format HH:mm:ss)] all replicates finished"

# Then aggregate:
#   Rscript code/03_controls/29_stability_analysis.R
