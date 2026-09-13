# Stage 5 of the vertical-focusing campaign: waits for the SC retune to finish (sc_retune.md),
# then levels the 8 mA bunch centroid with the exit truncation and the axial shift as knobs
# (SCLevel.py from the `triplet` worktree; the spiral voltage stays the centering knob), at the
# SC-retuned quad setting (sc_retune\sc_grid.json), with the full runs and hand-offs (--final).
#   powershell -File RunSCLevel.ps1 [-Run vfocus1_final] [-Package "<worktree path>"] [-BestJson sc_retune\sc_grid.json]
param([string]$Run = "vfocus1_final",
      [string]$Package = "D:\Dropbox (Personal)\Code\Python\spyral_inflector_triplet",
      [string]$BestJson = "sc_retune\sc_grid.json")
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
$rundir = Join-Path $deck "Results\final\$Run"
$log = Join-Path $deck "Results\vfocus\log_${Run}_sclevel.txt"
$env:PYTHONPATH = $Package
$env:SI_DECK = $deck
Set-Location (Join-Path $Package "Examples\study")
"[$(Get-Date -Format HH:mm:ss)] waiting for $rundir\sc_retune.md (pid $PID; package $Package)" | Out-File -Append -Encoding ascii $log
while (-not (Test-Path (Join-Path $rundir "sc_retune.md"))) { Start-Sleep -Seconds 60 }
"[$(Get-Date -Format HH:mm:ss)] SC retune finished; levelling $Run at the quad setting of $BestJson" | Out-File -Append -Encoding ascii $log
python -u SCLevel.py --run $rundir --base-run $rundir --best-json $BestJson --quads 16 24 --final 2>&1 | Out-File -Append -Encoding utf8 $log
"[$(Get-Date -Format HH:mm:ss)] RunSCLevel.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
