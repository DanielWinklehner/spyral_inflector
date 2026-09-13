# Second half of the vertical-focusing campaign: waits for the shape scan to finish
# (summary.json), re-ranks the coarse grids, refines the top shapes with a fine q2 grid,
# then runs the final push for the refined winner with fine quad grids of its own.
#   powershell -File RunVFocusFinal.ps1 [-Name vfocus1] [-Tol 0.05] [-Top 4]
param([string]$Name = "vfocus1", [double]$Tol = 0.05, [int]$Top = 4)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$out = Join-Path $deck "Results\vfocus\$Name"
$log = Join-Path $deck "Results\vfocus\log_${Name}_final.txt"
"[$(Get-Date -Format HH:mm:ss)] waiting for $out\summary.json (pid $PID)" | Out-File -Append -Encoding ascii $log
while (-not (Test-Path (Join-Path $out "summary.json"))) { Start-Sleep -Seconds 60 }
"[$(Get-Date -Format HH:mm:ss)] scan finished; re-ranking and refining the top $Top at tolerance $Tol" | Out-File -Append -Encoding ascii $log
python -u VFocusRank.py --name $Name --tol 0.03 0.05 0.08 --write 2>&1 | Out-File -Append -Encoding utf8 $log
python -u VFocusRefine.py --name $Name --top $Top --tol $Tol 2>&1 | Out-File -Append -Encoding utf8 $log
$knobs = Join-Path $out "winner_knobs_refined.json"
if (-not (Test-Path $knobs)) { "[$(Get-Date -Format HH:mm:ss)] no refined winner; stopping" | Out-File -Append -Encoding ascii $log; exit 1 }
"[$(Get-Date -Format HH:mm:ss)] final push ${Name}_final with $(Get-Content $knobs -Raw)" | Out-File -Append -Encoding ascii $log
python -u RunFinalPush.py --name "${Name}_final" --bfield (Join-Path $deck "Fields\final1_baseline_1mm.pickle") --knobs-json $knobs --free-exit-truncation --rank vert --rank-tol $Tol --no-q2-exit-plate --q1 5575 7575 5 --q2 -7450 -5450 5 --fine-step 250 2>&1 | Out-File -Append -Encoding utf8 $log
"[$(Get-Date -Format HH:mm:ss)] RunVFocusFinal.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
