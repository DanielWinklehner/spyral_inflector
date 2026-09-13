# Stage 4 of the vertical-focusing campaign: waits for the final push to finish (phase_11.done),
# then re-tunes the quads with space charge in the loop (SCRetuneV.py) and re-runs the final
# space-charge run at that setting.
#   powershell -File RunSCRetune.ps1 [-Run vfocus1_final] [-Tol 0.05]
param([string]$Run = "vfocus1_final", [double]$Tol = 0.05)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$rundir = Join-Path $deck "Results\final\$Run"
$log = Join-Path $deck "Results\vfocus\log_${Run}_scretune.txt"
"[$(Get-Date -Format HH:mm:ss)] waiting for $rundir\phase_11.done (pid $PID)" | Out-File -Append -Encoding ascii $log
while (-not (Test-Path (Join-Path $rundir "phase_11.done"))) { Start-Sleep -Seconds 60 }
"[$(Get-Date -Format HH:mm:ss)] final push finished; SC retune at tolerance $Tol" | Out-File -Append -Encoding ascii $log
python -u SCRetuneV.py --run $rundir --tol $Tol 2>&1 | Out-File -Append -Encoding utf8 $log
"[$(Get-Date -Format HH:mm:ss)] RunSCRetune.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
