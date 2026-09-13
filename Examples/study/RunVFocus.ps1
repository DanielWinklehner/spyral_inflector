# Vertical-focusing campaign: the shape scan (27 points, ~14 h) and then the final push for
# the winner (~2 h). Resumable: finished points and phases are reused on a relaunch.
#   powershell -File RunVFocus.ps1 [-Name vfocus1]
param([string]$Name = "vfocus1")
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$log = Join-Path $deck "Results\vfocus\log_$Name.txt"
New-Item -ItemType Directory -Force (Join-Path $deck "Results\vfocus") | Out-Null
"[$(Get-Date -Format HH:mm:ss)] launching VFocusScan --name $Name, scan only (pid $PID)" | Out-File -Append -Encoding ascii $log
$tol = 0.05
# scan only; RunVFocusFinal.ps1 waits for summary.json and does the re-ranking, refinement and final push
python -u VFocusScan.py --name $Name --sigma 0.0004 0.0012 0.0022 --angling 7 11 15 --gamma 5 8 11 --rank-tol $tol 2>&1 | Out-File -Append -Encoding utf8 $log
"[$(Get-Date -Format HH:mm:ss)] RunVFocus.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
