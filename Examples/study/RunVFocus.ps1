# Vertical-focusing campaign: the shape scan (sigma x angling x gamma points, ~17 min each), scan
# only; RunVFocusFinal.ps1 waits for summary.json and does the re-ranking, refinement and final
# push. Resumable: finished points are reused on a relaunch.
#   powershell -File RunVFocus.ps1 [-Name vfocus1] [-Sigma "0.0004 0.0012 0.0022"] [-Angling "7 11 15"] [-Gamma "5 8 11"] [-Tol 0.05]
#   e.g. the angling follow-up: -Name vfocus2 -Sigma "0.0004" -Angling "15 19 23" -Gamma "8"
param([string]$Name = "vfocus1",
      [string]$Sigma = "0.0004 0.0012 0.0022",
      [string]$Angling = "7 11 15",
      [string]$Gamma = "5 8 11",
      [double]$Tol = 0.05)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$log = Join-Path $deck "Results\vfocus\log_$Name.txt"
New-Item -ItemType Directory -Force (Join-Path $deck "Results\vfocus") | Out-Null
"[$(Get-Date -Format HH:mm:ss)] launching VFocusScan --name $Name --sigma $Sigma --angling $Angling --gamma $Gamma --rank-tol $Tol, scan only (pid $PID)" | Out-File -Append -Encoding ascii $log
$cmd = "python -u VFocusScan.py --name $Name --sigma $Sigma --angling $Angling --gamma $Gamma --rank-tol $Tol"
Invoke-Expression "$cmd 2>&1 | Out-File -Append -Encoding utf8 `"$log`""
"[$(Get-Date -Format HH:mm:ss)] RunVFocus.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
