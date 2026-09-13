# Follow-up 3 of the vfocus campaign: doublet quad-rotation scan by superposition on a final-push run.
#   powershell -File RunDoubletRot.ps1 [-Run vfocus1_final] [-Q "6575 -6700"] [-Alpha1 "-30 -20 -10 0 10 20 30"] [-Alpha2 "..."] [-N 1000]
param([string]$Run = "vfocus1_final",
      [string]$Q = "6575 -6700",
      [string]$Alpha1 = "-30 -20 -10 0 10 20 30",
      [string]$Alpha2 = "-30 -20 -10 0 10 20 30",
      [int]$N = 1000)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$rundir = Join-Path $deck "Results\final\$Run"
$log = Join-Path $deck "Results\vfocus\log_${Run}_rot.txt"
"[$(Get-Date -Format HH:mm:ss)] DoubletRotScan --run $rundir --q $Q --alpha1 $Alpha1 --alpha2 $Alpha2 --n $N (pid $PID)" | Out-File -Append -Encoding ascii $log
$cmd = "python -u DoubletRotScan.py --run `"$rundir`" --q $Q --alpha1 $Alpha1 --alpha2 $Alpha2 --n $N"
Invoke-Expression "$cmd 2>&1 | Out-File -Append -Encoding utf8 `"$log`""
"[$(Get-Date -Format HH:mm:ss)] RunDoubletRot.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
