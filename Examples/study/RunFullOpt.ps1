# Full DFO-LS optimization with 8 mA space charge in the loop (FullOpt.py), seeded at a final-push run.
#   powershell -File RunFullOpt.ps1 [-Run vfocus1_level] [-BestJson basis\scan2_coarse.json] [-MaxFun 80] [-N 8000]
param([string]$Run = "vfocus1_level",
      [string]$BestJson = "basis\scan2_coarse.json",
      [int]$MaxFun = 80,
      [int]$N = 8000,
      [int]$WaitForPid = 0)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$rundir = Join-Path $deck "Results\final\$Run"
$log = Join-Path $deck "Results\vfocus\log_${Run}_fullopt.txt"
if ($WaitForPid -gt 0) {
    # share the GPU with nobody: wait for another job (e.g. a central-region simulation) to finish first
    "[$(Get-Date -Format HH:mm:ss)] waiting for process $WaitForPid to exit before starting (pid $PID)" | Out-File -Append -Encoding ascii $log
    while (Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue) { Start-Sleep -Seconds 120 }
    "[$(Get-Date -Format HH:mm:ss)] process $WaitForPid is gone; starting" | Out-File -Append -Encoding ascii $log
}
"[$(Get-Date -Format HH:mm:ss)] FullOpt --run $rundir --best-json $BestJson --maxfun $MaxFun --n $N (pid $PID)" | Out-File -Append -Encoding ascii $log
python -u FullOpt.py --run $rundir --best-json $BestJson --maxfun $MaxFun --n $N 2>&1 | Out-File -Append -Encoding utf8 $log
"[$(Get-Date -Format HH:mm:ss)] RunFullOpt.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
