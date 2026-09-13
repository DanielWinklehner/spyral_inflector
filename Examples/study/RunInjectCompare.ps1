# Spatial bunch versus arrival-time injection of the RFQ core, on a finished final-push run.
# Waits for the full DFO-LS optimization to finish (its report full_opt.md) so the two do not
# share the GPU, then runs the injected-convention twins and writes <run>\inject_compare.md.
#   powershell -File RunInjectCompare.ps1 [-Run vfocus1_level] [-WaitFor full_opt.md] [-Cases "vac sc"]
param([string]$Run = "vfocus1_level",
      [string]$WaitFor = "full_opt.md",
      [string]$Cases = "vac sc")
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$rundir = Join-Path $deck "Results\final\$Run"
$log = Join-Path $deck "Results\vfocus\log_${Run}_injectcompare.txt"
if ($WaitFor -and $WaitFor -ne "none") {
    "[$(Get-Date -Format HH:mm:ss)] waiting for $rundir\$WaitFor (pid $PID)" | Out-File -Append -Encoding ascii $log
    while (-not (Test-Path (Join-Path $rundir $WaitFor))) { Start-Sleep -Seconds 120 }
    "[$(Get-Date -Format HH:mm:ss)] $WaitFor present; starting the comparison" | Out-File -Append -Encoding ascii $log
}
$cmd = "python -u InjectCompare.py --run `"$rundir`" --cases $Cases"
Invoke-Expression "$cmd 2>&1 | Out-File -Append -Encoding utf8 `"$log`""
"[$(Get-Date -Format HH:mm:ss)] RunInjectCompare.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
