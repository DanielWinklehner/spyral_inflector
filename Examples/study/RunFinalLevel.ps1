# Self-consistent final push at a space-charge-levelled setting (follow-up 4 of the vfocus campaign):
# the levelled knobs are held fixed (both truncations, dz) and the quads pinned at the SC-retuned
# pair, so the rotation R iterates to its own fixed point and the full runs, hand-offs and Baseline
# exports are made at that R. Waits for the doublet rotation scan's report first (RAM/GPU).
#   powershell -File RunFinalLevel.ps1 [-Name vfocus1_level] [-Trunc "0.34 0.9269377296"] [-Dz 0.0045673337] [-Q1 6575] [-Q2 -6700]
param([string]$Name = "vfocus1_level",
      [string]$Knobs = "Results\vfocus\vfocus1\winner_knobs_refined.json",
      [string]$Trunc = "0.34 0.9269377296074524",
      [double]$Dz = 0.004567333738100502,
      [double]$Q1 = 6575, [double]$Q2 = -6700,
      [string]$WaitFor = "Results\final\vfocus1_final\rot_scan.md")
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$log = Join-Path $deck "Results\vfocus\log_$Name.txt"
"[$(Get-Date -Format HH:mm:ss)] waiting for $deck\$WaitFor (pid $PID)" | Out-File -Append -Encoding ascii $log
while (-not (Test-Path (Join-Path $deck $WaitFor))) { Start-Sleep -Seconds 60 }
"[$(Get-Date -Format HH:mm:ss)] final push $Name at the levelled knobs: truncations $Trunc deg, dz $Dz m, quads $Q1 / $Q2 V" | Out-File -Append -Encoding ascii $log
$cmd = "python -u RunFinalPush.py --name $Name --bfield `"$deck\Fields\final1_baseline_1mm.pickle`" --knobs-json `"$deck\$Knobs`" --fix-truncations $Trunc --fix-dz $Dz --rank vert --rank-tol 0.05 --no-q2-exit-plate --q1 $Q1 $Q1 1 --q2 $Q2 $Q2 1 --fine-step 0"
Invoke-Expression "$cmd 2>&1 | Out-File -Append -Encoding utf8 `"$log`""
"[$(Get-Date -Format HH:mm:ss)] RunFinalLevel.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
