# Follow-up 2 of the vfocus campaign: where the transverse emittance grows along the line.
# Four bunch runs with 6-D snapshots every N steps (vfocus1_final at its levelled + SC-retuned
# setting, and final1 at its fine_best setting; vacuum and 8 mA; core particles only), then
# EmittanceAlongLine.py on the four.
#   powershell -File RunEmittance.ps1 [-N 5000] [-Every 10]
param([int]$N = 5000, [int]$Every = 10,
      [string]$Run = "vfocus1_final", [string]$LevelDir = "sc_level\t0.927",
      [string]$Superpose = "6575 0 -6700 0", [double]$ShiftZ = 0.0014123557415166178, [double]$Phi = 2.30894689391354,
      [string]$Ref = "final1", [string]$RefTag = "fine_best", [double]$RefPhi = 2.987565927210454)
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env
$deck = Split-Path -Parent $PSScriptRoot
Set-Location $PSScriptRoot
$log = Join-Path $deck "Results\vfocus\log_emittance.txt"
$rundir = Join-Path $deck "Results\final\$Run"
$refdir = Join-Path $deck "Results\final\$Ref"
$lvl = Join-Path $rundir $LevelDir
$out = Join-Path $rundir "emittance"
New-Item -ItemType Directory -Force $out | Out-Null
$particles = Join-Path $deck "Particles\MIT RFQ Beamdynamics\ext3_exit.dst"
$bfield = Join-Path $deck "Fields\final1_baseline_1mm.pickle"
$sc = "--sc --h 0.004 --resolve-every 16 --current-ma 8.0 --sc-min-particles 200"
$common = "--particles `"$particles`" --bfield `"$bfield`" --out-dir `"$out`" --n $N --snapshot-every $Every --handoff-frame machine --handoff-distance 0.0 --no-plot"
function Say($m) { "[$(Get-Date -Format HH:mm:ss)] $m" | Out-File -Append -Encoding ascii $log }
Say "EMITTANCE runs: $Run (levelled, superpose $Superpose, shift $ShiftZ m) and $Ref ($RefTag); N $N, snapshots every $Every steps (pid $PID)"
$runs = @(
  @("vf_vac", "python -m spyral_inflector.tracking.bunch --tag spiral --reload-dir `"$lvl\basis`" --step-dir `"$lvl\steps`" --basis-dir `"$lvl\basis`" --superpose $Superpose --shift-z $ShiftZ --phi $Phi --out-tag vf_vac $common"),
  @("vf_sc",  "python -m spyral_inflector.tracking.bunch --tag spiral --reload-dir `"$lvl\basis`" --step-dir `"$lvl\steps`" --basis-dir `"$lvl\basis`" --superpose $Superpose --shift-z $ShiftZ --phi $Phi --out-tag vf_sc $common $sc"),
  @("f1_vac", "python -m spyral_inflector.tracking.bunch --tag $RefTag --reload-dir `"$refdir\basis`" --step-dir `"$refdir\steps`" --phi $RefPhi --out-tag f1_vac $common"),
  @("f1_sc",  "python -m spyral_inflector.tracking.bunch --tag $RefTag --reload-dir `"$refdir\basis`" --step-dir `"$refdir\steps`" --phi $RefPhi --out-tag f1_sc $common $sc"))
foreach ($r in $runs) {
  $tag = $r[0]; $cmd = $r[1]
  if (Test-Path (Join-Path $out "bunch_$tag.npz")) { Say "  $tag exists"; continue }
  Say "  run $tag"
  Invoke-Expression "$cmd 2>&1 | Out-File -Encoding utf8 `"$out\log_$tag.txt`""
  Say "  $tag done (exit $LASTEXITCODE)"
}
Say "analysis"
$npz = "`"$out\bunch_vf_vac.npz`" `"$out\bunch_vf_sc.npz`" `"$out\bunch_f1_vac.npz`" `"$out\bunch_f1_sc.npz`""
$rj = "`"$lvl\basis\reload_spiral.json`" `"$lvl\basis\reload_spiral.json`" `"$refdir\basis\reload_spiral.json`" `"$refdir\basis\reload_spiral.json`""
Invoke-Expression "python -u EmittanceAlongLine.py --npz $npz --labels `"$Run vac`" `"$Run 8 mA`" `"$Ref vac`" `"$Ref 8 mA`" --reload-json $rj --out `"$out\emittance_vs_$Ref`" 2>&1 | Out-File -Append -Encoding utf8 `"$log`""
Say "RunEmittance.ps1 finished (exit $LASTEXITCODE); report $out\emittance_vs_$Ref.md"
