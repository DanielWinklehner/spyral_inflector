# Night chain v2 (2026-09-06, after the geometry candidates): Phase B re-optimizes the chosen
# geometry at the beam's real mean energy with 2.5 mm fields and adds the skew basis fields;
# Phase C tunes the cheap knobs (beam angle, quad rotations, 375 V quad grid, spiral voltage),
# checks with a 10k bunch and a direct BEM solve; Phase D shifts the inflector by minus the
# bunch's mean exit height and re-checks; Phase E tracks all core particles with and without
# PyAMG space charge on the final direct field; then the report.
#
# ---- settings (filled in from the candidate results) ----
$cand   = "bore18_slot19_pg5L"                      # geometry name (for the report)
$knobs  = @("--quad-bore", "0.018", "--aper-hole", "0.0175", "--slot-width", "0.019", "--exit-opening", "0.023", "0.040", "--gap", "0.019", "--plate-gap", "0.005", "--quad-z1", "-0.264", "--quad-len", "0.055", "--quad-len2", "0.060", "--shared-plates")                    # e.g. "--quad-bore","0.016","--aper-hole","0.0155","--slot-width","0.019"
$orient = @()                   # @("--swap-xy") or @()
$phiC   = 90                        # centre beam angle [deg] on top of the orientation
$qc1    = 6150                         # quad-1 voltage centre [V]
$qc2    = -7400                         # quad-2 voltage centre [V]
$qw     = 1500                          # Phase B inner grid half-width [V]
$qn     = 4                          # Phase B inner grid points per axis
$nAll   = 43969                           # all core particles
$imA    = 8.0                             # design current on the core distribution [mA] (Daniel, 2026-09-06)
# ---------------------------------------------------------
& "C:\Users\Daniel\anaconda3\shell\condabin\conda-hook.ps1"
conda activate accel-dev-env

$deck = "D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector"
$scripts = "$deck\Scripts"
$p = "$deck\Particles\MIT RFQ Beamdynamics\ext3_exit_core_as_txt.txt"
$b = "$deck\Fields\HCHC-60_CentralBField_z-40to5cm_1mm.pickle"
$fin = "$deck\Results\final_pg5L"
New-Item -ItemType Directory -Force $fin | Out-Null
Set-Location $fin

function Stamp($msg) { Write-Output ("[{0}] {1}" -f (Get-Date -Format "HH:mm:ss"), $msg) }
function Best($json) { return ((python "$scripts\scan_best.py" $json) -split " ") }   # phi a1 a2 q1 q2 vscale T
function JsonValue($json, $expr) { return (python -c "import json; d = json.load(open(r'$json')); print($expr)") }

# ---------------------------------------------------------------- Phase B: final geometry at 2.5 mm
Stamp "Phase B: geometry '$cand' re-optimized at the beam energy, 2.5 mm fields"
$g1 = "$deck\Results\wiggle\pg5La"; $s1 = "$deck\Geometry\wiggle\pg5La_steps"
if (Test-Path "$g1\bunch_bunch.npz") {
    Stamp "Phase B: pg5La already complete, skipping the geometry run"
} else {
    python -u "$scripts\GeometryPoint.py" --name pg5La @knobs @orient --phi $phiC --rotate-quads 0 0 --res 0.0025 `
        --q1 ($qc1 - $qw) ($qc1 + $qw) $qn --q2 ($qc2 - $qw) ($qc2 + $qw) $qn --n-scan 1500 --n-bunch 5000 2>&1 |
        Out-File -Encoding utf8 "$fin\log_pg5La.txt"
    Get-Content "$fin\log_pg5La.txt" | Select-String -Pattern "design energy|optimizer:|\] done:|failed" | ForEach-Object { Stamp $_.Line }
}
$common = @("--step-dir", $s1, "--voltages", "$g1\voltages.csv", "--state", "$g1\state.pickle", "--out-dir", $g1,
            "--res", "0.0025", "--bfield", $b, "--no-test", "--spiral-voltage", "0")
if (-not (Test-Path "$g1\ef_itp_q1skew.pickle")) {
    python -u "$scripts\TrackFromStep.py" --tag q1skew --quad-voltages 3500 0 --rotate-quads 45 0 @common 2>&1 | Out-File -Encoding utf8 "$fin\log_q1skew.txt"
}
if (-not (Test-Path "$g1\ef_itp_q2skew.pickle")) {
    python -u "$scripts\TrackFromStep.py" --tag q2skew --quad-voltages 0 3500 --rotate-quads 0 45 @common 2>&1 | Out-File -Encoding utf8 "$fin\log_q2skew.txt"
}
Stamp "skew basis fields solved"
$r = Best "$g1\scan2_inner.json"; $q1 = [double]$r[3]; $q2 = [double]$r[4]
Stamp "pg5La inner scan best: q1=$q1 q2=$q2 transmission=$($r[6])"
$V = [double](JsonValue "$g1\summary.json" "d['optimizer']['voltage_V']")
$dz1 = [double](JsonValue "$g1\summary.json" "d['optimizer']['dz_mm']")
Stamp "pg5La spiral voltage $V V, dz $dz1 mm"

# ---------------------------------------------------------------- Phase C: cheap knobs on the 2.5 mm fields
$scan = @("--reload-dir", $g1, "--out-dir", $g1, "--step-dir", $s1, "--basis", "spiral", "q1", "q1skew", "q2", "q2skew",
          "--particles", $p, "--bfield", $b, "--n", "1500") + $orient
Stamp "Phase C1: beam angle $($phiC - 10)..$($phiC + 10) deg in 5 deg steps, 4x4 quad grid at 667 V"
python "$scripts\scan_complete.py" "$g1\scan2_c1_angle.json" 80 | Out-Null
if (Test-Path "$g1\c1_best_override.txt") { Stamp "C1: result taken from c1_best_override.txt (the 22:34 scan; its json was overwritten by a relaunch)" }
elseif ($LASTEXITCODE -eq 0) { Stamp "C1: finished scan reused" } else {
python -u "$scripts\BunchScan2.py" --tag c1_angle @scan --phi ($phiC - 10) ($phiC - 5) $phiC ($phiC + 5) ($phiC + 10) `
    --q1 ($q1 - 1000) ($q1 + 1000) 4 --q2 ($q2 - 1000) ($q2 + 1000) 4 2>&1 | Out-File -Encoding utf8 "$fin\log_c1_angle.txt"
}
if (Test-Path "$g1\c1_best_override.txt") {
    $r = ((Get-Content "$g1\c1_best_override.txt" | Select-Object -First 1) -split " ")
    if ($r[3] -eq "B") { $r[3] = "$q1"; $r[4] = "$q2"; $r[6] = "n/a" }   # "B" = keep the Phase B quad voltages
} else { $r = Best "$g1\scan2_c1_angle.json" }
$phi = [double]$r[0]; $q1 = [double]$r[3]; $q2 = [double]$r[4]
Stamp "C1 best: phi=$phi q1=$q1 q2=$q2 transmission=$($r[6])"

Stamp "Phase C2: quad rotations 0 8 16 x 8 16 24 deg, 3x3 at 750 V"
python "$scripts\scan_complete.py" "$g1\scan2_c2_alpha.json" 81 | Out-Null
if ($LASTEXITCODE -eq 0) {
    Stamp "C2: finished scan reused"
} else {
    python -u "$scripts\BunchScan2.py" --tag c2_alpha @scan --phi $phi --alpha1 0 8 16 --alpha2 8 16 24 `
        --q1 ($q1 - 750) ($q1 + 750) 3 --q2 ($q2 - 750) ($q2 + 750) 3 2>&1 | Out-File -Encoding utf8 "$fin\log_c2_alpha.txt"
}
$r = Best "$g1\scan2_c2_alpha.json"; $a1 = [double]$r[1]; $a2 = [double]$r[2]; $q1 = [double]$r[3]; $q2 = [double]$r[4]
Stamp "C2 best: alpha=$a1/$a2 q1=$q1 q2=$q2 transmission=$($r[6])"

if ($true) {
    Copy-Item "$g1\scan2_c2_alpha.json" "$g1\scan2_c2_alpha_all.json" -Force
    Stamp "C2x skipped (rotation range already brackets the ridge)"
} else {
# The C2 optimum (8/16) sat at the edge of the range on a ridge alpha2 - alpha1 ~ 8 deg (the Larmor
# rotation between the quads); extend along and across the ridge, then merge with C2.
Stamp "Phase C2x: quad rotations 8/16/24 x 16/24/32 deg, 3x3 at 750 V around the C2 best"
python "$scripts\scan_complete.py" "$g1\scan2_c2x_alpha_ext.json" 81 | Out-Null
if ($LASTEXITCODE -eq 0) { Stamp "C2x: finished scan reused" } else {
python -u "$scripts\BunchScan2.py" --tag c2x_alpha_ext @scan --phi $phi --alpha1 8 16 24 --alpha2 16 24 32 `
    --q1 ($q1 - 750) ($q1 + 750) 3 --q2 ($q2 - 750) ($q2 + 750) 3 2>&1 | Out-File -Encoding utf8 "$fin\log_c2x_alpha_ext.txt"
}
python "$scripts\merge_scans.py" c2_alpha_all "$g1\scan2_c2_alpha.json" "$g1\scan2_c2x_alpha_ext.json" | ForEach-Object { Stamp $_ }
$r = Best "$g1\scan2_c2_alpha_all.json"; $a1 = [double]$r[1]; $a2 = [double]$r[2]; $q1 = [double]$r[3]; $q2 = [double]$r[4]
Stamp "C2x best (merged): alpha=$a1/$a2 q1=$q1 q2=$q2 transmission=$($r[6])"
}

# an override re-centres the fine grid (e.g. when the C3 optimum sat on the edge of its grid):
# 'phi alpha1 alpha2 q1 q2 ...'
if (Test-Path "$g1\c2_best_override.txt") {
    $r = ((Get-Content "$g1\c2_best_override.txt" | Select-Object -First 1) -split " ")
    $a1 = [double]$r[1]; $a2 = [double]$r[2]; $q1 = [double]$r[3]; $q2 = [double]$r[4]
    Stamp "C2 result overridden by c2_best_override.txt: alpha=$a1/$a2 q1=$q1 q2=$q2 (fine grid re-centred)"
}
Stamp "Phase C3: 5x5 quad grid at 375 V"
python "$scripts\scan_complete.py" "$g1\scan2_c3_fine.json" 25 | Out-Null
if ($LASTEXITCODE -eq 0) { Stamp "C3: finished scan reused" } else {
python -u "$scripts\BunchScan2.py" --tag c3_fine @scan --phi $phi --alpha1 $a1 --alpha2 $a2 `
    --q1 ($q1 - 750) ($q1 + 750) 5 --q2 ($q2 - 750) ($q2 + 750) 5 2>&1 | Out-File -Encoding utf8 "$fin\log_c3_fine.txt"
}
$fit = @(python "$scripts\scan_fit.py" "$g1\scan2_c3_fine.json"); $r = $fit[0] -split " "; $q1 = [double]$r[3]; $q2 = [double]$r[4]
Stamp "C3 best (quadratic surrogate of the 5x5 grid): q1=$q1 q2=$q2 transmission=$($r[6]); $($fit[1])"

Stamp "Phase C4: spiral voltage scale"
python "$scripts\scan_complete.py" "$g1\scan2_c4_vscale.json" 5 | Out-Null
if ($LASTEXITCODE -eq 0) { Stamp "C4: finished scan reused" } else {
python -u "$scripts\BunchScan2.py" --tag c4_vscale @scan --phi $phi --alpha1 $a1 --alpha2 $a2 --points $q1 $q2 `
    --vscale 0.97 0.985 1.0 1.015 1.03 2>&1 | Out-File -Encoding utf8 "$fin\log_c4_vscale.txt"
}
$r = Best "$g1\scan2_c4_vscale.json"; $vs = [double]$r[5]
Stamp "C4 best: vscale=$vs transmission=$($r[6])"
# The voltage trades < 1 point of transmission against the asymptotic vertical angle of the centroid
# (1.03: -1.4 deg, 1.00: +0.2 deg); an override file fixes the scale used from here on.
$vsOverride = $null
if (Test-Path "$g1\c4_vscale_override.txt") {
    $vsOverride = [double](Get-Content "$g1\c4_vscale_override.txt" | Select-Object -First 1)
    Stamp "C4: spiral voltage scale overridden to $vsOverride (level beam; the C4 optimum $vs tilts the centroid)"
    $vs = $vsOverride
}

if ($vsOverride -ne $null) {
    Stamp "Phase C5 skipped (voltage overridden; the direct solve of C6 is the check)"
} else {
    Stamp "Phase C5: 10k bunch on the superposed field"
    python -u "$scripts\BunchTrack.py" --tag c4_vscale_best --reload-dir $g1 --step-dir $s1 --particles $p --bfield $b @orient --phi $phi `
        --n 10000 --out-tag c5_bunch 2>&1 | Out-File -Encoding utf8 "$fin\log_c5_bunch.txt"
    python "$scripts\exit_metrics.py" "$g1\bunch_c5_bunch.npz" "$g1\si_state_c4_vscale_best.pickle" "$fin\exit_metrics_c5.json" | ForEach-Object { Stamp "C5 $_" }
}

$spiralV = $V * $vs
Stamp "Phase C6: direct BEM solve, quads $q1/$q2 V rotated $a1/$a2 deg, spiral $spiralV V"
if (Test-Path "$fin\exit_metrics_c6.json") { Stamp "C6: existing direct-solve result reused" } else {
python -u "$scripts\TrackFromStep.py" --tag c6_direct --step-dir $s1 --voltages "$g1\voltages.csv" --state "$g1\state.pickle" --out-dir $g1 `
    --res 0.0025 --bfield $b --quad-voltages $q1 $q2 --rotate-quads $a1 $a2 --spiral-voltage $spiralV 2>&1 | Out-File -Encoding utf8 "$fin\log_c6_direct.txt"
python -u "$scripts\BunchTrack.py" --tag c6_direct --reload-dir $g1 --step-dir $s1 --particles $p --bfield $b @orient --phi $phi `
    --n 10000 --out-tag c6_bunch 2>&1 | Out-File -Encoding utf8 "$fin\log_c6_bunch.txt"
python "$scripts\exit_metrics.py" "$g1\bunch_c6_bunch.npz" "$g1\si_state_c6_direct.pickle" "$fin\exit_metrics_c6.json" | ForEach-Object { Stamp "C6 $_" }
}
@{cand = $cand; phi = $phi; swap_xy = [bool]($orient -contains "--swap-xy"); alpha1 = $a1; alpha2 = $a2; q1 = $q1; q2 = $q2; vscale = $vs; spiral_V = $spiralV; dz1_mm = $dz1} |
    ConvertTo-Json | Out-File -Encoding utf8 "$fin\stageC_settings.json"

# 5 mm fields lost 18 points against 2.5 mm on the same geometry (19 mm gap, 4 cells across), so
# check 2.5 mm against 1.25 mm on the direct solve; the finer one is used for the final fields
# if the 10k transmissions differ by more than 2 points.
Stamp "Phase C7: resolution check, direct solve at 1.25 mm"
if (Test-Path "$fin\exit_metrics_c7.json") { Stamp "C7: existing 1.25 mm result reused" } else {
python -u "$scripts\TrackFromStep.py" --tag c7_direct_hr --step-dir $s1 --voltages "$g1\voltages.csv" --state "$g1\state.pickle" --out-dir $g1 `
    --res 0.00125 --bfield $b --quad-voltages $q1 $q2 --rotate-quads $a1 $a2 --spiral-voltage $spiralV 2>&1 | Out-File -Encoding utf8 "$fin\log_c7_direct_hr.txt"
python -u "$scripts\BunchTrack.py" --tag c7_direct_hr --reload-dir $g1 --step-dir $s1 --particles $p --bfield $b @orient --phi $phi `
    --n 10000 --out-tag c7_bunch 2>&1 | Out-File -Encoding utf8 "$fin\log_c7_bunch.txt"
python "$scripts\exit_metrics.py" "$g1\bunch_c7_bunch.npz" "$g1\si_state_c7_direct_hr.pickle" "$fin\exit_metrics_c7.json" | ForEach-Object { Stamp "C7 $_" }
}
$t25 = [double](JsonValue "$fin\exit_metrics_c6.json" "d['transmission']")
$t125 = [double](JsonValue "$fin\exit_metrics_c7.json" "d['transmission']")
$resFinal = 0.0025
if ([math]::Abs($t125 - $t25) -gt 0.02) { $resFinal = 0.00125 }
Stamp "resolution check: 2.5 mm $($t25 * 100) %, 1.25 mm $($t125 * 100) % -> final fields at $($resFinal * 1e3) mm"

# ---------------------------------------------------------------- Phase D: centroid on the median plane
$zc = [double](JsonValue "$fin\exit_metrics_c6.json" "d.get('z_asym_mean_mm', d['z_exit_mean_mm'])")
$dz2 = ($dz1 - $zc) * 1e-3
Stamp "Phase D: bunch mean asymptotic exit height (past the fringe) $zc mm -> inflector shift dz = $($dz2 * 1e3) mm (was $dz1 mm); re-export, re-solve, re-tune"
$g2 = "$deck\Results\wiggle\pg5Lb"; $s2 = "$deck\Geometry\wiggle\pg5Lb_steps"
if (Test-Path "$g2\bunch_bunch.npz") { Stamp "Phase D: pg5Lb geometry already complete, reused" } else {
python -u "$scripts\GeometryPoint.py" --name pg5Lb @knobs @orient --phi $phi --rotate-quads $a1 $a2 --fix-dz $dz2 --res 0.0025 `
    --q1 ($q1 - 375) ($q1 + 375) 3 --q2 ($q2 - 375) ($q2 + 375) 3 --n-scan 1500 --n-bunch 5000 2>&1 | Out-File -Encoding utf8 "$fin\log_pg5Lb.txt"
Get-Content "$fin\log_pg5Lb.txt" | Select-String -Pattern "optimizer:|\] done:|failed" | ForEach-Object { Stamp $_.Line }
}
$r = Best "$g2\scan2_inner.json"; $q1 = [double]$r[3]; $q2 = [double]$r[4]
$V2 = [double](JsonValue "$g2\summary.json" "d['optimizer']['voltage_V']")
$spiralV2 = $V2 * $vs
Stamp "pg5Lb: quads $q1/$q2 V, spiral $V2 V x $vs = $spiralV2 V; direct solve + 10k bunch"
python -u "$scripts\TrackFromStep.py" --tag d_direct --step-dir $s2 --voltages "$g2\voltages.csv" --state "$g2\state.pickle" --out-dir $g2 `
    --res $resFinal --bfield $b --quad-voltages $q1 $q2 --rotate-quads $a1 $a2 --spiral-voltage $spiralV2 2>&1 | Out-File -Encoding utf8 "$fin\log_d_direct.txt"
python -u "$scripts\BunchTrack.py" --tag d_direct --reload-dir $g2 --step-dir $s2 --particles $p --bfield $b @orient --phi $phi `
    --n 10000 --out-tag d_bunch 2>&1 | Out-File -Encoding utf8 "$fin\log_d_bunch.txt"
python "$scripts\exit_metrics.py" "$g2\bunch_d_bunch.npz" "$g2\si_state_d_direct.pickle" "$fin\exit_metrics_d.json" | ForEach-Object { Stamp "D $_" }
@{cand = $cand; phi = $phi; swap_xy = [bool]($orient -contains "--swap-xy"); alpha1 = $a1; alpha2 = $a2; q1 = $q1; q2 = $q2; vscale = $vs; spiral_V = $spiralV2; dz2_mm = $dz2 * 1e3; steps = $s2; res_m = $resFinal} |
    ConvertTo-Json | Out-File -Encoding utf8 "$fin\final_settings.json"

# ---------------------------------------------------------------- Phase E: all core particles, with and without space charge
Stamp "Phase E: $nAll particles on the final direct field, no space charge and PyAMG space charge at $imA mA (parallel)"
$py = "C:\Users\Daniel\anaconda3\envs\accel-dev-env\python.exe"
$q = [char]34
$argsNo = @("-u", "$q$scripts\BunchTrack.py$q", "--tag", "d_direct", "--reload-dir", "$q$g2$q", "--step-dir", "$q$s2$q",
            "--particles", "$q$p$q", "--bfield", "$q$b$q") + $orient + @("--phi", "$phi", "--n", "$nAll", "--out-tag", "e_nosc_all")
$argsSc = @("-u", "$q$scripts\BunchTrackSC.py$q", "--tag", "d_direct", "--reload-dir", "$q$g2$q", "--step-dir", "$q$s2$q", "--out-dir", "$q$fin$q",
            "--particles", "$q$p$q", "--bfield", "$q$b$q") + $orient + @("--phi", "$phi", "--n", "$nAll", "--current-ma", "$imA", "--out-tag", "e_sc_all")
$pNo = Start-Process -FilePath $py -ArgumentList $argsNo -RedirectStandardOutput "$fin\log_e_nosc.txt" -RedirectStandardError "$fin\err_e_nosc.txt" -WindowStyle Hidden -PassThru
$pSc = Start-Process -FilePath $py -ArgumentList $argsSc -RedirectStandardOutput "$fin\log_e_sc.txt" -RedirectStandardError "$fin\err_e_sc.txt" -WindowStyle Hidden -PassThru
$pNo | Wait-Process
Stamp "no-SC run finished (exit $($pNo.ExitCode))"
python "$scripts\exit_metrics.py" "$g2\bunch_e_nosc_all.npz" "$g2\si_state_d_direct.pickle" "$fin\exit_metrics_e_nosc.json" | ForEach-Object { Stamp "E no-SC $_" }
$pSc | Wait-Process
Stamp "SC run finished (exit $($pSc.ExitCode))"
Get-Content "$fin\log_e_sc.txt" | Select-String -Pattern "TRANSMITTED|Traceback|Error" | Select-Object -Last 2 | ForEach-Object { Stamp "E SC $($_.Line)" }
python "$scripts\exit_metrics.py" "$fin\bunch_e_sc_all.npz" "$g2\si_state_d_direct.pickle" "$fin\exit_metrics_e_sc.json" | ForEach-Object { Stamp "E SC $_" }

Stamp "report"
python "$scripts\make_final_report.py" --fin "final_pg5L" --run1 "pg5La" --run2 "pg5Lb" --title "longer quads 55/60 mm, shared plate, plate gap 5 mm, SI gap 19 mm, housing exit opening 23 mm (2026-09-07)" 2>&1 | ForEach-Object { Stamp $_ }
Stamp "done"
