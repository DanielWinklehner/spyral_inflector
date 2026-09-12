# The laptop's night: triplet scan on the final1 shape, space-charge levelling of its best
# setting with the full runs, then (when the desktop's refined vfocus1 winner has synced) the
# same on the winner's shape. Resumable: finished stages are reused.
#   powershell -File RunTripletNight.ps1 -Deck "<deck path>" [-WaitHours 6]
param([string]$Deck = "D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector",
      [double]$WaitHours = 6)
$hooks = @("$env:USERPROFILE\anaconda3\shell\condabin\conda-hook.ps1", "$env:USERPROFILE\miniconda3\shell\condabin\conda-hook.ps1",
           "$env:LOCALAPPDATA\anaconda3\shell\condabin\conda-hook.ps1", "C:\ProgramData\anaconda3\shell\condabin\conda-hook.ps1")
if ($env:CONDA_EXE) { $hooks = @((Join-Path (Split-Path -Parent (Split-Path -Parent $env:CONDA_EXE)) "shell\condabin\conda-hook.ps1")) + $hooks }
$hook = $hooks | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $hook) { throw "conda-hook.ps1 not found; edit RunTripletNight.ps1" }
& $hook
conda activate accel-dev-env
$study = $PSScriptRoot
$repo = Split-Path -Parent (Split-Path -Parent $study)
$env:PYTHONPATH = $repo
$env:SI_DECK = $Deck
Set-Location $study
New-Item -ItemType Directory -Force (Join-Path $Deck "Results\triplet") | Out-Null
$log = Join-Path $Deck "Results\triplet\log_night.txt"
function Say($m) { "[$(Get-Date -Format HH:mm:ss)] $m" | Out-File -Append -Encoding utf8 $log; Write-Host $m }

Say "night chain started (package $repo, pid $PID)"
Say "stage 1: triplet1 scan (final1 shape)"
python -u TripletScan.py --name triplet1 2>&1 | Out-File -Append -Encoding utf8 $log
Say "stage 2: space-charge levelling of triplet1 + full runs"
python -u SCLevel.py --run (Join-Path $Deck "Results\triplet\triplet1") --final 2>&1 | Out-File -Append -Encoding utf8 $log

$winner = Join-Path $Deck "Results\vfocus\vfocus1\winner_knobs_refined.json"
$deadline = (Get-Date).AddHours($WaitHours)
while (-not (Test-Path $winner) -and (Get-Date) -lt $deadline) { Start-Sleep -Seconds 300 }
if (Test-Path $winner) {
    $k = Get-Content $winner -Raw | ConvertFrom-Json
    Say "stage 3: triplet2 scan on the vfocus1 winner shape sigma $($k.sigma) angling $($k.angling) gamma $($k.gamma)"
    python -u TripletScan.py --name triplet2 --sigma $k.sigma --angling $k.angling --gamma $k.gamma 2>&1 | Out-File -Append -Encoding utf8 $log
    Say "stage 4: space-charge levelling of triplet2 + full runs"
    python -u SCLevel.py --run (Join-Path $Deck "Results\triplet\triplet2") --final 2>&1 | Out-File -Append -Encoding utf8 $log
} else {
    Say "no refined vfocus1 winner within $WaitHours h; stages 3-4 skipped"
}
Say "night chain finished"
