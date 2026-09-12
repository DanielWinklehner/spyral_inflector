# Triplet scan on the `triplet` branch (worktree): the package is taken from this worktree
# (PYTHONPATH), the deck from SI_DECK. Runs on any machine with the accel-dev-env environment,
# the deck synced and this branch checked out.
#   powershell -File RunTriplet.ps1 [-Name triplet1] [-Deck "D:\MIT Dropbox\...\Spiral_inflector"] [-Extra "--sigma 0.0004 --angling 15"]
param([string]$Name = "triplet1",
      [string]$Deck = "D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector",
      [string]$Extra = "")
# conda: the first hook found (CONDA_EXE, or anaconda3 / miniconda3 in the usual places)
$hooks = @("$env:USERPROFILE\anaconda3\shell\condabin\conda-hook.ps1", "$env:USERPROFILE\miniconda3\shell\condabin\conda-hook.ps1",
           "$env:LOCALAPPDATA\anaconda3\shell\condabin\conda-hook.ps1", "C:\ProgramData\anaconda3\shell\condabin\conda-hook.ps1")
if ($env:CONDA_EXE) { $hooks = @((Join-Path (Split-Path -Parent (Split-Path -Parent $env:CONDA_EXE)) "shell\condabin\conda-hook.ps1")) + $hooks }
$hook = $hooks | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $hook) { throw "conda-hook.ps1 not found; edit RunTriplet.ps1" }
& $hook
conda activate accel-dev-env
$study = $PSScriptRoot
$repo = Split-Path -Parent (Split-Path -Parent $study)
$env:PYTHONPATH = $repo
$env:SI_DECK = $Deck
Set-Location $study
New-Item -ItemType Directory -Force (Join-Path $Deck "Results\triplet") | Out-Null
$log = Join-Path $Deck "Results\triplet\log_$Name.txt"
"[$(Get-Date -Format HH:mm:ss)] TripletScan --name $Name $Extra (package $repo, pid $PID)" | Out-File -Append -Encoding ascii $log
$cmd = "python -u TripletScan.py --name $Name $Extra"
Invoke-Expression "$cmd >> `"$log`" 2>&1"
"[$(Get-Date -Format HH:mm:ss)] RunTriplet.ps1 finished (exit $LASTEXITCODE)" | Out-File -Append -Encoding ascii $log
