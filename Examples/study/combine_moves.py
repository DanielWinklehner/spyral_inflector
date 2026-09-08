"""Pick the single-knob moves that raised the transmission by more than 1.5 points without
raising the losses on the spiral or quad electrodes, combine up to three of them (one per
knob) and run that geometry as the point 'combo'.
"""
import json
import os
import subprocess
import sys

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
W = os.path.join(DECK, "Results", "wiggle")
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
KNOB_OF = {"sigma_m": ("--sigma", "0.0012"), "sigma_p": ("--sigma", "0.0032"),
           "angling_m": ("--angling", "7.0"), "angling_p": ("--angling", "15.0"),
           "gamma_m": ("--gamma", "2.0"), "gamma_p": ("--gamma", "8.0"),
           "aspect_m": ("--aspect", "2.0"), "aspect_p": ("--aspect", "2.8"),
           "gap_p": ("--gap", "0.022"), "bore_p": ("--quad-bore", "0.016", "--aper-hole", "0.0155"),
           "quad2_up": ("--quad-z2", "-0.21"), "quad2_down": ("--quad-z2", "-0.17"),
           "hole_m": ("--aper-hole", "0.0105")}

base = json.load(open(os.path.join(W, "base", "summary.json")))["bunch"]
gains = []
for name, argv in KNOB_OF.items():
    fn = os.path.join(W, name, "summary.json")
    if not os.path.exists(fn):
        continue
    s = json.load(open(fn))
    if "bunch" not in s:
        continue
    b = s["bunch"]
    gain = 100 * (b["transmission"] - base["transmission"])
    arc = 100 * ((b["lost_spiral"] + b["lost_quads"]) - (base["lost_spiral"] + base["lost_quads"]))
    print("  {:12s} transmission {:+.1f} points, arc-prone losses {:+.1f} points".format(name, gain, arc))
    if gain > 1.5 and arc <= 0.5:
        gains.append((gain, name, argv))
gains.sort(reverse=True)
chosen, knobs_used, extra = [], set(), []
for gain, name, argv in gains:
    knob = argv[0]
    if knob in knobs_used:
        continue
    knobs_used.add(knob)
    chosen.append(name)
    extra += list(argv)
    if len(chosen) == 3:
        break
if not chosen:
    print("no move improved the transmission by more than 1.5 points without adding electrode losses; no combination run")
    sys.exit(0)
print("combining: {}".format(chosen))
st = json.load(open(os.path.join(W, "stage1_settings.json"), encoding="utf-8-sig"))   # written by PowerShell, with BOM
cmd = [sys.executable, "-u", os.path.join(SCRIPTS, "GeometryPoint.py"), "--name", "combo",
       "--phi", str(st["phi"]), "--rotate-quads", str(st["alpha1"]), str(st["alpha2"]),
       "--q1", str(st["q1"] - 1500), str(st["q1"] + 1500), "3", "--q2", str(st["q2"] - 1500), str(st["q2"] + 1500), "3"] + extra
with open(os.path.join(W, "log_combo.txt"), "w", encoding="utf-8") as fh:
    fh.write("moves: {}\n".format(chosen))
    p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=W)
print("combo exit code {}".format(p.returncode))
