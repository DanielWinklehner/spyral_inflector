"""Merge BunchScan2 jsons (e.g. the coarse and the fine beam-rotation scans) into one json
with 'best' recomputed, and draw the per-angle best-transmission figure.

The json is written BEFORE the png: RunOvernight.ps1 waits for the png, then reads the json.

    python merge_scans.py OUT_TAG IN1.json IN2.json ...
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

out_tag, ins = sys.argv[1], sys.argv[2:]
out_dir = os.path.dirname(os.path.abspath(ins[0]))
merged, sources = None, []
for f in ins:
    d = json.load(open(f))
    sources.append((d["tag"], d["results"]))
    if merged is None:
        merged = dict(d)
        merged["results"] = list(d["results"])
    else:
        merged["results"] += d["results"]
merged["tag"] = out_tag
merged["sources"] = [t for t, _ in sources]
merged["best"] = max(merged["results"], key=lambda r: r["transmission"])
with open(os.path.join(out_dir, "scan2_{}.json".format(out_tag)), "w") as fh:
    json.dump(merged, fh, indent=2)

fig, ax = plt.subplots(figsize=(9, 5))
for tag, res in sources:
    by = {}
    for r in res:
        if r["phi"] not in by or r["transmission"] > by[r["phi"]]["transmission"]:
            by[r["phi"]] = r
    ph = sorted(by)
    npts = max(sum(1 for r in res if r["phi"] == p) for p in ph)
    ax.plot(ph, [100 * by[p]["transmission"] for p in ph], "o-", label="{} ({} quad points per angle)".format(tag, npts))
b = merged["best"]
ax.set_xlabel("beam rotation phi [deg]")
ax.set_ylabel("best transmission over the quad grid [%]")
ax.set_title("scan '{}': {:,d} particles per point, best {:.1f} % at phi {:.0f}, q1 {:+.0f} / q2 {:+.0f} V".format(
    out_tag, merged["n"], 100 * b["transmission"], b["phi"], b["q1"], b["q2"]), fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "scan2_{}.png".format(out_tag)), dpi=140)
print("merged {} points from {} into scan2_{}.json/.png; best {:.1f} % at phi {:.0f}, q1 {:+.0f} / q2 {:+.0f} V".format(
    len(merged["results"]), " + ".join(t for t, _ in sources), out_tag, 100 * b["transmission"], b["phi"], b["q1"], b["q2"]))
