"""Print the best point of a BunchScan2 json as 'phi alpha1 alpha2 q1 q2 vscale transmission'.

Best = highest transmission; among points within TIE of it (the binomial noise of 1,500
particles is 1.1 points) the one with the lowest spiral-electrode loss is preferred --
Daniel (2026-09-06): reduce anode/cathode losses rather than grounded-aperture losses when
the transmission is the same.
"""
import json
import sys

TIE = 0.01

d = json.load(open(sys.argv[1]))
res = d["results"]
top = max(r["transmission"] for r in res)
tied = [r for r in res if r["transmission"] >= top - TIE]
b = min(tied, key=lambda r: (r["lost_spiral"], -r["transmission"]))
print("{:.1f} {:.1f} {:.1f} {:.0f} {:.0f} {:.4f} {:.4f}".format(
    b.get("phi", 0.0), b.get("alpha1", 0.0), b.get("alpha2", 0.0), b["q1"], b["q2"], b.get("vscale", 1.0), b["transmission"]))
