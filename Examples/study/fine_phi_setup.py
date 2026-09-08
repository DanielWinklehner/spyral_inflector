"""Decide the fine beam-rotation scan from the coarse scan2_phi_newB.json.

  --check   exit 0 once the last wanted coarse angle (default 90 deg) has all its quad
            points, 1 while it is incomplete, 2 if the json is unreadable (mid-write)
  default   print 'phi q1 q2 transmission' of the best point among the COMPLETED angles
"""
import argparse
import collections
import json
import sys

p = argparse.ArgumentParser()
p.add_argument("json")
p.add_argument("--check", action="store_true")
p.add_argument("--last-phi", type=float, default=90.0)
p.add_argument("--points-per-phi", type=int, default=25)
a = p.parse_args()

try:
    d = json.load(open(a.json))
except Exception as e:  # noqa: BLE001
    print("unreadable: {}".format(e))
    sys.exit(2)
per = collections.defaultdict(list)
for r in d["results"]:
    per[r["phi"]].append(r)

if a.check:
    n = len(per.get(a.last_phi, []))
    print("phi {:.0f}: {} of {} points".format(a.last_phi, n, a.points_per_phi))
    sys.exit(0 if n >= a.points_per_phi else 1)

done = [r for lst in per.values() if len(lst) >= a.points_per_phi for r in lst]
best = max(done, key=lambda r: r["transmission"])
print("{:.1f} {:.0f} {:.0f} {:.4f}".format(best["phi"], best["q1"], best["q2"], best["transmission"]))
