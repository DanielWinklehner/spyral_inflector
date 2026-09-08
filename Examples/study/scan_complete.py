"""Exit 0 if a BunchScan2 json exists and holds at least N results (a finished scan), else 1.

    python scan_complete.py scan2_x.json 80
"""
import json
import os
import sys

path, n = sys.argv[1], int(sys.argv[2])
if not os.path.exists(path):
    sys.exit(1)
try:
    k = len(json.load(open(path))["results"])
except Exception:  # noqa: BLE001
    sys.exit(1)
print("{}: {} of {} points".format(os.path.basename(path), k, n))
sys.exit(0 if k >= n else 1)
