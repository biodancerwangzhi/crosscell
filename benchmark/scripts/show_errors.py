#!/usr/bin/env python3
"""Show actual error messages from benchmark results."""
import json
import re
import sys

tool = sys.argv[1] if len(sys.argv) > 1 else "zellkonverter"
fpath = f"benchmark/results/{tool}.json"

with open(fpath, encoding="utf-8") as f:
    data = json.load(f)

for direction in ["rds_to_h5ad", "h5ad_to_rds"]:
    results = data.get(direction, [])
    for r in results:
        if r.get("status") != "failed":
            continue
        err = r.get("error", "")
        # Extract the actual R error line(s)
        lines = err.split("\n")
        error_lines = [l for l in lines if "Error" in l or "error" in l or "cannot" in l or "TIMEOUT" in l]
        if not error_lines:
            # Take last 3 non-empty lines
            error_lines = [l for l in lines if l.strip()][-3:]
        print(f"{direction} | {r['test_id']}: {' | '.join(error_lines[:3])}")
