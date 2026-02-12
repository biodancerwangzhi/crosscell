#!/usr/bin/env python3
"""Quick check of benchmark results."""
import json
import os
import sys

results_dir = sys.argv[1] if len(sys.argv) > 1 else "benchmark/results"

for fname in sorted(os.listdir(results_dir)):
    if not fname.endswith(".json"):
        continue
    tool = fname.replace(".json", "")
    if tool in ("dataset_summary", "tool_versions", "validation"):
        continue
    fpath = os.path.join(results_dir, fname)
    try:
        with open(fpath, encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        print(f"=== {tool} === ERROR reading: {e}")
        continue

    print(f"\n=== {tool} ===")
    for direction in ["rds_to_h5ad", "h5ad_to_rds"]:
        results = data.get(direction, [])
        if not results:
            print(f"  {direction}: no data")
            continue
        passed = sum(1 for r in results if r.get("status") == "success")
        failed = sum(1 for r in results if r.get("status") == "failed")
        skipped = sum(1 for r in results if r.get("status") == "skipped")
        total = len(results)
        print(f"  {direction}: {passed}/{total} success, {failed} failed, {skipped} skipped")
        for r in results:
            if r.get("status") == "failed":
                err = r.get("error", "")[:300].replace("\n", " | ")
                print(f"    FAIL: {r['test_id']} - {err}")
