#!/usr/bin/env python3
"""Show full error for a specific failed test."""
import json
import sys

tool = sys.argv[1] if len(sys.argv) > 1 else "zellkonverter"
direction = sys.argv[2] if len(sys.argv) > 2 else "rds_to_h5ad"
idx = int(sys.argv[3]) if len(sys.argv) > 3 else 0

with open(f"benchmark/results/{tool}.json", encoding="utf-8") as f:
    data = json.load(f)

failed = [r for r in data.get(direction, []) if r.get("status") == "failed"]
if idx >= len(failed):
    print(f"Only {len(failed)} failed tests")
    sys.exit(1)

r = failed[idx]
print(f"=== {r['test_id']} ({direction}) ===")
print(f"Status: {r['status']}")
err = r.get("error", "")
# Print last 2000 chars to see the actual error
print("--- ERROR (last 2000 chars) ---")
print(err[-2000:])
