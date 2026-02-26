#!/usr/bin/env python3
"""Check convert2anndata results from fig3_notebook_results.json"""
import json

with open('benchmark/results/fig3_notebook_results.json', encoding='utf-8') as f:
    data = json.load(f)

# Flatten nested structure
flat = {}
for key, val in data.items():
    if isinstance(val, dict) and any(d in val for d in ['rds_to_h5ad', 'h5ad_to_rds']):
        flat[key] = val
    elif isinstance(val, dict):
        for sk, sv in val.items():
            if isinstance(sv, dict):
                flat[sk] = sv

c2a = flat.get('convert2anndata', {})
for direction in ['rds_to_h5ad', 'h5ad_to_rds']:
    results = c2a.get(direction, [])
    if not results:
        print(f"{direction}: no data")
        continue
    ok = sum(1 for r in results if r.get('status') == 'success')
    fail = sum(1 for r in results if r.get('status') == 'failed')
    print(f"\n{direction}: {ok} success, {fail} failed, {len(results)} total")
    for r in results:
        tid = r.get('test_id', '?')
        status = r.get('status', '?')
        err = r.get('error', '')[:120]
        sv = r.get('seurat_version', '?')
        tp = r.get('type', '?')
        if status == 'failed':
            print(f"  FAIL  {tid:40s} {sv}/{tp}  {err}")
        else:
            t = r.get('conversion_time_seconds', 0)
            print(f"  OK    {tid:40s} {sv}/{tp}  {t:.1f}s")
