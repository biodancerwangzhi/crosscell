#!/usr/bin/env python3
"""Compare crosscell vs zellkonverter performance."""
import json

tools = {}
for name in ["crosscell", "zellkonverter"]:
    with open(f"benchmark/results/{name}.json", encoding="utf-8") as f:
        tools[name] = json.load(f)

for direction in ["rds_to_h5ad", "h5ad_to_rds"]:
    print(f"\n{'='*80}")
    print(f"  {direction}")
    print(f"{'='*80}")
    print(f"{'Dataset':<30} {'CrossCell':>12} {'Zellkonverter':>14} {'Speedup':>10} {'CC Mem':>10} {'ZK Mem':>10} {'Mem Ratio':>10}")
    print("-" * 96)

    cc_results = {r["test_id"]: r for r in tools["crosscell"].get(direction, []) if r["status"] == "success"}
    zk_results = {r["test_id"]: r for r in tools["zellkonverter"].get(direction, []) if r["status"] == "success"}

    all_ids = sorted(set(list(cc_results.keys()) + list(zk_results.keys())),
                     key=lambda x: cc_results.get(x, zk_results.get(x, {})).get("n_cells", 0))

    speedups = []
    mem_ratios = []

    for tid in all_ids:
        cc = cc_results.get(tid)
        zk = zk_results.get(tid)
        if not cc or not zk:
            status = "CC only" if cc else "ZK only/fail"
            print(f"{tid:<30} {'--':>12} {'--':>14} {status:>10}")
            continue

        cc_time = cc["conversion_time_seconds"]
        zk_time = zk["conversion_time_seconds"]
        speedup = zk_time / cc_time if cc_time > 0 else 0
        speedups.append(speedup)

        cc_mem = cc["peak_memory_mb"]
        zk_mem = zk["peak_memory_mb"]
        mem_ratio = zk_mem / cc_mem if cc_mem > 0 else 0
        mem_ratios.append(mem_ratio)

        print(f"{tid:<30} {cc_time:>10.2f}s {zk_time:>12.2f}s {speedup:>9.1f}x {cc_mem:>8.0f}MB {zk_mem:>8.0f}MB {mem_ratio:>9.1f}x")

    if speedups:
        avg_speedup = sum(speedups) / len(speedups)
        avg_mem = sum(mem_ratios) / len(mem_ratios)
        min_speedup = min(speedups)
        max_speedup = max(speedups)
        print("-" * 96)
        print(f"{'AVERAGE':<30} {'':>12} {'':>14} {avg_speedup:>9.1f}x {'':>10} {'':>10} {avg_mem:>9.1f}x")
        print(f"{'RANGE':<30} {'':>12} {'':>14} {min_speedup:.1f}-{max_speedup:.1f}x")
