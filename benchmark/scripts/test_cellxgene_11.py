#!/usr/bin/env python3
"""
系统测试 11 个 CELLxGENE benchmark 数据集
测试内容：
  1. crosscell inspect — 获取数据集详细信息
  2. crosscell convert H5AD→RDS — 测试转换能力和性能
  3. crosscell convert RDS→H5AD (roundtrip) — 往返转换测试

在 benchmark 容器中运行：
  python3 /benchmark/scripts/test_cellxgene_11.py
"""

import subprocess
import time
import json
import os
import sys
import traceback

DATA_DIR = "/benchmark/data/generated"
OUTPUT_DIR = "/benchmark/results/cellxgene_11"

# 11 个数据集，按规模从小到大排列
DATASETS = [
    ("cellxgene_skin_bcc_10k.h5ad", "~10k cells, skin/tumor"),
    ("cellxgene_tabula_liver_22k.h5ad", "~22k cells, liver"),
    ("cellxgene_kidney_atacseq_37k.h5ad", "~37k cells, kidney scATAC"),
    ("cellxgene_brain_multiome_102k.h5ad", "~102k cells, brain multiome"),
    ("cellxgene_pancreas_122k.h5ad", "~122k cells, pancreas"),
    ("cellxgene_brain_dlpfc_172k.h5ad", "~172k cells, brain dlPFC"),
    ("cellxgene_retina_244k.h5ad", "~244k cells, retina"),
    ("cellxgene_gut_428k.h5ad", "~428k cells, gut"),
    ("cellxgene_heart_486k.h5ad", "~486k cells, heart"),
    ("cellxgene_hlca_core_585k.h5ad", "~585k cells, lung HLCA"),
    ("cellxgene_combat_pbmc_836k.h5ad", "~836k cells, blood COMBAT"),
]


def run_cmd(cmd, timeout=1800):
    """运行命令，返回 (returncode, stdout, stderr, elapsed_seconds)"""
    start = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )
        elapsed = time.time() - start
        return result.returncode, result.stdout, result.stderr, elapsed
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        return -1, "", f"TIMEOUT after {timeout}s", elapsed
    except Exception as e:
        elapsed = time.time() - start
        return -2, "", str(e), elapsed


def test_inspect(h5ad_path, name):
    """测试 crosscell inspect"""
    print(f"  [inspect] {name}...")
    cmd = ["crosscell", "inspect", "-i", h5ad_path, "--detailed"]
    rc, stdout, stderr, elapsed = run_cmd(cmd, timeout=300)
    
    result = {
        "test": "inspect",
        "returncode": rc,
        "elapsed_seconds": round(elapsed, 2),
        "success": rc == 0,
    }
    
    if rc == 0:
        result["output"] = stdout[:2000]  # 截断保存
        print(f"    ✅ 成功 ({elapsed:.1f}s)")
        # 打印关键信息
        for line in stdout.split("\n"):
            if any(k in line.lower() for k in ["cells", "genes", "sparse", "size", "format", "层", "细胞", "基因"]):
                print(f"    {line.strip()}")
    else:
        result["error"] = stderr[:1000]
        print(f"    ❌ 失败 (rc={rc}, {elapsed:.1f}s)")
        print(f"    {stderr[:300]}")
    
    return result


def test_convert_h5ad_to_rds(h5ad_path, rds_path, name):
    """测试 H5AD → RDS 转换"""
    print(f"  [h5ad→rds] {name}...")
    
    # 删除已有输出
    if os.path.exists(rds_path):
        os.remove(rds_path)
    
    cmd = ["crosscell", "convert", "-i", h5ad_path, "-o", rds_path, "-f", "seurat", "--verbose"]
    rc, stdout, stderr, elapsed = run_cmd(cmd, timeout=1800)
    
    result = {
        "test": "h5ad_to_rds",
        "returncode": rc,
        "elapsed_seconds": round(elapsed, 2),
        "success": rc == 0,
    }
    
    if rc == 0:
        rds_size = os.path.getsize(rds_path) / 1024 / 1024
        result["output_size_mb"] = round(rds_size, 2)
        print(f"    ✅ 成功 ({elapsed:.1f}s, output={rds_size:.1f} MB)")
        # 提取关键性能信息
        for line in (stdout + stderr).split("\n"):
            if any(k in line.lower() for k in ["peak memory", "time", "cells", "genes", "converted"]):
                print(f"    {line.strip()}")
    else:
        result["error"] = (stderr + stdout)[:1000]
        print(f"    ❌ 失败 (rc={rc}, {elapsed:.1f}s)")
        print(f"    {(stderr + stdout)[:500]}")
    
    return result


def test_convert_rds_to_h5ad(rds_path, h5ad_out_path, name):
    """测试 RDS → H5AD 往返转换"""
    if not os.path.exists(rds_path):
        print(f"  [rds→h5ad] {name}... ⏭️ 跳过（无 RDS 输入）")
        return {"test": "rds_to_h5ad", "success": False, "error": "no rds input"}
    
    print(f"  [rds→h5ad] {name}...")
    
    if os.path.exists(h5ad_out_path):
        os.remove(h5ad_out_path)
    
    cmd = ["crosscell", "convert", "-i", rds_path, "-o", h5ad_out_path, "-f", "anndata", "--verbose"]
    rc, stdout, stderr, elapsed = run_cmd(cmd, timeout=1800)
    
    result = {
        "test": "rds_to_h5ad",
        "returncode": rc,
        "elapsed_seconds": round(elapsed, 2),
        "success": rc == 0,
    }
    
    if rc == 0:
        h5ad_size = os.path.getsize(h5ad_out_path) / 1024 / 1024
        result["output_size_mb"] = round(h5ad_size, 2)
        print(f"    ✅ 成功 ({elapsed:.1f}s, output={h5ad_size:.1f} MB)")
    else:
        result["error"] = (stderr + stdout)[:1000]
        print(f"    ❌ 失败 (rc={rc}, {elapsed:.1f}s)")
        print(f"    {(stderr + stdout)[:500]}")
    
    return result


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("CrossCell 系统测试 — 11 个 CELLxGENE Benchmark 数据集")
    print("=" * 70)
    print(f"数据目录: {DATA_DIR}")
    print(f"输出目录: {OUTPUT_DIR}")
    print()
    
    all_results = {}
    summary = {"total": 0, "inspect_pass": 0, "h2r_pass": 0, "r2h_pass": 0, "inspect_fail": 0, "h2r_fail": 0, "r2h_fail": 0}
    
    for filename, desc in DATASETS:
        h5ad_path = os.path.join(DATA_DIR, filename)
        base = filename.replace(".h5ad", "")
        rds_path = os.path.join(OUTPUT_DIR, f"{base}.rds")
        h5ad_roundtrip_path = os.path.join(OUTPUT_DIR, f"{base}_roundtrip.h5ad")
        
        print(f"\n{'─' * 70}")
        print(f"📊 {filename} ({desc})")
        print(f"   输入大小: {os.path.getsize(h5ad_path) / 1024 / 1024:.1f} MB")
        print(f"{'─' * 70}")
        
        if not os.path.exists(h5ad_path):
            print(f"  ❌ 文件不存在，跳过")
            all_results[base] = {"error": "file not found"}
            continue
        
        summary["total"] += 1
        dataset_results = {"filename": filename, "desc": desc}
        dataset_results["input_size_mb"] = round(os.path.getsize(h5ad_path) / 1024 / 1024, 2)
        
        # 1. Inspect
        try:
            r = test_inspect(h5ad_path, base)
            dataset_results["inspect"] = r
            if r["success"]:
                summary["inspect_pass"] += 1
            else:
                summary["inspect_fail"] += 1
        except Exception as e:
            dataset_results["inspect"] = {"success": False, "error": str(e)}
            summary["inspect_fail"] += 1
            traceback.print_exc()
        
        # 2. H5AD → RDS
        try:
            r = test_convert_h5ad_to_rds(h5ad_path, rds_path, base)
            dataset_results["h5ad_to_rds"] = r
            if r["success"]:
                summary["h2r_pass"] += 1
            else:
                summary["h2r_fail"] += 1
        except Exception as e:
            dataset_results["h5ad_to_rds"] = {"success": False, "error": str(e)}
            summary["h2r_fail"] += 1
            traceback.print_exc()
        
        # 3. RDS → H5AD (roundtrip)
        try:
            r = test_convert_rds_to_h5ad(rds_path, h5ad_roundtrip_path, base)
            dataset_results["rds_to_h5ad"] = r
            if r["success"]:
                summary["r2h_pass"] += 1
            else:
                summary["r2h_fail"] += 1
        except Exception as e:
            dataset_results["rds_to_h5ad"] = {"success": False, "error": str(e)}
            summary["r2h_fail"] += 1
            traceback.print_exc()
        
        all_results[base] = dataset_results
        
        # 清理 roundtrip 文件节省磁盘（保留 RDS）
        if os.path.exists(h5ad_roundtrip_path):
            rt_size = os.path.getsize(h5ad_roundtrip_path) / 1024 / 1024
            os.remove(h5ad_roundtrip_path)
            print(f"  🗑️ 已清理 roundtrip 文件 ({rt_size:.1f} MB)")
    
    # 汇总
    print(f"\n{'=' * 70}")
    print("📋 测试汇总")
    print(f"{'=' * 70}")
    print(f"数据集总数: {summary['total']}")
    print(f"Inspect:    ✅ {summary['inspect_pass']} / ❌ {summary['inspect_fail']}")
    print(f"H5AD→RDS:   ✅ {summary['h2r_pass']} / ❌ {summary['h2r_fail']}")
    print(f"RDS→H5AD:   ✅ {summary['r2h_pass']} / ❌ {summary['r2h_fail']}")
    
    # 性能汇总表
    print(f"\n{'─' * 70}")
    print(f"{'数据集':<40} {'大小(MB)':>10} {'H2R(s)':>10} {'R2H(s)':>10} {'RDS(MB)':>10}")
    print(f"{'─' * 70}")
    for base, res in all_results.items():
        if "error" in res and isinstance(res["error"], str):
            continue
        size = res.get("input_size_mb", 0)
        h2r_time = res.get("h5ad_to_rds", {}).get("elapsed_seconds", "—")
        r2h_time = res.get("rds_to_h5ad", {}).get("elapsed_seconds", "—")
        rds_size = res.get("h5ad_to_rds", {}).get("output_size_mb", "—")
        h2r_ok = "✅" if res.get("h5ad_to_rds", {}).get("success") else "❌"
        r2h_ok = "✅" if res.get("rds_to_h5ad", {}).get("success") else "❌"
        
        h2r_str = f"{h2r_time}" if isinstance(h2r_time, str) else f"{h2r_time:.1f}"
        r2h_str = f"{r2h_time}" if isinstance(r2h_time, str) else f"{r2h_time:.1f}"
        rds_str = f"{rds_size}" if isinstance(rds_size, str) else f"{rds_size:.1f}"
        
        print(f"{base:<40} {size:>10.1f} {h2r_ok}{h2r_str:>8} {r2h_ok}{r2h_str:>8} {rds_str:>10}")
    
    # 保存结果
    output_file = os.path.join(OUTPUT_DIR, "test_results.json")
    with open(output_file, "w") as f:
        json.dump({"summary": summary, "datasets": all_results}, f, indent=2, default=str)
    print(f"\n结果已保存: {output_file}")


if __name__ == "__main__":
    main()
