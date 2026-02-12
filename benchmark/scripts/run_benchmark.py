#!/usr/bin/env python3
"""
CrossCell Performance Benchmark Runner

读取 test_matrix.json，对每个工具 × 数据集组合运行转换测试，
测量转换时间、峰值内存、输出文件大小。

容器内路径:
  /benchmark/data/generated/  — RDS + H5AD 数据文件
  /benchmark/config/          — test_matrix.json, datasets.json
  /benchmark/results/         — 输出结果 JSON
  /benchmark/scripts/         — 本脚本 (只读)

用法:
  python3 /benchmark/scripts/run_benchmark.py --tool crosscell
  python3 /benchmark/scripts/run_benchmark.py --tool zellkonverter --direction rds_to_h5ad
  python3 /benchmark/scripts/run_benchmark.py --tool all --scale small
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ============================================
# 容器内路径常量
# ============================================
DATA_DIR = "/benchmark/data/generated"
CONFIG_DIR = "/benchmark/config"
RESULTS_DIR = "/benchmark/results"
MATRIX_FILE = os.path.join(CONFIG_DIR, "test_matrix.json")
OUTPUT_TMP_DIR = "/tmp/benchmark_output"

# ============================================
# 工具转换命令定义
# ============================================

def get_output_path(test_id, tool, direction):
    """生成输出文件路径"""
    ext = "h5ad" if direction == "rds_to_h5ad" else "rds"
    return os.path.join(OUTPUT_TMP_DIR, f"{tool}_{test_id}.{ext}")


def cmd_crosscell(input_path, output_path, direction):
    """CrossCell CLI 转换命令"""
    fmt = "anndata" if direction == "rds_to_h5ad" else "seurat"
    return ["crosscell", "convert", "-i", input_path, "-o", output_path, "-f", fmt]


def cmd_zellkonverter(input_path, output_path, direction):
    """Zellkonverter (R) 转换命令
    RDS→H5AD: UpdateSeuratObject → as.SingleCellExperiment → writeH5AD (basilisk)
    H5AD→RDS: readH5AD(reader="R") → as.Seurat
      - 先尝试 reader="R" (纯 R，不需要 basilisk，依赖 rhdf5+HDF5Array)
      - 如果 as.Seurat 因 dimnames 缺失而失败，自动补 dimnames 后重试
      - 如果 R reader 本身失败，fallback 到 reticulate 直连 Python reader
    """
    if direction == "rds_to_h5ad":
        script = f"""
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
obj <- readRDS("{input_path}")
obj <- UpdateSeuratObject(obj)
sce <- as.SingleCellExperiment(obj)
writeH5AD(sce, "{output_path}")
"""
    else:
        # H5AD → RDS: reader="R" + dimnames 修复 + reticulate fallback
        script = f"""
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)

fix_dimnames <- function(sce) {{
  # 部分 H5AD (如 scanpy_pbmc3k) 用旧格式，R reader 解析后缺 dimnames
  # as.Seurat 需要 rownames (基因名) 和 colnames (细胞 barcode)
  if (is.null(rownames(sce)) || length(rownames(sce)) == 0) {{
    if ("_index" %in% names(rowData(sce))) {{
      rownames(sce) <- rowData(sce)[["_index"]]
    }} else if ("gene_ids" %in% names(rowData(sce))) {{
      rownames(sce) <- rowData(sce)[["gene_ids"]]
    }} else {{
      rownames(sce) <- paste0("Gene_", seq_len(nrow(sce)))
    }}
  }}
  if (is.null(colnames(sce)) || length(colnames(sce)) == 0) {{
    if ("_index" %in% names(colData(sce))) {{
      colnames(sce) <- colData(sce)[["_index"]]
    }} else if ("obs_names" %in% names(colData(sce))) {{
      colnames(sce) <- colData(sce)[["obs_names"]]
    }} else {{
      colnames(sce) <- paste0("Cell_", seq_len(ncol(sce)))
    }}
  }}
  sce
}}

convert_ok <- FALSE

# 方法 1: reader="R" (纯 R，不需要 basilisk)
tryCatch({{
  sce <- readH5AD("{input_path}", reader = "R")
  sce <- fix_dimnames(sce)
  obj <- as.Seurat(sce, counts = "X", data = NULL)
  saveRDS(obj, "{output_path}")
  convert_ok <- TRUE
}}, error = function(e) {{
  cat("R reader failed:", conditionMessage(e), "\\n")
}})

# 方法 2: fallback 到 reticulate 直连 Python (绕过 basilisk)
if (!convert_ok) {{
  cat("Falling back to reticulate Python reader...\\n")
  Sys.setenv(RETICULATE_PYTHON = "/opt/venv/bin/python3")
  library(reticulate)
  use_python("/opt/venv/bin/python3", required = TRUE)
  ad <- import("anndata")
  adata <- ad$read_h5ad("{input_path}")
  sce <- AnnData2SCE(adata)
  sce <- fix_dimnames(sce)
  obj <- as.Seurat(sce, counts = "X", data = NULL)
  saveRDS(obj, "{output_path}")
}}
"""
    return ["Rscript", "-e", script.strip()]


def cmd_anndatar(input_path, output_path, direction):
    """anndataR (R) 转换命令
    API (v1.0.0):
      - as_AnnData(seurat_obj) → AnnData R6 object
      - read_h5ad(path) → AnnData R6 object
      - ad$write_h5ad(path) / write_h5ad(ad, path)
      - ad$as_Seurat() → Seurat object
    Requires: rhdf5 (Bioconductor)
    """
    if direction == "rds_to_h5ad":
        script = f"""
library(anndataR)
library(Seurat)
seurat_obj <- readRDS("{input_path}")
seurat_obj <- UpdateSeuratObject(seurat_obj)
ad <- as_AnnData(seurat_obj)
write_h5ad(ad, "{output_path}")
"""
    else:
        script = f"""
library(anndataR)
library(Seurat)
ad <- read_h5ad("{input_path}")
seurat_obj <- ad$as_Seurat()
saveRDS(seurat_obj, "{output_path}")
"""
    return ["Rscript", "-e", script.strip()]


def cmd_convert2anndata(input_path, output_path, direction):
    """convert2anndata (R) — 只支持 RDS→H5AD
    API (settylab/convert2anndata):
      - convert_seurat_to_sce(seurat_obj) → SingleCellExperiment
      - convert_to_anndata(sce) → AnnData (via reticulate)
      - anndata::write_h5ad(ad, path)
    """
    if direction != "rds_to_h5ad":
        return None  # 不支持反向转换
    script = f"""
library(convert2anndata)
library(Seurat)
obj <- readRDS("{input_path}")
obj <- UpdateSeuratObject(obj)
sce <- convert_seurat_to_sce(obj)
ad <- convert_to_anndata(sce)
anndata::write_h5ad(ad, "{output_path}")
"""
    return ["Rscript", "-e", script.strip()]


def cmd_easyscf(input_path, output_path, direction):
    """easySCF 转换命令
    使用 R 包 (easySCFr)，因为它直接处理 Seurat 对象。
    API:
      - saveH5(seurat_obj, path) — Seurat → H5 (easySCF 自有格式，非标准 H5AD)
      - loadH5(path) → Seurat object (读取 easySCF H5 格式)
    注意: easySCF 使用自己的 H5 格式，不是标准 H5AD。
          RDS→H5AD 方向: 输出为 easySCF H5 格式
          H5AD→RDS 方向: 不可用。虽然理论上可以 Python saveH5 → R loadH5 中转，
          但 easySCFpy 的 saveH5 使用 anndata._io.specs.write_elem 写 var DataFrame，
          R 端 loadH5 读取 var/rawvar/_index 时因 anndata 编码版本不兼容而失败
          (notSubsettableError: object of type 'environment' is not subsettable)
    """
    if direction == "rds_to_h5ad":
        script = f"""
library(easySCFr)
library(Seurat)
obj <- readRDS("{input_path}")
obj <- UpdateSeuratObject(obj)
saveH5(obj, "{output_path}")
"""
        return ["Rscript", "-e", script.strip()]
    else:
        # easySCF h5ad_to_rds 不可用：Python saveH5 写的 var 编码与 R loadH5 不兼容
        return None


TOOL_COMMANDS = {
    "crosscell": cmd_crosscell,
    "zellkonverter": cmd_zellkonverter,
    "anndataR": cmd_anndatar,
    "convert2anndata": cmd_convert2anndata,
    "easySCF": cmd_easyscf,
}

# ============================================
# 性能测量
# ============================================

def split_time_stderr(stderr):
    """
    将 /usr/bin/time -v 的统计信息从 stderr 中分离出来。
    /usr/bin/time -v 的输出以 "Command being timed:" 开头。
    返回 (tool_stderr, time_stats)
    """
    marker = "\tCommand being timed:"
    idx = stderr.find(marker)
    if idx >= 0:
        return stderr[:idx].rstrip(), stderr[idx:]
    return stderr, ""


def run_with_measurement(cmd, timeout=600):
    """
    运行命令并测量时间和峰值内存。
    使用 /usr/bin/time -v 获取峰值 RSS。
    返回 (success, time_sec, peak_mem_mb, tool_stderr)
    """
    time_cmd = ["/usr/bin/time", "-v"] + cmd

    start = time.time()
    try:
        result = subprocess.run(
            time_cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.time() - start
    except subprocess.TimeoutExpired:
        return False, timeout, 0, "TIMEOUT"

    # 分离工具 stderr 和 /usr/bin/time 统计
    tool_stderr, time_stats = split_time_stderr(result.stderr)

    # 从 /usr/bin/time -v 的统计中提取峰值内存
    peak_mem_mb = 0
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", time_stats)
    if match:
        peak_mem_mb = int(match.group(1)) / 1024.0

    success = result.returncode == 0
    # 失败时返回工具自身的 stderr + stdout（有些工具把错误写到 stdout）
    error_output = tool_stderr
    if not success and result.stdout:
        error_output = result.stdout.rstrip() + "\n" + tool_stderr
    return success, elapsed, peak_mem_mb, error_output


def get_file_size_mb(path):
    """获取文件大小 (MB)"""
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return 0


# ============================================
# 主逻辑
# ============================================

def load_matrix(matrix_path):
    """加载测试矩阵"""
    with open(matrix_path) as f:
        return json.load(f)


def resolve_input_path(raw_path):
    """
    将 test_matrix.json 中的相对路径转换为容器内绝对路径。
    例: data/generated/seurat_v4_pbmc3k_raw.rds → /benchmark/data/generated/seurat_v4_pbmc3k_raw.rds
    """
    # 去掉 "data/" 前缀，因为容器内 /benchmark/data 对应 host 的 data/
    if raw_path.startswith("data/"):
        return os.path.join("/benchmark", raw_path)
    return raw_path


def filter_cases(cases, scale=None):
    """按规模过滤测试用例"""
    if scale and scale != "all":
        cases = [c for c in cases if c.get("scale") == scale]
    return cases


def check_tool_supports(matrix, tool_name, direction):
    """检查工具是否支持指定方向"""
    for tool in matrix["tools"]:
        if tool["name"] == tool_name:
            return direction in tool.get("directions", [])
    return False


def run_single_test(tool_name, test_case, direction, runs=3, timeout=600):
    """
    对单个测试用例运行多次，返回中位数结果。
    """
    test_id = test_case["id"]
    input_path = resolve_input_path(test_case["input"])
    output_path = get_output_path(test_id, tool_name, direction)

    # 检查输入文件存在
    if not os.path.exists(input_path):
        return {
            "test_id": test_id,
            "status": "skipped",
            "reason": f"input not found: {input_path}",
        }

    # 获取命令生成函数
    cmd_fn = TOOL_COMMANDS.get(tool_name)
    if not cmd_fn:
        return {
            "test_id": test_id,
            "status": "error",
            "reason": f"unknown tool: {tool_name}",
        }

    cmd = cmd_fn(input_path, output_path, direction)
    if cmd is None:
        return {
            "test_id": test_id,
            "status": "skipped",
            "reason": f"{tool_name} does not support {direction}",
        }

    times = []
    mems = []
    last_stderr = ""

    for i in range(runs):
        # 清理上次输出
        if os.path.exists(output_path):
            os.remove(output_path)

        success, elapsed, peak_mem, stderr = run_with_measurement(cmd, timeout)
        last_stderr = stderr

        if not success:
            # 如果任何一次运行失败，记录失败并停止
            error_msg = last_stderr[:5000] if len(last_stderr) > 5000 else last_stderr
            return {
                "test_id": test_id,
                "dataset": test_case.get("dataset", ""),
                "direction": direction,
                "status": "failed",
                "error": error_msg,
                "run": i + 1,
            }

        times.append(elapsed)
        mems.append(peak_mem)

    # 取中位数
    times.sort()
    mems.sort()
    median_time = times[len(times) // 2]
    median_mem = mems[len(mems) // 2]

    # 输出文件大小
    output_size = get_file_size_mb(output_path)

    return {
        "test_id": test_id,
        "dataset": test_case.get("dataset", ""),
        "direction": direction,
        "status": "success",
        "conversion_time_seconds": round(median_time, 3),
        "peak_memory_mb": round(median_mem, 1),
        "output_file_size_mb": round(output_size, 2),
        "all_times": [round(t, 3) for t in times],
        "all_mems": [round(m, 1) for m in mems],
        "n_cells": test_case.get("n_cells"),
        "n_genes": test_case.get("n_genes"),
        "scale": test_case.get("scale"),
        "seurat_version": test_case.get("seurat_version", ""),
        "type": test_case.get("type", ""),
    }


def run_benchmark(tool_name, matrix, direction="all", scale=None, runs=3, timeout=600):
    """运行指定工具的完整基准测试"""
    results = {
        "tool": tool_name,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": {
            "runs_per_test": runs,
            "timeout_seconds": timeout,
            "scale_filter": scale or "all",
            "direction_filter": direction,
        },
        "rds_to_h5ad": [],
        "h5ad_to_rds": [],
        "summary": {},
    }

    directions = []
    if direction in ("all", "rds_to_h5ad"):
        directions.append("rds_to_h5ad")
    if direction in ("all", "h5ad_to_rds"):
        directions.append("h5ad_to_rds")

    total = 0
    success_count = 0
    failed_count = 0
    skipped_count = 0

    for d in directions:
        if not check_tool_supports(matrix, tool_name, d):
            print(f"  ⏭️  {tool_name} 不支持 {d}，跳过")
            continue

        cases = filter_cases(matrix["test_cases"].get(d, []), scale)
        print(f"\n{'='*60}")
        print(f"  {tool_name} — {d} ({len(cases)} 个测试用例)")
        print(f"{'='*60}")

        for i, case in enumerate(cases):
            total += 1
            test_id = case["id"]
            n_cells = case.get("n_cells", "?")
            n_genes = case.get("n_genes", "?")
            print(f"\n  [{i+1}/{len(cases)}] {test_id}")
            print(f"    {n_cells} cells × {n_genes} genes ({case.get('scale', '?')})")

            result = run_single_test(tool_name, case, d, runs, timeout)

            if result["status"] == "success":
                success_count += 1
                t = result["conversion_time_seconds"]
                m = result["peak_memory_mb"]
                s = result["output_file_size_mb"]
                print(f"    ✅ {t:.2f}s | {m:.0f}MB | {s:.1f}MB output")
            elif result["status"] == "skipped":
                skipped_count += 1
                print(f"    ⏭️  跳过: {result.get('reason', '')}")
            else:
                failed_count += 1
                print(f"    ❌ 失败: {result.get('error', '')[:500]}")

            results[d].append(result)

    results["summary"] = {
        "total": total,
        "success": success_count,
        "failed": failed_count,
        "skipped": skipped_count,
    }

    return results


def main():
    parser = argparse.ArgumentParser(description="CrossCell Performance Benchmark Runner")
    parser.add_argument(
        "--tool",
        required=True,
        choices=list(TOOL_COMMANDS.keys()) + ["all"],
        help="要测试的工具名称，或 'all' 测试所有工具",
    )
    parser.add_argument(
        "--direction",
        default="all",
        choices=["rds_to_h5ad", "h5ad_to_rds", "all"],
        help="转换方向 (默认: all)",
    )
    parser.add_argument(
        "--scale",
        default=None,
        choices=["small", "medium", "large"],
        help="只测试指定规模的数据集",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=3,
        help="每个测试运行次数，取中位数 (默认: 3)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="单次运行超时秒数 (默认: 600)",
    )
    parser.add_argument(
        "--matrix",
        default=MATRIX_FILE,
        help=f"test_matrix.json 路径 (默认: {MATRIX_FILE})",
    )
    parser.add_argument(
        "--output-dir",
        default=RESULTS_DIR,
        help=f"结果输出目录 (默认: {RESULTS_DIR})",
    )
    args = parser.parse_args()

    # 加载测试矩阵
    if not os.path.exists(args.matrix):
        print(f"❌ 找不到测试矩阵: {args.matrix}")
        sys.exit(1)

    matrix = load_matrix(args.matrix)
    print(f"📋 加载测试矩阵: {matrix['summary']['total_test_cases']} 个测试用例")

    # 创建临时输出目录
    os.makedirs(OUTPUT_TMP_DIR, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    # 确定要测试的工具列表
    if args.tool == "all":
        tools = list(TOOL_COMMANDS.keys())
    else:
        tools = [args.tool]

    for tool_name in tools:
        print(f"\n{'#'*60}")
        print(f"# 开始测试: {tool_name}")
        print(f"{'#'*60}")

        results = run_benchmark(
            tool_name, matrix, args.direction, args.scale, args.runs, args.timeout
        )

        # 保存结果（如果只跑了单方向，合并已有结果）
        output_file = os.path.join(args.output_dir, f"{tool_name}.json")
        if args.direction != "all" and os.path.exists(output_file):
            try:
                with open(output_file, "r") as f:
                    existing = json.load(f)
                # 合并：保留已有方向的结果，更新新跑的方向
                for d in ["rds_to_h5ad", "h5ad_to_rds"]:
                    if not results[d] and existing.get(d):
                        results[d] = existing[d]
                # 重新计算 summary
                total = len(results["rds_to_h5ad"]) + len(results["h5ad_to_rds"])
                success = sum(1 for r in results["rds_to_h5ad"] + results["h5ad_to_rds"] if r["status"] == "success")
                failed = sum(1 for r in results["rds_to_h5ad"] + results["h5ad_to_rds"] if r["status"] == "failed")
                skipped = sum(1 for r in results["rds_to_h5ad"] + results["h5ad_to_rds"] if r["status"] == "skipped")
                results["summary"] = {"total": total, "success": success, "failed": failed, "skipped": skipped}
            except (json.JSONDecodeError, KeyError):
                pass  # 已有文件损坏，直接覆盖

        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"\n💾 结果已保存: {output_file}")

        # 打印汇总
        s = results["summary"]
        print(f"\n📊 {tool_name} 汇总: {s['success']}/{s['total']} 成功, "
              f"{s['failed']} 失败, {s['skipped']} 跳过")

    print(f"\n{'='*60}")
    print("✅ 基准测试完成!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
