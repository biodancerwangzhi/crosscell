#!/usr/bin/env python3
"""
下载/导出 Python 生态 H5AD 数据集用于性能基准测试。

数据源:
  1. scanpy 内置数据集 (pbmc3k, pbmc3k_processed)
  2. scvelo 数据集 (pancreas, dentategyrus) — 直接下载 h5ad
  3. squidpy 空间数据集 (imc, seqfish, visium_hne, merfish, slideseqv2, mibitof)

使用方法 (在 benchmark 容器内):
  python3 /benchmark/data/05_download_h5ad_datasets.py --output-dir /benchmark/data/generated

  或分步执行:
  python3 /benchmark/data/05_download_h5ad_datasets.py --output-dir /benchmark/data/generated --source scanpy
  python3 /benchmark/data/05_download_h5ad_datasets.py --output-dir /benchmark/data/generated --source scvelo
  python3 /benchmark/data/05_download_h5ad_datasets.py --output-dir /benchmark/data/generated --source squidpy
"""

import argparse
import os
import subprocess
import sys
import json
from datetime import datetime


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def file_info(path):
    """返回文件大小 (MB)"""
    if os.path.exists(path):
        return os.path.getsize(path) / 1024 / 1024
    return 0


def curl_download(url, outpath):
    """用 curl 下载文件，支持断点续传和自动重试"""
    tmp_path = outpath + ".tmp"
    try:
        # 使用 -C - 支持断点续传，curl 会自动从上次中断处继续
        max_attempts = 30
        for attempt in range(1, max_attempts + 1):
            cmd = [
                "curl", "-L", "-C", "-",
                "--retry", "5",
                "--retry-delay", "3",
                "--connect-timeout", "30",
                "--max-time", "600",
                "-o", tmp_path,
                url,
            ]
            result = subprocess.run(cmd, timeout=900, capture_output=True, text=True)
            if result.returncode == 0:
                os.rename(tmp_path, outpath)
                return True
            # 如果是 SSL 错误 (56) 或超时 (28)，继续重试（断点续传）
            if result.returncode in (56, 28, 18) and attempt < max_attempts:
                print(f"    ⚠️  连接中断 (尝试 {attempt}/{max_attempts})，断点续传...")
                continue
            else:
                print(f"    ❌ curl 退出码 {result.returncode}: {result.stderr[-200:] if result.stderr else ''}")
                break

        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False
    except Exception as e:
        print(f"    ❌ curl 下载失败: {e}")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False


def download_and_verify(filename, url, desc, output_dir):
    """下载 h5ad 文件并验证可读性"""
    outpath = os.path.join(output_dir, filename)
    if os.path.exists(outpath):
        size = file_info(outpath)
        print(f"  ⏭️  已存在: {filename} ({size:.1f} MB)")
        return {"file": filename, "status": "skipped", "size_mb": round(size, 2)}

    print(f"  📥 下载 {desc}...")
    try:
        if not curl_download(url, outpath):
            raise RuntimeError(f"下载失败: {url}")

        import anndata
        adata = anndata.read_h5ad(outpath)
        size = file_info(outpath)
        print(f"  ✅ {filename}: {adata.n_obs} cells × {adata.n_vars} genes ({size:.1f} MB)")
        return {
            "file": filename,
            "status": "success",
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "size_mb": round(size, 2),
        }
    except Exception as e:
        print(f"  ❌ {filename}: {e}")
        if os.path.exists(outpath):
            os.remove(outpath)
        return {"file": filename, "status": "failed", "error": str(e)}


def download_scanpy_datasets(output_dir):
    """下载 scanpy 数据集"""
    results = []

    # pbmc3k raw — 直接下载
    results.append(download_and_verify(
        "scanpy_pbmc3k.h5ad",
        "https://falexwolf.de/data/pbmc3k_raw.h5ad",
        "pbmc3k raw (2700 cells × 32738 genes)",
        output_dir,
    ))

    # pbmc3k processed — 直接下载
    results.append(download_and_verify(
        "scanpy_pbmc3k_processed.h5ad",
        "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad",
        "pbmc3k processed (2638 cells × 1838 genes)",
        output_dir,
    ))

    return results


def download_scvelo_datasets(output_dir):
    """下载 scvelo 数据集"""
    datasets = [
        (
            "scvelo_pancreas.h5ad",
            "https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad",
            "pancreas endocrinogenesis (~3696 cells)",
        ),
        (
            "scvelo_dentategyrus.h5ad",
            "https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/10X43_1.h5ad",
            "dentate gyrus (~2930 cells)",
        ),
    ]
    return [download_and_verify(f, u, d, output_dir) for f, u, d in datasets]


def download_squidpy_datasets(output_dir):
    """下载 squidpy 空间数据集 (直接从 figshare 下载，无需安装 squidpy)"""
    # 从 squidpy 源码提取的 figshare 下载 URL
    datasets = [
        ("squidpy_imc.h5ad", "https://ndownloader.figshare.com/files/26098406", "IMC (4668 cells × 34 genes)"),
        ("squidpy_seqfish.h5ad", "https://ndownloader.figshare.com/files/26098403", "seqFISH"),
        ("squidpy_visium_hne.h5ad", "https://ndownloader.figshare.com/files/26098397", "Visium H&E"),
        ("squidpy_merfish.h5ad", "https://ndownloader.figshare.com/files/28169379", "MERFISH"),
        ("squidpy_slideseqv2.h5ad", "https://ndownloader.figshare.com/files/28242783", "Slide-seqV2"),
        ("squidpy_mibitof.h5ad", "https://ndownloader.figshare.com/files/28241139", "MIBI-TOF"),
    ]
    return [download_and_verify(f, u, d, output_dir) for f, u, d in datasets]


def download_cellxgene_datasets(output_dir):
    """从 CELLxGENE Discover 下载中/大规模 H5AD 数据集"""
    datasets = [
        (
            "cellxgene_pbmc_15k.h5ad",
            "https://datasets.cellxgene.cziscience.com/5f93ba67-b129-4bd7-9d46-df4804225843.h5ad",
            "PBMC scRNA-seq (~14783 cells, blood)",
        ),
        (
            "cellxgene_heart_23k.h5ad",
            "https://datasets.cellxgene.cziscience.com/6fa831d7-b8ab-41f1-9134-b25ff31263d3.h5ad",
            "Atrial cardiomyocytes (~23483 cells, heart)",
        ),
        (
            "cellxgene_brain_40k.h5ad",
            "https://datasets.cellxgene.cziscience.com/7ad15494-dfd0-43f7-b858-34be24f4df1f.h5ad",
            "Claustrum (~40165 cells, brain)",
        ),
    ]
    return [download_and_verify(f, u, d, output_dir) for f, u, d in datasets]


def main():
    parser = argparse.ArgumentParser(description="下载 H5AD 基准测试数据集")
    parser.add_argument("--output-dir", default="/benchmark/data/generated")
    parser.add_argument("--source", choices=["all", "scanpy", "scvelo", "squidpy", "cellxgene"], default="all")
    args = parser.parse_args()

    ensure_dir(args.output_dir)

    print("=" * 50)
    print("CrossCell Benchmark - H5AD 数据集下载")
    print("=" * 50)
    print(f"输出目录: {args.output_dir}")
    print(f"数据源: {args.source}")
    print()

    all_results = {}

    if args.source in ("all", "scanpy"):
        print("--- scanpy 数据集 ---")
        all_results["scanpy"] = download_scanpy_datasets(args.output_dir)
        print()

    if args.source in ("all", "scvelo"):
        print("--- scvelo 数据集 ---")
        all_results["scvelo"] = download_scvelo_datasets(args.output_dir)
        print()

    if args.source in ("all", "squidpy"):
        print("--- squidpy 空间数据集 ---")
        all_results["squidpy"] = download_squidpy_datasets(args.output_dir)
        print()

    if args.source in ("all", "cellxgene"):
        print("--- CELLxGENE 中/大规模数据集 ---")
        all_results["cellxgene"] = download_cellxgene_datasets(args.output_dir)
        print()

    # 保存下载结果摘要
    summary = {
        "timestamp": datetime.now().isoformat(),
        "output_dir": args.output_dir,
        "results": all_results,
    }
    summary_path = os.path.join(args.output_dir, "download_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # 打印汇总
    print("=" * 50)
    print("下载汇总")
    print("=" * 50)
    total = sum(len(v) for v in all_results.values())
    success = sum(1 for v in all_results.values() for r in v if r["status"] == "success")
    skipped = sum(1 for v in all_results.values() for r in v if r["status"] == "skipped")
    failed = sum(1 for v in all_results.values() for r in v if r["status"] == "failed")
    print(f"总计: {total}")
    print(f"✅ 成功: {success}")
    print(f"⏭️  跳过: {skipped}")
    print(f"❌ 失败: {failed}")
    print(f"\n结果已保存到: {summary_path}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
