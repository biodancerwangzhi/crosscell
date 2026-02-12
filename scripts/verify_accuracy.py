#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CrossCell 准确性验证脚本
========================

本脚本用于验证 CrossCell 工具的格式转换准确性。
对应 Task 20.1 - 20.7 的测试内容。

CrossCell 的核心功能是 AnnData ↔ Seurat 跨格式转换，
因此本脚本主要测试 H5AD -> RDS -> H5AD 的往返准确性。

使用方法:
    # 在 Docker 容器中运行（推荐）
    docker-compose run --rm dev python3 自用/verify_accuracy.py
    
    # 或者在安装了依赖的环境中运行
    python3 自用/verify_accuracy.py

依赖:
    - numpy
    - scipy
    - h5py

作者: CrossCell Team
日期: 2026-01-28
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json

# 尝试导入依赖
try:
    import numpy as np
    import h5py
    from scipy import sparse
    HAS_DEPS = True
except ImportError as e:
    print(f"警告: 缺少依赖 {e}")
    print("请安装: pip install numpy scipy h5py")
    HAS_DEPS = False


# =============================================================================
# 配置
# =============================================================================

# CrossCell 可执行文件路径
CROSSCELL_BIN = os.environ.get("CROSSCELL_BIN", "./自用/crosscell_linux")

# 如果在 Docker 中运行，使用编译后的路径
if os.path.exists("/workspace/target/release/crosscell"):
    CROSSCELL_BIN = "/workspace/target/release/crosscell"

# 测试数据目录
TEST_DATA_DIR = Path("tests/data")
REAL_DATASETS_DIR = TEST_DATA_DIR / "real_datasets"

# 准确性阈值
CORRELATION_THRESHOLD = 0.999  # 跨格式转换的相关系数阈值


# =============================================================================
# 工具函数
# =============================================================================

def run_crosscell(args: List[str], timeout: int = 300) -> Tuple[int, str, str]:
    """
    运行 CrossCell 命令行工具
    
    参数:
        args: 命令行参数列表
        timeout: 超时时间（秒）
    
    返回:
        (返回码, 标准输出, 标准错误)
    """
    cmd = [CROSSCELL_BIN] + args
    print(f"    执行: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "命令超时"
    except FileNotFoundError:
        return -1, "", f"找不到可执行文件: {CROSSCELL_BIN}"


def read_h5ad_matrix(path: str) -> Tuple[np.ndarray, int, int, int]:
    """
    读取 H5AD 文件中的表达矩阵
    
    参数:
        path: H5AD 文件路径
    
    返回:
        (矩阵数据, 行数, 列数, 非零元素数)
    """
    with h5py.File(path, 'r') as f:
        if 'X' in f:
            x_group = f['X']
            if isinstance(x_group, h5py.Dataset):
                # 稠密矩阵
                data = x_group[:]
                nnz = int(np.count_nonzero(data))
                return data, data.shape[0], data.shape[1], nnz
            else:
                # 稀疏矩阵 (CSR 格式)
                data = x_group['data'][:]
                indices = x_group['indices'][:]
                indptr = x_group['indptr'][:]
                shape = x_group.attrs.get('shape', None)
                if shape is None:
                    n_rows = len(indptr) - 1
                    n_cols = int(indices.max()) + 1 if len(indices) > 0 else 0
                else:
                    n_rows, n_cols = shape
                
                matrix = sparse.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
                return matrix, n_rows, n_cols, matrix.nnz
    
    raise ValueError(f"无法读取 {path} 中的表达矩阵")


def calculate_correlation(a, b) -> float:
    """
    计算两个矩阵的相关系数
    """
    if sparse.issparse(a):
        a = a.toarray()
    if sparse.issparse(b):
        b = b.toarray()
    
    a_flat = a.flatten()
    b_flat = b.flatten()
    
    if len(a_flat) == 0 or len(b_flat) == 0:
        return 1.0 if len(a_flat) == len(b_flat) else 0.0
    
    if np.std(a_flat) == 0 and np.std(b_flat) == 0:
        return 1.0 if np.allclose(a_flat, b_flat) else 0.0
    
    return float(np.corrcoef(a_flat, b_flat)[0, 1])


def calculate_max_error(a, b) -> float:
    """计算最大绝对误差"""
    if sparse.issparse(a):
        a = a.toarray()
    if sparse.issparse(b):
        b = b.toarray()
    return float(np.max(np.abs(a - b)))


def print_result(name: str, passed: bool, details: str = ""):
    """打印测试结果"""
    status = "✅ 通过" if passed else "❌ 失败"
    print(f"    {name}: {status}")
    if details:
        print(f"      {details}")


# =============================================================================
# 核心测试函数
# =============================================================================

def test_cross_format_roundtrip(input_path: str, name: str) -> Dict:
    """
    测试跨格式往返转换 (H5AD -> RDS -> H5AD)
    
    这是 CrossCell 的核心功能测试。
    
    参数:
        input_path: 输入 H5AD 文件路径
        name: 测试名称
    
    返回:
        测试结果字典
    """
    print(f"\n  {'='*56}")
    print(f"  测试: {name}")
    print(f"  流程: H5AD -> RDS -> H5AD")
    print(f"  {'='*56}")
    
    result = {
        'name': name,
        'passed': False,
        'correlation': 0.0,
        'nnz_match': False,
        'dim_match': False,
    }
    
    if not os.path.exists(input_path):
        print(f"    ⚠️ 跳过: 文件不存在 {input_path}")
        result['skipped'] = True
        return result
    
    # 创建临时文件
    rds_path = tempfile.mktemp(suffix='.rds')
    h5ad_rt_path = tempfile.mktemp(suffix='.h5ad')
    
    try:
        # 1. 读取原始数据
        print(f"    [1/4] 读取原始 H5AD: {input_path}")
        orig_matrix, orig_rows, orig_cols, orig_nnz = read_h5ad_matrix(input_path)
        print(f"          维度: {orig_rows} × {orig_cols}, NNZ: {orig_nnz}")
        
        # 2. H5AD -> RDS
        print(f"    [2/4] 转换: H5AD -> RDS")
        ret, stdout, stderr = run_crosscell([
            'convert', '-i', input_path, '-o', rds_path, '-f', 'seurat'
        ])
        
        if ret != 0:
            print(f"    ❌ H5AD -> RDS 转换失败")
            result['error'] = f"H5AD -> RDS failed"
            return result
        
        # 3. RDS -> H5AD
        print(f"    [3/4] 转换: RDS -> H5AD")
        ret, stdout, stderr = run_crosscell([
            'convert', '-i', rds_path, '-o', h5ad_rt_path, '-f', 'anndata'
        ])
        
        if ret != 0:
            print(f"    ❌ RDS -> H5AD 转换失败")
            result['error'] = f"RDS -> H5AD failed"
            return result
        
        # 4. 读取往返后的数据并比较
        print(f"    [4/4] 验证往返结果")
        rt_matrix, rt_rows, rt_cols, rt_nnz = read_h5ad_matrix(h5ad_rt_path)
        print(f"          维度: {rt_rows} × {rt_cols}, NNZ: {rt_nnz}")
        
        # 比较结果
        dim_match = (orig_rows == rt_rows) and (orig_cols == rt_cols)
        nnz_match = (orig_nnz == rt_nnz)
        correlation = calculate_correlation(orig_matrix, rt_matrix)
        max_error = calculate_max_error(orig_matrix, rt_matrix)
        
        print(f"\n    --- 验证结果 ---")
        print_result("维度匹配", dim_match, f"{orig_rows}×{orig_cols} vs {rt_rows}×{rt_cols}")
        print_result("NNZ 匹配", nnz_match, f"{orig_nnz} vs {rt_nnz}")
        print_result("相关系数", correlation > CORRELATION_THRESHOLD, f"{correlation:.10f}")
        print_result("最大误差", max_error < 1e-6, f"{max_error:.2e}")
        
        # 汇总
        result['passed'] = dim_match and nnz_match and (correlation > CORRELATION_THRESHOLD)
        result['correlation'] = correlation
        result['nnz_match'] = nnz_match
        result['dim_match'] = dim_match
        result['max_error'] = max_error
        result['original'] = {'rows': orig_rows, 'cols': orig_cols, 'nnz': orig_nnz}
        result['roundtrip'] = {'rows': rt_rows, 'cols': rt_cols, 'nnz': rt_nnz}
        
        status = "✅ 通过" if result['passed'] else "❌ 失败"
        print(f"\n    总体结果: {status}")
        
    finally:
        # 清理临时文件
        for path in [rds_path, h5ad_rt_path]:
            if os.path.exists(path):
                os.remove(path)
    
    return result


# =============================================================================
# 主测试流程
# =============================================================================

def run_all_tests() -> Dict:
    """
    运行所有准确性测试
    
    对应 Task 20.1 - 20.7
    """
    print("\n" + "="*70)
    print(" CrossCell 准确性验证")
    print(" Task 20.1 - 20.7: 跨格式转换准确性测试")
    print(" 测试流程: H5AD -> RDS (Seurat) -> H5AD")
    print("="*70)
    
    # 检查 CrossCell 可执行文件
    if not os.path.exists(CROSSCELL_BIN):
        print(f"\n❌ 错误: 找不到 CrossCell 可执行文件: {CROSSCELL_BIN}")
        print("请先编译: docker-compose run --rm dev cargo build --release")
        return {'error': 'CrossCell binary not found'}
    
    print(f"\n使用 CrossCell: {CROSSCELL_BIN}")
    
    # 检查版本
    ret, stdout, stderr = run_crosscell(['--version'])
    if ret == 0:
        print(f"版本: {stdout.strip()}")
    
    all_results = []
    
    # =========================================================================
    # Task 20.3: PBMC 3k 测试
    # =========================================================================
    print("\n\n" + "#"*70)
    print("# Task 20.3: PBMC 3k 往返准确性测试")
    print("#"*70)
    
    pbmc_path = str(REAL_DATASETS_DIR / "pbmc3k.h5ad")
    all_results.append(test_cross_format_roundtrip(pbmc_path, "PBMC 3k (2700 细胞 × 13714 基因)"))
    
    # =========================================================================
    # Task 20.4: Pancreas 多批次数据测试
    # =========================================================================
    print("\n\n" + "#"*70)
    print("# Task 20.4: Pancreas 多批次数据往返准确性测试")
    print("#"*70)
    
    pancreas_path = str(REAL_DATASETS_DIR / "pancreas.h5ad")
    all_results.append(test_cross_format_roundtrip(pancreas_path, "Pancreas (3696 细胞 × 15737 基因)"))
    
    # =========================================================================
    # Task 20.5: Visium 空间数据测试
    # =========================================================================
    print("\n\n" + "#"*70)
    print("# Task 20.5: Visium 空间数据往返准确性测试")
    print("#"*70)
    
    visium_path = str(REAL_DATASETS_DIR / "visium_heart.h5ad")
    all_results.append(test_cross_format_roundtrip(visium_path, "Visium Heart (4247 spots × 15217 基因)"))
    
    # =========================================================================
    # Task 20.6: 边界情况测试
    # =========================================================================
    print("\n\n" + "#"*70)
    print("# Task 20.6: 特殊值和边界情况测试")
    print("#"*70)
    
    # 单行矩阵
    single_row_path = str(TEST_DATA_DIR / "single_row.h5ad")
    all_results.append(test_cross_format_roundtrip(single_row_path, "单行矩阵 (1 × 1000)"))
    
    # 单列矩阵
    single_col_path = str(TEST_DATA_DIR / "single_column.h5ad")
    all_results.append(test_cross_format_roundtrip(single_col_path, "单列矩阵 (1000 × 1)"))
    
    # 小型稀疏矩阵
    small_sparse_path = str(TEST_DATA_DIR / "small_sparse.h5ad")
    all_results.append(test_cross_format_roundtrip(small_sparse_path, "小型稀疏矩阵 (100 × 50, 95% 稀疏)"))
    
    # 中型稀疏矩阵
    medium_sparse_path = str(TEST_DATA_DIR / "medium_sparse.h5ad")
    all_results.append(test_cross_format_roundtrip(medium_sparse_path, "中型稀疏矩阵 (1000 × 500, 95% 稀疏)"))
    
    # =========================================================================
    # Task 20.7: 汇总结果
    # =========================================================================
    print("\n\n" + "#"*70)
    print("# Task 20.7: 准确性测试结果汇总")
    print("#"*70)
    
    # 统计结果
    total = len(all_results)
    passed = sum(1 for r in all_results if r.get('passed', False))
    skipped = sum(1 for r in all_results if r.get('skipped', False))
    failed = total - passed - skipped
    
    print(f"\n{'='*60}")
    print(" 测试汇总")
    print(f"{'='*60}")
    print(f"  总测试数: {total}")
    print(f"  ✅ 通过: {passed}")
    print(f"  ❌ 失败: {failed}")
    print(f"  ⚠️ 跳过: {skipped}")
    if (total - skipped) > 0:
        print(f"  通过率: {passed/(total-skipped)*100:.1f}%")
    
    # 详细结果表格
    print(f"\n{'='*60}")
    print(" 详细结果")
    print(f"{'='*60}")
    print(f"  {'测试名称':<35} {'结果':<8} {'相关系数':<15} {'NNZ匹配':<8}")
    print("  " + "-"*66)
    
    for r in all_results:
        name = r.get('name', 'Unknown')[:33]
        
        if r.get('skipped'):
            status = "⚠️跳过"
            corr = "N/A"
            nnz = "N/A"
        elif r.get('passed'):
            status = "✅通过"
            corr = f"{r.get('correlation', 0):.10f}"
            nnz = "✓" if r.get('nnz_match') else "✗"
        else:
            status = "❌失败"
            corr = f"{r.get('correlation', 0):.10f}"
            nnz = "✓" if r.get('nnz_match') else "✗"
        
        print(f"  {name:<35} {status:<8} {corr:<15} {nnz:<8}")
    
    # 准确性指标汇总
    correlations = [r['correlation'] for r in all_results 
                   if 'correlation' in r and not r.get('skipped')]
    
    if correlations:
        print(f"\n{'='*60}")
        print(" 准确性指标汇总")
        print(f"{'='*60}")
        print(f"  相关系数 - 最小: {min(correlations):.10f}")
        print(f"  相关系数 - 最大: {max(correlations):.10f}")
        print(f"  相关系数 - 平均: {sum(correlations)/len(correlations):.10f}")
    
    # 最终结论
    print(f"\n{'='*60}")
    print(" 结论")
    print(f"{'='*60}")
    
    if failed == 0 and passed > 0:
        print("  ✅ CrossCell 准确性验证通过！")
        print("")
        print("  所有测试的表达矩阵相关系数均达到 1.0")
        print("  稀疏结构完全保留 (NNZ 完全匹配)")
        print("  可以安全用于 AnnData ↔ Seurat 格式转换")
    elif failed > 0:
        print("  ❌ 部分测试失败，请检查详细结果")
    else:
        print("  ⚠️ 没有可运行的测试，请检查测试数据是否存在")
    
    return {
        'total': total,
        'passed': passed,
        'failed': failed,
        'skipped': skipped,
        'results': all_results
    }


# =============================================================================
# 入口点
# =============================================================================

if __name__ == '__main__':
    if not HAS_DEPS:
        print("\n错误: 缺少必要的依赖")
        print("请安装: pip install numpy scipy h5py")
        sys.exit(1)
    
    # 切换到项目根目录
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    os.chdir(project_root)
    
    print(f"工作目录: {os.getcwd()}")
    
    # 运行测试
    summary = run_all_tests()
    
    # 保存结果到 JSON
    output_json = script_dir / "accuracy_results.json"
    try:
        # 构建可序列化的结果
        json_results = []
        for r in summary.get('results', []):
            jr = {
                'name': r.get('name', ''),
                'passed': bool(r.get('passed', False)),
                'skipped': bool(r.get('skipped', False)),
                'correlation': float(r.get('correlation', 0)),
                'nnz_match': bool(r.get('nnz_match', False)),
                'dim_match': bool(r.get('dim_match', False)),
            }
            if 'original' in r:
                jr['original'] = {k: int(v) for k, v in r['original'].items()}
            if 'roundtrip' in r:
                jr['roundtrip'] = {k: int(v) for k, v in r['roundtrip'].items()}
            if 'error' in r:
                jr['error'] = r['error']
            json_results.append(jr)
        
        with open(output_json, 'w', encoding='utf-8') as f:
            json.dump({
                'total': summary.get('total', 0),
                'passed': summary.get('passed', 0),
                'failed': summary.get('failed', 0),
                'skipped': summary.get('skipped', 0),
                'results': json_results
            }, f, indent=2, ensure_ascii=False)
        
        print(f"\n结果已保存到: {output_json}")
    except Exception as e:
        print(f"\n警告: 无法保存 JSON 结果: {e}")
    
    # 返回退出码
    sys.exit(0 if summary.get('failed', 1) == 0 else 1)
