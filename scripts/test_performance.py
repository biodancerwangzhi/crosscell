#!/usr/bin/env python3
"""
CrossCell 性能优化测试脚本

本脚本用于测试 CrossCell 的性能优化功能，包括：
1. Lazy Loading 和分块处理 (Task 21)
2. 并行稀疏矩阵转换
3. 内存使用估算
4. 流式处理大文件 (Task 30)
5. CLI 参数测试

使用方法：
    # 在 Docker 容器中运行
    docker-compose run --rm dev python3 scripts/test_performance.py
    
    # 或者直接运行（需要安装依赖）
    python3 scripts/test_performance.py

依赖：
    - anndata
    - numpy
    - scipy
    - pandas
"""

import os
import sys
import time
import subprocess
import tempfile
import json
from pathlib import Path

# 尝试导入依赖
try:
    import numpy as np
    import anndata as ad
    from scipy import sparse
    import pandas as pd
    HAS_DEPS = True
except ImportError as e:
    print(f"⚠️ 缺少依赖: {e}")
    print("请安装: pip install anndata numpy scipy pandas")
    HAS_DEPS = False


def print_header(title: str):
    """打印测试标题"""
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def print_result(name: str, passed: bool, details: str = ""):
    """打印测试结果"""
    status = "✅ 通过" if passed else "❌ 失败"
    print(f"  {status} - {name}")
    if details:
        print(f"         {details}")


def create_test_h5ad(n_cells: int, n_genes: int, sparsity: float, output_path: str):
    """
    创建测试用的 .h5ad 文件
    
    参数:
        n_cells: 细胞数量
        n_genes: 基因数量
        sparsity: 稀疏度 (0-1)
        output_path: 输出文件路径
    
    返回:
        文件大小（字节）
    """
    print(f"  创建测试数据: {n_cells} 细胞 × {n_genes} 基因, 稀疏度 {sparsity*100:.0f}%")
    
    # 生成稀疏矩阵
    nnz = int(n_cells * n_genes * (1 - sparsity))
    row_indices = np.random.randint(0, n_cells, nnz)
    col_indices = np.random.randint(0, n_genes, nnz)
    data = np.random.rand(nnz).astype(np.float32)
    
    X = sparse.csr_matrix((data, (row_indices, col_indices)), shape=(n_cells, n_genes))
    
    # 创建 AnnData 对象
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['A', 'B', 'C'], n_cells),
        'n_counts': np.random.randint(1000, 10000, n_cells),
    }, index=[f'cell_{i}' for i in range(n_cells)])
    
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_genes)],
    }, index=[f'gene_{i}' for i in range(n_genes)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # 添加降维结果
    if n_cells >= 50:
        adata.obsm['X_pca'] = np.random.rand(n_cells, min(50, n_cells)).astype(np.float32)
        adata.obsm['X_umap'] = np.random.rand(n_cells, 2).astype(np.float32)
    
    # 保存文件
    adata.write_h5ad(output_path)
    
    file_size = os.path.getsize(output_path)
    print(f"  文件大小: {file_size / 1024 / 1024:.2f} MB")
    
    return file_size


def test_cli_help():
    """测试 CLI 帮助信息"""
    print_header("测试 1: CLI 帮助信息")
    
    try:
        result = subprocess.run(
            ['./target/debug/crosscell', 'convert', '--help'],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        # 检查新增的参数是否存在
        help_text = result.stdout
        
        has_lazy = '--lazy' in help_text
        has_chunk_size = '--chunk-size' in help_text
        has_estimate_memory = '--estimate-memory' in help_text
        has_streaming = '--streaming' in help_text
        
        print_result("--lazy 参数", has_lazy)
        print_result("--chunk-size 参数", has_chunk_size)
        print_result("--estimate-memory 参数", has_estimate_memory)
        print_result("--streaming 参数", has_streaming)
        
        return has_lazy and has_chunk_size and has_estimate_memory and has_streaming
        
    except FileNotFoundError:
        print("  ⚠️ 未找到 crosscell 可执行文件，请先编译: cargo build")
        return False
    except Exception as e:
        print(f"  ❌ 错误: {e}")
        return False


def test_memory_estimation():
    """测试内存估算功能"""
    print_header("测试 2: 内存估算功能")
    
    if not HAS_DEPS:
        print("  ⚠️ 跳过（缺少依赖）")
        return True
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建测试文件
        h5ad_path = os.path.join(tmpdir, 'test.h5ad')
        create_test_h5ad(1000, 500, 0.95, h5ad_path)
        
        try:
            # 运行带 --estimate-memory 的转换
            result = subprocess.run(
                ['./target/debug/crosscell', 'convert',
                 '-i', h5ad_path,
                 '-o', os.path.join(tmpdir, 'test.rds'),
                 '-f', 'seurat',
                 '--estimate-memory'],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            output = result.stdout + result.stderr
            
            # 检查输出是否包含内存估算信息
            has_memory_info = 'Memory Estimation' in output or '内存' in output or 'MB' in output or 'GB' in output
            has_file_size = 'File size' in output or '文件大小' in output
            
            print_result("内存估算输出", has_memory_info, 
                        "包含内存估算信息" if has_memory_info else "未找到内存估算信息")
            print_result("文件大小显示", has_file_size)
            
            return has_memory_info
            
        except FileNotFoundError:
            print("  ⚠️ 未找到 crosscell 可执行文件")
            return False
        except Exception as e:
            print(f"  ❌ 错误: {e}")
            return False


def test_lazy_loading():
    """测试 Lazy Loading 功能"""
    print_header("测试 3: Lazy Loading 功能")
    
    if not HAS_DEPS:
        print("  ⚠️ 跳过（缺少依赖）")
        return True
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建测试文件
        h5ad_path = os.path.join(tmpdir, 'test.h5ad')
        create_test_h5ad(2000, 1000, 0.95, h5ad_path)
        
        rds_path = os.path.join(tmpdir, 'test.rds')
        
        try:
            # 测试普通模式
            print("\n  普通模式转换...")
            start_time = time.time()
            result_normal = subprocess.run(
                ['./target/debug/crosscell', 'convert',
                 '-i', h5ad_path,
                 '-o', rds_path,
                 '-f', 'seurat'],
                capture_output=True,
                text=True,
                timeout=120
            )
            normal_time = time.time() - start_time
            normal_success = result_normal.returncode == 0
            
            print_result("普通模式转换", normal_success, f"耗时: {normal_time:.2f}s")
            
            # 删除输出文件
            if os.path.exists(rds_path):
                os.remove(rds_path)
            
            # 测试 Lazy 模式
            print("\n  Lazy 模式转换...")
            start_time = time.time()
            result_lazy = subprocess.run(
                ['./target/debug/crosscell', 'convert',
                 '-i', h5ad_path,
                 '-o', rds_path,
                 '-f', 'seurat',
                 '--lazy',
                 '--chunk-size', '500'],
                capture_output=True,
                text=True,
                timeout=120
            )
            lazy_time = time.time() - start_time
            lazy_success = result_lazy.returncode == 0
            
            print_result("Lazy 模式转换", lazy_success, f"耗时: {lazy_time:.2f}s")
            
            # 比较结果
            if normal_success and lazy_success:
                print(f"\n  📊 性能对比:")
                print(f"     普通模式: {normal_time:.2f}s")
                print(f"     Lazy 模式: {lazy_time:.2f}s")
                print(f"     时间差异: {((lazy_time - normal_time) / normal_time * 100):+.1f}%")
            
            return normal_success and lazy_success
            
        except FileNotFoundError:
            print("  ⚠️ 未找到 crosscell 可执行文件")
            return False
        except Exception as e:
            print(f"  ❌ 错误: {e}")
            return False


def test_rust_unit_tests():
    """运行 Rust 单元测试"""
    print_header("测试 4: Rust 单元测试")
    
    tests = [
        ("LazyMatrix 测试", "ir::expression"),
        ("ChunkedReader 测试", "backend::tests"),
        ("稀疏矩阵转换测试", "sparse::convert"),
        ("内存估算测试", "sparse::memory"),
    ]
    
    all_passed = True
    
    for name, test_filter in tests:
        try:
            result = subprocess.run(
                ['cargo', 'test', test_filter, '--', '--nocapture'],
                capture_output=True,
                text=True,
                timeout=120
            )
            
            # 解析测试结果
            output = result.stdout + result.stderr
            passed = 'test result: ok' in output or result.returncode == 0
            
            # 提取测试数量
            import re
            match = re.search(r'(\d+) passed', output)
            count = match.group(1) if match else "?"
            
            print_result(name, passed, f"{count} 个测试通过")
            
            if not passed:
                all_passed = False
                
        except Exception as e:
            print_result(name, False, str(e))
            all_passed = False
    
    return all_passed


def test_benchmark_exists():
    """检查基准测试是否存在"""
    print_header("测试 5: 基准测试套件")
    
    bench_file = Path("benches/sparse_conversion.rs")
    
    if not bench_file.exists():
        print_result("基准测试文件", False, "benches/sparse_conversion.rs 不存在")
        return False
    
    # 读取文件内容检查测试函数
    content = bench_file.read_text()
    
    benchmarks = [
        ("CSR → CSC 转换", "bench_csr_to_csc_conversion"),
        ("并行 vs 串行", "bench_parallel_vs_serial"),
        ("行子集操作", "bench_subset_rows"),
        ("列子集操作", "bench_subset_contiguous_columns"),
        ("内存估算", "bench_memory_estimation"),
        ("往返转换", "bench_roundtrip_conversion"),
    ]
    
    all_found = True
    for name, func_name in benchmarks:
        found = func_name in content
        print_result(name, found)
        if not found:
            all_found = False
    
    return all_found


def test_performance_doc():
    """检查性能文档"""
    print_header("测试 6: 性能文档")
    
    doc_file = Path("docs/PERFORMANCE.md")
    
    if not doc_file.exists():
        print_result("性能文档", False, "docs/PERFORMANCE.md 不存在")
        return False
    
    content = doc_file.read_text()
    
    sections = [
        ("Lazy Loading 文档", "LazyMatrix"),
        ("分块处理文档", "ChunkedReader"),
        ("并行转换文档", "parallel"),
        ("CLI 参数文档", "--lazy"),
        ("内存优化文档", "Memory"),
        ("流式处理文档", "streaming"),
    ]
    
    all_found = True
    for name, keyword in sections:
        found = keyword.lower() in content.lower()
        print_result(name, found)
        if not found:
            all_found = False
    
    return all_found


def test_streaming_mode():
    """测试流式处理功能 (Task 30)"""
    print_header("测试 7: 流式处理功能")
    
    if not HAS_DEPS:
        print("  ⚠️ 跳过（缺少依赖）")
        return True
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建测试文件
        h5ad_path = os.path.join(tmpdir, 'test.h5ad')
        create_test_h5ad(3000, 1500, 0.95, h5ad_path)
        
        rds_path = os.path.join(tmpdir, 'test.rds')
        h5ad_out_path = os.path.join(tmpdir, 'test_out.h5ad')
        
        try:
            # 测试 H5AD -> RDS 流式转换
            print("\n  流式模式 H5AD -> RDS...")
            start_time = time.time()
            result_streaming = subprocess.run(
                ['./target/debug/crosscell', 'convert',
                 '-i', h5ad_path,
                 '-o', rds_path,
                 '-f', 'seurat',
                 '--streaming'],
                capture_output=True,
                text=True,
                timeout=180
            )
            streaming_time = time.time() - start_time
            streaming_success = result_streaming.returncode == 0
            
            print_result("H5AD -> RDS 流式转换", streaming_success, f"耗时: {streaming_time:.2f}s")
            
            if not streaming_success:
                print(f"    错误: {result_streaming.stderr[:200]}")
            
            # 测试 RDS -> H5AD 流式转换
            if streaming_success and os.path.exists(rds_path):
                print("\n  流式模式 RDS -> H5AD...")
                start_time = time.time()
                result_reverse = subprocess.run(
                    ['./target/debug/crosscell', 'convert',
                     '-i', rds_path,
                     '-o', h5ad_out_path,
                     '-f', 'anndata',
                     '--streaming'],
                    capture_output=True,
                    text=True,
                    timeout=180
                )
                reverse_time = time.time() - start_time
                reverse_success = result_reverse.returncode == 0
                
                print_result("RDS -> H5AD 流式转换", reverse_success, f"耗时: {reverse_time:.2f}s")
                
                if not reverse_success:
                    print(f"    错误: {result_reverse.stderr[:200]}")
                
                return streaming_success and reverse_success
            
            return streaming_success
            
        except FileNotFoundError:
            print("  ⚠️ 未找到 crosscell 可执行文件")
            return False
        except Exception as e:
            print(f"  ❌ 错误: {e}")
            return False


def test_streaming_benchmark_exists():
    """检查流式处理基准测试是否存在"""
    print_header("测试 8: 流式处理基准测试套件")
    
    bench_file = Path("benches/streaming_conversion.rs")
    
    if not bench_file.exists():
        print_result("基准测试文件", False, "benches/streaming_conversion.rs 不存在")
        return False
    
    # 读取文件内容检查测试函数
    content = bench_file.read_text()
    
    benchmarks = [
        ("流式读取基准", "streaming"),
        ("分块处理基准", "chunk"),
        ("Criterion 框架", "criterion"),
    ]
    
    all_found = True
    for name, keyword in benchmarks:
        found = keyword.lower() in content.lower()
        print_result(name, found)
        if not found:
            all_found = False
    
    return all_found


def main():
    """主函数"""
    print("\n" + "🚀" * 20)
    print("  CrossCell 性能优化测试")
    print("  Task 21 (Lazy Loading) + Task 30 (Streaming)")
    print("🚀" * 20)
    
    results = {}
    
    # Task 21 测试
    print("\n" + "=" * 60)
    print("  Part 1: Task 21 - Lazy Loading 和内存优化")
    print("=" * 60)
    
    results['CLI 帮助'] = test_cli_help()
    results['内存估算'] = test_memory_estimation()
    results['Lazy Loading'] = test_lazy_loading()
    results['Rust 单元测试'] = test_rust_unit_tests()
    results['基准测试套件'] = test_benchmark_exists()
    results['性能文档'] = test_performance_doc()
    
    # Task 30 测试
    print("\n" + "=" * 60)
    print("  Part 2: Task 30 - 流式处理大文件")
    print("=" * 60)
    
    results['流式处理'] = test_streaming_mode()
    results['流式基准测试'] = test_streaming_benchmark_exists()
    
    # 打印总结
    print_header("测试总结")
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, result in results.items():
        status = "✅" if result else "❌"
        print(f"  {status} {name}")
    
    print(f"\n  总计: {passed}/{total} 通过")
    
    if passed == total:
        print("\n  🎉 所有测试通过！性能优化功能正常。")
        return 0
    else:
        print(f"\n  ⚠️ {total - passed} 个测试失败，请检查。")
        return 1


if __name__ == "__main__":
    sys.exit(main())
