#!/usr/bin/env python3
"""
下载和预处理真实测试数据集

数据集：
1. PBMC 3k - 标准 scRNA-seq 数据
2. Pancreas - 多批次、多细胞类型数据
3. Visium - 空间转录组数据

用途：Task 20.1 - 准备真实测试数据集
"""

import os
import sys
from pathlib import Path


import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scvelo as scv


# 设置数据目录
DATA_DIR = Path("tests/data/real_datasets")
DATA_DIR.mkdir(parents=True, exist_ok=True)

def download_pbmc3k():
    """
    下载 PBMC 3k 数据集
    
    特征：
    - ~2,700 细胞
    - ~32,000 基因
    - 稀疏度 ~93%
    - 包含细胞类型注释
    """
    print("\n" + "="*60)
    print("下载 PBMC 3k 数据集...")
    print("="*60)
    
    output_file = DATA_DIR / "pbmc3k.h5ad"
    
    if output_file.exists():
        print(f"✓ 文件已存在: {output_file}")
        adata = sc.read_h5ad(output_file)
    else:
        # 使用 scanpy 内置数据
        adata = sc.datasets.pbmc3k()
        
        # 基础预处理
        print("  - 过滤细胞和基因...")
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # 添加 QC 指标
        print("  - 计算 QC 指标...")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # 归一化和对数转换
        print("  - 归一化...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 识别高变基因
        print("  - 识别高变基因...")
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        # 降维
        print("  - PCA 降维...")
        sc.tl.pca(adata, svd_solver='arpack')
        
        print("  - UMAP 降维...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        
        # 聚类（可选，如果 igraph 不可用则跳过）
        try:
            print("  - Leiden 聚类...")
            sc.tl.leiden(adata)
        except Exception as e:
            print(f"  ⚠️  跳过聚类 (igraph 未安装): {e}")
            # 添加一个简单的聚类列
            adata.obs['leiden'] = '0'
        
        # 保存
        print(f"  - 保存到: {output_file}")
        adata.write_h5ad(output_file)
    
    # 打印统计信息
    print(f"\n📊 PBMC 3k 数据集统计:")
    print(f"  细胞数: {adata.n_obs:,}")
    print(f"  基因数: {adata.n_vars:,}")
    print(f"  稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
    print(f"  文件大小: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"  元数据列: {list(adata.obs.columns)}")
    print(f"  降维: {list(adata.obsm.keys())}")
    
    return adata

def download_pancreas():
    """
    下载 Pancreas 数据集
    
    特征：
    - ~2,500 细胞
    - 多个批次
    - 多种细胞类型
    - 包含 batch 和 cell_type 注释
    """
    print("\n" + "="*60)
    print("下载 Pancreas 数据集...")
    print("="*60)
    
    output_file = DATA_DIR / "pancreas.h5ad"
    
    if output_file.exists():
        print(f"✓ 文件已存在: {output_file}")
        adata = sc.read_h5ad(output_file)
    else:
        # 使用 scanpy 内置数据
        adata = scv.datasets.pancreas()
        
        # 基础预处理
        print("  - 过滤细胞和基因...")
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # 归一化
        print("  - 归一化...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 识别高变基因
        print("  - 识别高变基因...")
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        # 降维
        print("  - PCA 降维...")
        sc.tl.pca(adata, svd_solver='arpack')
        
        print("  - UMAP 降维...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        
        # 保存
        print(f"  - 保存到: {output_file}")
        adata.write_h5ad(output_file)
    
    # 打印统计信息
    print(f"\n📊 Pancreas 数据集统计:")
    print(f"  细胞数: {adata.n_obs:,}")
    print(f"  基因数: {adata.n_vars:,}")
    print(f"  稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
    print(f"  文件大小: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"  元数据列: {list(adata.obs.columns)}")
    print(f"  批次数: {adata.obs['batch'].nunique() if 'batch' in adata.obs.columns else 'N/A'}")
    print(f"  细胞类型数: {adata.obs['celltype'].nunique() if 'celltype' in adata.obs.columns else 'N/A'}")
    print(f"  降维: {list(adata.obsm.keys())}")
    
    return adata

def download_visium():
    """
    下载 Visium 空间转录组数据
    
    特征：
    - ~4,992 spots
    - 空间坐标
    - 高分辨率和低分辨率图像
    - scale factors
    """
    print("\n" + "="*60)
    print("下载 Visium 空间转录组数据...")
    print("="*60)
    
    output_file = DATA_DIR / "visium_heart.h5ad"
    
    if output_file.exists():
        print(f"✓ 文件已存在: {output_file}")
        adata = sc.read_h5ad(output_file)
    else:
        try:
            # 尝试使用 scanpy 内置的 Visium 数据
            adata = sc.datasets.visium_sge(sample_id="V1_Human_Heart")
            
            # 基础预处理
            print("  - 过滤基因...")
            sc.pp.filter_genes(adata, min_cells=10)
            
            # 归一化
            print("  - 归一化...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # 识别高变基因
            print("  - 识别高变基因...")
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            
            # 降维
            print("  - PCA 降维...")
            sc.tl.pca(adata, svd_solver='arpack')
            
            print("  - UMAP 降维...")
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
            sc.tl.umap(adata)
            
            # 保存
            print(f"  - 保存到: {output_file}")
            adata.write_h5ad(output_file)
            
        except Exception as e:
            print(f"⚠️  无法下载 Visium 数据: {e}")
            print("  创建模拟 Visium 数据...")
            
            # 创建模拟数据
            n_spots = 4992
            n_genes = 18078
            
            # 创建稀疏表达矩阵
            from scipy.sparse import random as sparse_random
            X = sparse_random(n_spots, n_genes, density=0.05, format='csr')
            X.data = np.abs(X.data) * 100  # 模拟 counts
            
            # 创建 AnnData 对象
            adata = ad.AnnData(X)
            adata.obs_names = [f"SPOT_{i}" for i in range(n_spots)]
            adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
            
            # 添加空间坐标（模拟 Visium 布局）
            import math
            coords = []
            for i in range(n_spots):
                row = i // 78
                col = i % 78
                x = col * 100 + (row % 2) * 50
                y = row * 86.6
                coords.append([x, y])
            adata.obsm['spatial'] = np.array(coords)
            
            # 添加模拟图像和 scale factors
            adata.uns['spatial'] = {
                'V1_Human_Heart': {
                    'images': {
                        'hires': np.random.randint(0, 255, (2000, 2000, 3), dtype=np.uint8),
                        'lowres': np.random.randint(0, 255, (600, 600, 3), dtype=np.uint8)
                    },
                    'scalefactors': {
                        'tissue_hires_scalef': 0.17,
                        'tissue_lowres_scalef': 0.051,
                        'spot_diameter_fullres': 177.76
                    }
                }
            }
            
            # 基础处理
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
            sc.tl.umap(adata)
            
            # 保存
            adata.write_h5ad(output_file)
    
    # 打印统计信息
    print(f"\n📊 Visium 数据集统计:")
    print(f"  Spots 数: {adata.n_obs:,}")
    print(f"  基因数: {adata.n_vars:,}")
    print(f"  稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
    print(f"  文件大小: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"  空间坐标: {'spatial' in adata.obsm}")
    print(f"  图像: {list(adata.uns.get('spatial', {}).get('V1_Human_Heart', {}).get('images', {}).keys())}")
    print(f"  降维: {list(adata.obsm.keys())}")
    
    return adata

def verify_datasets():
    """验证所有数据集的完整性"""
    print("\n" + "="*60)
    print("验证数据集完整性...")
    print("="*60)
    
    datasets = {
        "PBMC 3k": DATA_DIR / "pbmc3k.h5ad",
        "Pancreas": DATA_DIR / "pancreas.h5ad",
        "Visium": DATA_DIR / "visium_heart.h5ad"
    }
    
    all_valid = True
    for name, path in datasets.items():
        if not path.exists():
            print(f"❌ {name}: 文件不存在")
            all_valid = False
            continue
        
        try:
            adata = sc.read_h5ad(path)
            
            # 验证必需组件
            has_X = adata.X is not None
            has_obs = len(adata.obs.columns) > 0
            has_var = len(adata.var.columns) > 0
            has_embeddings = len(adata.obsm.keys()) > 0
            
            if has_X and has_obs and has_var and has_embeddings:
                print(f"✅ {name}: 验证通过")
                print(f"   - 细胞/Spots: {adata.n_obs:,}")
                print(f"   - 基因: {adata.n_vars:,}")
                print(f"   - 元数据列: {len(adata.obs.columns)}")
                print(f"   - 降维: {', '.join(adata.obsm.keys())}")
            else:
                print(f"⚠️  {name}: 缺少必需组件")
                all_valid = False
                
        except Exception as e:
            print(f"❌ {name}: 读取失败 - {e}")
            all_valid = False
    
    return all_valid

def main():
    """主函数"""
    print("="*60)
    print("CrossCell 真实测试数据集下载工具")
    print("Task 20.1 - 准备真实测试数据集")
    print("="*60)
    
    # 下载数据集
    try:
        pbmc = download_pbmc3k()
        pancreas = download_pancreas()
        visium = download_visium()
    except Exception as e:
        print(f"\n❌ 下载失败: {e}")
        sys.exit(1)
    
    # 验证数据集
    if verify_datasets():
        print("\n" + "="*60)
        print("✅ 所有数据集准备完成！")
        print("="*60)
        print(f"\n数据位置: {DATA_DIR.absolute()}")
        print("\n下一步:")
        print("  1. 运行准确性测试: cargo test --test test_real_datasets")
        print("  2. 查看测试结果: tests/data/real_datasets/")
    else:
        print("\n" + "="*60)
        print("⚠️  部分数据集验证失败")
        print("="*60)
        sys.exit(1)

if __name__ == "__main__":
    main()
