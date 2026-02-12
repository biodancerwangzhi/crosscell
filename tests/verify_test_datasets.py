#!/usr/bin/env python3
"""
验证真实测试数据集的完整性和质量

用途：Task 20.1 - 验证下载的数据集
"""

import sys
from pathlib import Path

try:
    import scanpy as sc
    import numpy as np
    import pandas as pd
except ImportError:
    print("错误: 需要安装 scanpy")
    sys.exit(1)

DATA_DIR = Path("tests/data/real_datasets")

def verify_pbmc3k():
    """验证 PBMC 3k 数据集"""
    print("\n" + "="*60)
    print("验证 PBMC 3k 数据集...")
    print("="*60)
    
    file_path = DATA_DIR / "pbmc3k.h5ad"
    if not file_path.exists():
        print(f"❌ 文件不存在: {file_path}")
        return False
    
    try:
        adata = sc.read_h5ad(file_path)
        
        # 验证基本结构
        checks = {
            "表达矩阵存在": adata.X is not None,
            "细胞数 > 2000": adata.n_obs > 2000,
            "基因数 > 10000": adata.n_vars > 10000,
            "稀疏度 > 90%": (1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) > 0.9,
            "包含元数据": len(adata.obs.columns) > 0,
            "包含 PCA": 'X_pca' in adata.obsm,
            "包含 UMAP": 'X_umap' in adata.obsm,
        }
        
        all_passed = True
        for check, result in checks.items():
            status = "✅" if result else "❌"
            print(f"  {status} {check}")
            if not result:
                all_passed = False
        
        # 打印详细信息
        print(f"\n  详细信息:")
        print(f"    细胞数: {adata.n_obs:,}")
        print(f"    基因数: {adata.n_vars:,}")
        print(f"    稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
        print(f"    元数据列: {list(adata.obs.columns)}")
        print(f"    降维: {list(adata.obsm.keys())}")
        
        return all_passed
        
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return False

def verify_pancreas():
    """验证 Pancreas 数据集"""
    print("\n" + "="*60)
    print("验证 Pancreas 数据集...")
    print("="*60)
    
    file_path = DATA_DIR / "pancreas.h5ad"
    if not file_path.exists():
        print(f"❌ 文件不存在: {file_path}")
        return False
    
    try:
        adata = sc.read_h5ad(file_path)
        
        # 验证基本结构
        checks = {
            "表达矩阵存在": adata.X is not None,
            "细胞数 > 2000": adata.n_obs > 2000,
            "基因数 > 10000": adata.n_vars > 10000,
            "包含 Categorical 列": any(adata.obs[col].dtype.name == 'category' for col in adata.obs.columns),
            "包含 PCA": 'X_pca' in adata.obsm,
            "包含 UMAP": 'X_umap' in adata.obsm,
        }
        
        all_passed = True
        for check, result in checks.items():
            status = "✅" if result else "❌"
            print(f"  {status} {check}")
            if not result:
                all_passed = False
        
        # 打印详细信息
        print(f"\n  详细信息:")
        print(f"    细胞数: {adata.n_obs:,}")
        print(f"    基因数: {adata.n_vars:,}")
        print(f"    稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
        print(f"    元数据列: {list(adata.obs.columns)}")
        
        # 检查 Categorical 列
        categorical_cols = [col for col in adata.obs.columns if adata.obs[col].dtype.name == 'category']
        if categorical_cols:
            print(f"    Categorical 列: {categorical_cols}")
            for col in categorical_cols[:2]:  # 只显示前两个
                print(f"      - {col}: {adata.obs[col].nunique()} 类别")
        
        return all_passed
        
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return False

def verify_visium():
    """验证 Visium 数据集"""
    print("\n" + "="*60)
    print("验证 Visium 数据集...")
    print("="*60)
    
    file_path = DATA_DIR / "visium_heart.h5ad"
    if not file_path.exists():
        print(f"❌ 文件不存在: {file_path}")
        return False
    
    try:
        adata = sc.read_h5ad(file_path)
        
        # 验证基本结构
        checks = {
            "表达矩阵存在": adata.X is not None,
            "Spots 数 > 4000": adata.n_obs > 4000,
            "基因数 > 10000": adata.n_vars > 10000,
            "包含空间坐标": 'spatial' in adata.obsm,
            "空间坐标是 2D": adata.obsm['spatial'].shape[1] == 2 if 'spatial' in adata.obsm else False,
            "包含空间数据": 'spatial' in adata.uns,
            "包含 PCA": 'X_pca' in adata.obsm,
            "包含 UMAP": 'X_umap' in adata.obsm,
        }
        
        all_passed = True
        for check, result in checks.items():
            status = "✅" if result else "❌"
            print(f"  {status} {check}")
            if not result:
                all_passed = False
        
        # 打印详细信息
        print(f"\n  详细信息:")
        print(f"    Spots 数: {adata.n_obs:,}")
        print(f"    基因数: {adata.n_vars:,}")
        print(f"    稀疏度: {(1 - adata.X.nnz / (adata.n_obs * adata.n_vars)) * 100:.2f}%")
        
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            print(f"    空间坐标范围:")
            print(f"      X: [{coords[:, 0].min():.2f}, {coords[:, 0].max():.2f}]")
            print(f"      Y: [{coords[:, 1].min():.2f}, {coords[:, 1].max():.2f}]")
        
        if 'spatial' in adata.uns:
            spatial_data = adata.uns['spatial']
            sample_keys = list(spatial_data.keys())
            if sample_keys:
                sample_key = sample_keys[0]
                print(f"    空间数据键: {sample_key}")
                if 'images' in spatial_data[sample_key]:
                    images = spatial_data[sample_key]['images']
                    print(f"    图像: {list(images.keys())}")
                    for img_key, img in images.items():
                        print(f"      - {img_key}: {img.shape}")
                if 'scalefactors' in spatial_data[sample_key]:
                    scalefactors = spatial_data[sample_key]['scalefactors']
                    print(f"    Scale factors: {list(scalefactors.keys())}")
        
        return all_passed
        
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return False

def main():
    """主函数"""
    print("="*60)
    print("CrossCell 真实测试数据集验证工具")
    print("Task 20.1 - 验证数据集完整性")
    print("="*60)
    
    results = {
        "PBMC 3k": verify_pbmc3k(),
        "Pancreas": verify_pancreas(),
        "Visium": verify_visium(),
    }
    
    # 汇总结果
    print("\n" + "="*60)
    print("验证结果汇总")
    print("="*60)
    
    for dataset, passed in results.items():
        status = "✅ 通过" if passed else "❌ 失败"
        print(f"  {dataset}: {status}")
    
    all_passed = all(results.values())
    
    if all_passed:
        print("\n" + "="*60)
        print("✅ 所有数据集验证通过！")
        print("="*60)
        print("\n数据集已准备就绪，可以进行准确性测试。")
        print("\n下一步:")
        print("  - Task 20.2: 实现准确性测试框架")
        print("  - Task 20.3: PBMC 3k 往返准确性测试")
        print("  - Task 20.4: Pancreas 多批次数据往返准确性测试")
        print("  - Task 20.5: Visium 空间数据往返准确性测试")
        return 0
    else:
        print("\n" + "="*60)
        print("⚠️  部分数据集验证失败")
        print("="*60)
        print("\n请检查失败的数据集并重新下载。")
        return 1

if __name__ == "__main__":
    sys.exit(main())
