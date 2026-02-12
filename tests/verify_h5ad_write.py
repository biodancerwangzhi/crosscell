#!/usr/bin/env python3
"""
验证 Rust 写入的 .h5ad 文件是否可以被 Python anndata 正确读取
"""

import anndata as ad
import numpy as np
import sys
import tempfile
import os

def verify_metadata_file():
    """验证包含元数据的 .h5ad 文件"""
    print("🔍 验证元数据写入...")
    
    # 这个文件应该由 Rust 测试生成
    # 我们需要手动运行测试并保存文件
    test_file = "/tmp/test_metadata.h5ad"
    
    if not os.path.exists(test_file):
        print(f"⚠️  测试文件不存在: {test_file}")
        print("   请先运行 Rust 测试生成文件")
        return False
    
    try:
        # 读取文件
        adata = ad.read_h5ad(test_file)
        
        print(f"✓ 成功读取文件")
        print(f"  形状: {adata.shape}")
        print(f"  细胞元数据列: {list(adata.obs.columns)}")
        print(f"  基因元数据列: {list(adata.var.columns)}")
        
        # 验证细胞元数据
        assert adata.n_obs == 3, f"Expected 3 cells, got {adata.n_obs}"
        assert adata.n_vars == 2, f"Expected 2 genes, got {adata.n_vars}"
        
        # 验证列存在
        expected_obs_cols = ['cell_id', 'cell_type', 'n_genes']
        for col in expected_obs_cols:
            assert col in adata.obs.columns, f"Missing column: {col}"
        
        # 验证 Categorical 列
        assert adata.obs['cell_type'].dtype.name == 'category', "cell_type should be categorical"
        print(f"  cell_type categories: {list(adata.obs['cell_type'].cat.categories)}")
        
        # 验证数值
        assert np.allclose(adata.obs['n_genes'].values, [1.5, 2.5, 3.5]), "n_genes values mismatch"
        
        print("✅ 元数据验证通过！")
        return True
        
    except Exception as e:
        print(f"❌ 验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def create_test_file_for_rust():
    """创建一个测试文件供 Rust 读取"""
    print("\n📝 创建测试文件供 Rust 读取...")
    
    # 创建简单的 AnnData 对象
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    obs = {
        'n_genes': [1.5, 2.5, 3.5],
        'cell_id': ['A', 'B', 'C'],
        'cell_type': pd.Categorical(['type1', 'type2', 'type1'])
    }
    var = {
        'gene_name': ['gene1', 'gene2']
    }
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    output_file = "tests/data/test_metadata_write.h5ad"
    adata.write_h5ad(output_file)
    print(f"✓ 已创建: {output_file}")
    print(f"  形状: {adata.shape}")
    print(f"  obs 列: {list(adata.obs.columns)}")
    print(f"  var 列: {list(adata.var.columns)}")

if __name__ == "__main__":
    import pandas as pd
    
    print("=" * 60)
    print("验证 CrossCell .h5ad 元数据写入")
    print("=" * 60)
    
    # 创建测试文件
    create_test_file_for_rust()
    
    print("\n" + "=" * 60)
    print("提示：运行以下命令生成 Rust 写入的文件进行验证：")
    print("  docker run -it --rm -v ${PWD}:/workspace crosscell-dev bash")
    print("  cd /workspace")
    print("  cargo test test_write_metadata -- --nocapture")
    print("=" * 60)
