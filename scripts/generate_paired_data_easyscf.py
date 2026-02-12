#!/usr/bin/env python3
"""
使用 easySCF 生成配对的 h5ad 和 rds 数据
easySCF 是一个成熟的转换工具，支持 Seurat V5，可以生成高质量的配对数据

工作流程：
1. Python: *.h5ad → *.h5 (easySCF 统一格式)
2. R: *.h5 → *_easyscf.rds (SimplifiedSeurat 格式，CrossCell 兼容)

用途：Task 20.1.5 - 生成金标准配对数据用于验证 CrossCell

输入：tests/data/real_datasets/*.h5ad
输出：tests/data/real_datasets/*_easyscf.rds
"""

import sys
import subprocess
from pathlib import Path

DATA_DIR = Path("tests/data/real_datasets")

# 数据集定义
DATASETS = [
    {
        "name": "PBMC 3k",
        "h5ad": DATA_DIR / "pbmc3k.h5ad",
        "h5": DATA_DIR / "pbmc3k_easyscf.h5",
        "rds": DATA_DIR / "pbmc3k_easyscf.rds",
        "is_spatial": False,
    },
    {
        "name": "Pancreas",
        "h5ad": DATA_DIR / "pancreas.h5ad",
        "h5": DATA_DIR / "pancreas_easyscf.h5",
        "rds": DATA_DIR / "pancreas_easyscf.rds",
        "is_spatial": False,
    },
    {
        "name": "Visium Heart",
        "h5ad": DATA_DIR / "visium_heart.h5ad",
        "h5": DATA_DIR / "visium_easyscf.h5",
        "rds": DATA_DIR / "visium_easyscf.rds",
        "is_spatial": True,
    },
]


def check_dependencies():
    """检查依赖"""
    print("检查依赖...")
    
    try:
        import scanpy as sc
        print(f"  ✓ scanpy {sc.__version__}")
    except ImportError:
        print("  ❌ 缺少 scanpy: pip install scanpy")
        return False
    
    try:
        from easySCFpy import saveH5, loadH5
        print("  ✓ easySCFpy")
    except ImportError:
        print("  ❌ 缺少 easySCFpy: pip install easySCFpy")
        return False
    
    # 检查 R 和 easySCFr
    try:
        result = subprocess.run(
            ["Rscript", "-e", "library(easySCFr); cat('OK')"],
            capture_output=True,
            text=True,
            timeout=30
        )
        if "OK" in result.stdout:
            print("  ✓ easySCFr (R)")
        else:
            print("  ❌ 缺少 easySCFr: install.packages('easySCFr')")
            print(f"     错误: {result.stderr}")
            return False
    except Exception as e:
        print(f"  ❌ R 检查失败: {e}")
        return False
    
    print()
    return True


def preprocess_for_easyscf(adata, name: str):
    """预处理数据以兼容 easySCF"""
    import numpy as np
    
    # 1. 确保 var_names 唯一
    if not adata.var_names.is_unique:
        print(f"  ⚠️ 发现重复基因名，正在去重...")
        adata.var_names_make_unique()
    
    # 2. 将 obs 中的 category 列转换为字符串
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == 'category':
            adata.obs[col] = adata.obs[col].astype(str)
            print(f"  ⚠️ 将 obs.{col} 从 category 转换为 string")
    
    # 3. 将 var 中的 category 列转换为字符串
    for col in adata.var.columns:
        if adata.var[col].dtype.name == 'category':
            adata.var[col] = adata.var[col].astype(str)
            print(f"  ⚠️ 将 var.{col} 从 category 转换为 string")
    
    # 4. 处理 var 中的布尔列
    for col in adata.var.columns:
        if adata.var[col].dtype == bool:
            adata.var[col] = adata.var[col].astype(str)
            print(f"  ⚠️ 将 var.{col} 从 bool 转换为 string")
    
    return adata


def convert_h5ad_to_h5(dataset: dict) -> bool:
    """步骤 1: h5ad → h5 (easySCF 统一格式)"""
    import scanpy as sc
    from easySCFpy import saveH5
    
    name = dataset["name"]
    h5ad_path = dataset["h5ad"]
    h5_path = dataset["h5"]
    is_spatial = dataset["is_spatial"]
    
    print(f"\n[1/2] {name}: h5ad → h5")
    
    if not h5ad_path.exists():
        print(f"  ❌ 找不到文件: {h5ad_path}")
        return False
    
    try:
        print(f"  读取: {h5ad_path}")
        adata = sc.read_h5ad(str(h5ad_path))
        
        # 预处理
        adata = preprocess_for_easyscf(adata, name)
        
        print(f"  数据: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
        print(f"  元数据: {list(adata.obs.columns)[:5]}...")
        print(f"  降维: {list(adata.obsm.keys())}")
        
        # 检测空间数据的 image_name
        image_name = "slice"
        if is_spatial and "spatial" in adata.uns:
            spatial_keys = list(adata.uns["spatial"].keys())
            if spatial_keys:
                image_name = spatial_keys[0]
                print(f"  检测到空间数据，image_name: {image_name}")
        
        print(f"  保存: {h5_path}")
        if is_spatial:
            saveH5(adata, str(h5_path), image_name=image_name)
        else:
            saveH5(adata, str(h5_path))
        
        if h5_path.exists():
            size_mb = h5_path.stat().st_size / 1024 / 1024
            print(f"  ✓ 成功 ({size_mb:.2f} MB)")
            return True
        else:
            print("  ❌ 文件未生成")
            return False
            
    except Exception as e:
        print(f"  ❌ 失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def convert_h5_to_rds(dataset: dict) -> bool:
    """步骤 2: h5 → rds (SimplifiedSeurat 格式，CrossCell 兼容)"""
    name = dataset["name"]
    h5_path = dataset["h5"]
    rds_path = dataset["rds"]
    is_spatial = dataset["is_spatial"]
    
    print(f"\n[2/2] {name}: h5 → rds")
    
    if not h5_path.exists():
        print(f"  ❌ 找不到文件: {h5_path}")
        return False
    
    # 检测空间数据的 image_name
    image_name = "slice"
    if is_spatial:
        import h5py
        try:
            with h5py.File(str(h5_path), 'r') as f:
                if 'uns' in f and 'spatial' in f['uns']:
                    spatial_keys = list(f['uns']['spatial'].keys())
                    if spatial_keys:
                        image_name = spatial_keys[0]
        except:
            pass
    
    # R 脚本文件
    r_script_file = DATA_DIR / "temp_convert.R"
    
    # R 脚本内容 - 转换为 SimplifiedSeurat 格式
    r_script = f'''
library(easySCFr)
library(Seurat)

cat("  读取: {h5_path}\\n")
seurat <- tryCatch({{
    loadH5("{h5_path}", image_name = "{image_name}")
}}, error = function(e) {{
    cat("  ⚠️ 加载失败，尝试不加载图像...\\n")
    loadH5("{h5_path}", image_name = NULL)
}})

cat(sprintf("  数据: %d 细胞 x %d 基因\\n", ncol(seurat), nrow(seurat)))
cat(sprintf("  Assays: %s\\n", paste(names(seurat@assays), collapse=", ")))

if (length(seurat@reductions) > 0) {{
    cat(sprintf("  Reductions: %s\\n", paste(names(seurat@reductions), collapse=", ")))
}}

# 转换为 SimplifiedSeurat 格式（CrossCell 兼容）
cat("  转换为 SimplifiedSeurat 格式...\\n")

# 获取 counts 矩阵
counts <- GetAssayData(seurat, layer = 'counts')

# 创建简化的 assay 结构（匹配 CrossCell 期望的格式）
simplified_assay <- list(
    class = "Assay5",
    layers = list(
        counts = counts,
        data = NULL
    ),
    features = rownames(counts),
    cells = colnames(counts),
    meta_features = data.frame(row.names = rownames(counts))
)

# 创建简化的 Seurat 结构
simplified <- list(
    class = "SimplifiedSeurat",
    project_name = "CrossCell",
    active_assay = "RNA",
    active_ident = factor(rep("0", ncol(counts))),
    assays = list(RNA = simplified_assay),
    meta_data = seurat@meta.data
)

# 添加降维结果
if (length(seurat@reductions) > 0) {{
    simplified$reductions <- list()
    for (red_name in names(seurat@reductions)) {{
        red <- seurat@reductions[[red_name]]
        simplified$reductions[[red_name]] <- list(
            cell_embeddings = Embeddings(red),
            feature_loadings = tryCatch(Loadings(red), error = function(e) matrix(nrow=0, ncol=0)),
            key = Key(red)
        )
    }}
}}

class(simplified) <- "SimplifiedSeurat"

cat("  保存: {rds_path}\\n")
saveRDS(simplified, file = "{rds_path}")

if (file.exists("{rds_path}")) {{
    size_mb <- file.info("{rds_path}")$size / 1024 / 1024
    cat(sprintf("  成功 (%.2f MB)\\n", size_mb))
}} else {{
    cat("  文件未生成\\n")
    quit(status = 1)
}}
'''
    
    try:
        # 写入临时 R 脚本文件
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        result = subprocess.run(
            ["Rscript", str(r_script_file)],
            capture_output=True,
            text=True,
            timeout=300
        )
        
        print(result.stdout)
        if result.returncode != 0:
            print(f"  R 错误: {result.stderr}")
            return False
        
        return rds_path.exists()
        
    except subprocess.TimeoutExpired:
        print("  ❌ 超时")
        return False
    except Exception as e:
        print(f"  ❌ 失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # 清理临时文件
        if r_script_file.exists():
            r_script_file.unlink()


def cleanup_temp_files(dataset: dict):
    """清理临时文件"""
    h5_path = dataset["h5"]
    if h5_path.exists():
        h5_path.unlink()
        print(f"  清理: {h5_path.name}")


def fallback_to_crosscell(dataset: dict) -> bool:
    """使用 CrossCell 作为回退方案"""
    h5ad_path = dataset["h5ad"]
    rds_path = dataset["rds"]
    
    print(f"\n[CrossCell] h5ad → rds")
    print(f"  输入: {h5ad_path}")
    print(f"  输出: {rds_path}")
    
    cmd = [
        "cargo", "run", "--release", "--bin", "crosscell", "--",
        "convert",
        "-i", str(h5ad_path),
        "-o", str(rds_path),
        "-f", "seurat"
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300
        )
        
        if result.returncode != 0:
            print(f"  ❌ CrossCell 转换失败:")
            print(result.stderr)
            return False
        
        if rds_path.exists():
            size_mb = rds_path.stat().st_size / 1024 / 1024
            print(f"  ✓ CrossCell 转换成功 ({size_mb:.2f} MB)")
            print(f"  ⚠️ 注意: 此数据由 CrossCell 生成，不是独立金标准")
            return True
        else:
            print("  ❌ RDS 文件未生成")
            return False
        
    except subprocess.TimeoutExpired:
        print("  ❌ 转换超时")
        return False
    except Exception as e:
        print(f"  ❌ 转换失败: {e}")
        return False


def convert_dataset(dataset: dict) -> bool:
    """转换单个数据集"""
    name = dataset["name"]
    print(f"\n{'='*60}")
    print(f"转换 {name}")
    print(f"{'='*60}")
    
    # 步骤 1: h5ad → h5
    if not convert_h5ad_to_h5(dataset):
        return False
    
    # 步骤 2: h5 → rds
    success = convert_h5_to_rds(dataset)
    
    # 清理临时文件
    cleanup_temp_files(dataset)
    
    if success:
        print(f"\n✅ {name} 转换成功！(easySCF 金标准)")
    else:
        # 对于空间数据，尝试使用 CrossCell 作为回退
        if dataset["is_spatial"]:
            print(f"\n⚠️ easySCF 不支持 Seurat V5 空间数据格式")
            print(f"   尝试使用 CrossCell 作为回退...")
            success = fallback_to_crosscell(dataset)
            if success:
                print(f"\n✅ {name} 转换成功！(CrossCell 回退)")
            else:
                print(f"\n❌ {name} 转换失败")
        else:
            print(f"\n❌ {name} 转换失败")
    
    return success


def main():
    print("\n" + "="*60)
    print("easySCF 配对数据生成脚本")
    print("生成金标准配对数据用于验证 CrossCell")
    print("="*60 + "\n")
    
    # 检查依赖
    if not check_dependencies():
        print("\n❌ 依赖检查失败，请安装缺少的包")
        sys.exit(1)
    
    # 转换所有数据集
    results = {}
    for dataset in DATASETS:
        results[dataset["name"]] = convert_dataset(dataset)
    
    # 汇总
    print("\n" + "="*60)
    print("转换汇总")
    print("="*60 + "\n")
    
    success_count = sum(results.values())
    total_count = len(results)
    
    for name, success in results.items():
        status = "✅" if success else "❌"
        print(f"  {status} {name}")
    
    print(f"\n总计: {success_count}/{total_count} 成功")
    
    if success_count == total_count:
        print("\n🎉 所有数据集转换成功！")
        print("\n生成的配对数据:")
        for dataset in DATASETS:
            rds_path = dataset["rds"]
            if rds_path.exists():
                size_mb = rds_path.stat().st_size / 1024 / 1024
                print(f"  - {rds_path} ({size_mb:.2f} MB)")
    else:
        print("\n⚠️ 部分数据集转换失败")
        sys.exit(1)


if __name__ == "__main__":
    main()
