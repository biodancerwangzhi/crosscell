#!/usr/bin/env python3
"""检查下载的 H5AD 文件是否完整"""

import os
import sys

def check_file(filepath, expected_name):
    """检查单个 H5AD 文件"""
    if not os.path.exists(filepath):
        print(f"❌ 文件不存在: {expected_name}")
        return False
    
    size_mb = os.path.getsize(filepath) / 1024 / 1024
    
    try:
        import h5py
        with h5py.File(filepath, "r") as f:
            keys = list(f.keys())
            # 检查基本结构
            has_x = "X" in keys
            has_obs = "obs" in keys
            has_var = "var" in keys
            
            # 尝试获取细胞数和基因数
            n_obs = f["obs"].attrs.get("_index", "unknown")
            
        if has_x and has_obs and has_var:
            print(f"✅ {expected_name}: {size_mb:.1f} MB, keys={keys}")
            return True
        else:
            print(f"⚠️  {expected_name}: {size_mb:.1f} MB, 结构不完整 keys={keys}")
            return False
            
    except Exception as e:
        print(f"❌ {expected_name}: {size_mb:.1f} MB, 读取失败: {e}")
        return False

def main():
    data_dir = "/benchmark/data/generated"
    
    # 10 个 benchmark 数据集（排除 retina_244k：异常 dense 结构）
    files = [
        ("687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad", "cellxgene_combat_pbmc_836k.h5ad"),
        ("4cb45d80-499a-48ae-a056-c71ac3552c94.h5ad", "cellxgene_hlca_core_585k.h5ad"),
        ("a7f6822d-0e0e-451e-9858-af81392fcb9b.h5ad", "cellxgene_heart_486k.h5ad"),
        ("fa0f9f5e-e7fb-42a4-a801-3410b0923ccc.h5ad", "cellxgene_kidney_atacseq_37k.h5ad"),
        ("f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad", "cellxgene_gut_428k.h5ad"),
        ("62e4e805-9a9a-42ca-a6e7-85df258b382c.h5ad", "cellxgene_tabula_liver_22k.h5ad"),
        ("412352dd-a919-4d8e-9f74-e210627328b5.h5ad", "cellxgene_brain_multiome_102k.h5ad"),
        ("9198437d-f781-4215-bf4f-32c3177e57df.h5ad", "cellxgene_brain_dlpfc_172k.h5ad"),
        ("72e6a74d-8804-4305-81d7-a12b93f2ac0d.h5ad", "cellxgene_pancreas_122k.h5ad"),
        ("a6670528-1e94-4920-a78f-49a76b96ca9e.h5ad", "cellxgene_skin_bcc_10k.h5ad"),
    ]
    
    print("=" * 60)
    print("检查 CELLxGENE H5AD 文件完整性 (10 datasets)")
    print("=" * 60)
    print("注：retina_244k 已排除（异常 dense 结构，18 GB 文件需 62 GB 内存）")
    print()
    
    success = 0
    failed = 0
    
    for original_name, target_name in files:
        # 先检查目标文件名，再检查原始文件名
        target_path = os.path.join(data_dir, target_name)
        original_path = os.path.join(data_dir, original_name)
        
        if os.path.exists(target_path):
            if check_file(target_path, target_name):
                success += 1
            else:
                failed += 1
        elif os.path.exists(original_path):
            if check_file(original_path, f"{original_name} (需重命名为 {target_name})"):
                success += 1
            else:
                failed += 1
        else:
            print(f"❌ 文件不存在: {target_name} 或 {original_name}")
            failed += 1
    
    print("=" * 60)
    print(f"总计: {success + failed}, 成功: {success}, 失败: {failed}")
    
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
