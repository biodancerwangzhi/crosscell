#!/usr/bin/env python3
"""
easySCF H5AD→RDS 测试 — 步骤1: H5AD → .h5 (自定义格式)

使用 easySCFpy.saveH5() 将标准 H5AD 文件转换为 easySCF 的自定义 .h5 格式。
"""

import time
import traceback
import json
from pathlib import Path

DATA_DIR = Path('/benchmark/data/generated')
TMP_DIR = Path('/benchmark/results/easyscf_h2r_tmp')
TMP_DIR.mkdir(parents=True, exist_ok=True)

# 13 个 H5AD 测试数据集
H5AD_FILES = [
    'scanpy_pbmc3k.h5ad',
    'scanpy_pbmc3k_processed.h5ad',
    'scvelo_dentategyrus.h5ad',
    'scvelo_pancreas.h5ad',
    'cellxgene_pbmc_15k.h5ad',
    'cellxgene_heart_23k.h5ad',
    'cellxgene_brain_40k.h5ad',
    'squidpy_visium_hne.h5ad',
    'squidpy_mibitof.h5ad',
    'squidpy_slideseqv2.h5ad',
    'squidpy_seqfish.h5ad',
    'squidpy_merfish.h5ad',
    'squidpy_imc.h5ad',
]

results = []

for fname in H5AD_FILES:
    h5ad_path = DATA_DIR / fname
    h5_path = TMP_DIR / fname.replace('.h5ad', '_easyscf.h5')
    test_id = fname.replace('.h5ad', '')

    print(f'\n{"="*60}')
    print(f'Testing: {test_id}')
    print(f'  Input:  {h5ad_path}')
    print(f'  Output: {h5_path}')

    if not h5ad_path.exists():
        print(f'  ❌ File not found')
        results.append({'test_id': test_id, 'step1_ok': False, 'step1_error': 'file not found'})
        continue

    try:
        import scanpy as sc
        from easySCFpy import saveH5

        t0 = time.time()
        adata = sc.read_h5ad(str(h5ad_path))
        print(f'  read_h5ad: {adata.n_obs} cells × {adata.n_vars} genes ({time.time()-t0:.1f}s)')

        t1 = time.time()
        saveH5(adata, str(h5_path))
        elapsed = time.time() - t1
        
        ok = h5_path.exists() and h5_path.stat().st_size > 0
        print(f'  saveH5: {"✅" if ok else "❌"} ({elapsed:.1f}s, {h5_path.stat().st_size/1024/1024:.1f} MB)')
        results.append({
            'test_id': test_id,
            'step1_ok': ok,
            'step1_time': round(elapsed, 2),
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'h5_size_mb': round(h5_path.stat().st_size / 1024 / 1024, 1),
        })
    except Exception as e:
        print(f'  ❌ Error: {e}')
        traceback.print_exc()
        results.append({
            'test_id': test_id,
            'step1_ok': False,
            'step1_error': str(e),
        })

# 汇总
print(f'\n{"="*60}')
print(f'Step 1 Summary (H5AD → .h5):')
ok_count = sum(1 for r in results if r.get('step1_ok'))
print(f'  {ok_count}/{len(results)} succeeded')
for r in results:
    status = '✅' if r.get('step1_ok') else '❌'
    err = f" — {r.get('step1_error', '')}" if not r.get('step1_ok') else ''
    print(f'  {status} {r["test_id"]}{err}')

# 保存中间结果
out_file = TMP_DIR / 'step1_results.json'
with open(out_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f'\nResults saved to {out_file}')
