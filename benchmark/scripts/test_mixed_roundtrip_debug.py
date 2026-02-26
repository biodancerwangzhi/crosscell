#!/usr/bin/env python3
"""
Debug mixed roundtrip failures for anndataR and easySCF.

Tests on a single small dataset (pbmc3k V4 raw):
1. Each R tool converts RDS → H5AD
2. Inspect the generated H5AD structure (h5py)
3. Try reading with anndata
4. Try CrossCell H5AD → RDS conversion
5. Report detailed errors at each step

Run in Docker:
  docker-compose run --rm dev python3 benchmark/scripts/test_mixed_roundtrip_debug.py
"""

import subprocess, json, sys, os, traceback
from pathlib import Path

DATA_DIR = Path('/benchmark/data/generated')
TMP = Path('/tmp/mixed_rt_debug')
TMP.mkdir(parents=True, exist_ok=True)

RDS_FILE = DATA_DIR / 'seurat_v4_pbmc3k_raw.rds'

def run_cmd(cmd, timeout=300):
    """Run command, return (success, stdout, stderr)."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode == 0, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return False, '', 'TIMEOUT'
    except Exception as e:
        return False, '', str(e)

def inspect_h5ad(path, label):
    """Inspect H5AD/H5 file structure with h5py."""
    import h5py
    print(f'\n  --- h5py inspection of {label} ---')
    if not Path(path).exists():
        print(f'  FILE NOT FOUND: {path}')
        return False
    
    fsize = Path(path).stat().st_size
    print(f'  File size: {fsize:,} bytes')
    
    try:
        with h5py.File(path, 'r') as f:
            print(f'  Top-level keys: {list(f.keys())}')
            for key in f.keys():
                obj = f[key]
                if hasattr(obj, 'shape'):
                    print(f'    {key}: dataset, shape={obj.shape}, dtype={obj.dtype}')
                elif hasattr(obj, 'keys'):
                    subkeys = list(obj.keys())[:10]
                    print(f'    {key}: group, subkeys={subkeys}{"..." if len(obj.keys()) > 10 else ""}')
                    # Check X structure
                    if key == 'X':
                        for sk in obj.keys():
                            sub = obj[sk]
                            if hasattr(sub, 'shape'):
                                print(f'      {sk}: shape={sub.shape}, dtype={sub.dtype}')
                else:
                    print(f'    {key}: {type(obj)}')
            
            # Check encoding
            if 'encoding-type' in f.attrs:
                print(f'  encoding-type: {f.attrs["encoding-type"]}')
            if 'encoding-version' in f.attrs:
                print(f'  encoding-version: {f.attrs["encoding-version"]}')
            
            # Check X encoding
            if 'X' in f:
                x = f['X']
                if 'encoding-type' in x.attrs:
                    print(f'  X encoding-type: {x.attrs["encoding-type"]}')
        return True
    except Exception as e:
        print(f'  h5py ERROR: {e}')
        return False

def try_anndata_read(path, label):
    """Try reading with anndata."""
    import anndata as ad
    print(f'\n  --- anndata read of {label} ---')
    try:
        adata = ad.read_h5ad(path)
        print(f'  Success: {adata.shape[0]} cells × {adata.shape[1]} genes')
        print(f'  X type: {type(adata.X).__name__}', end='')
        if hasattr(adata.X, 'format'):
            print(f' ({adata.X.format})', end='')
        if hasattr(adata.X, 'dtype'):
            print(f' dtype={adata.X.dtype}', end='')
        print()
        print(f'  obs columns: {list(adata.obs.columns)[:5]}...')
        print(f'  var columns: {list(adata.var.columns)[:5]}...')
        print(f'  layers: {list(adata.layers.keys())}')
        return True
    except Exception as e:
        print(f'  anndata read FAILED: {e}')
        return False

def try_crosscell_h2r(h5ad_path, rds_out, label):
    """Try CrossCell H5AD → RDS."""
    print(f'\n  --- CrossCell H5AD→RDS of {label} ---')
    ok, stdout, stderr = run_cmd(
        ['crosscell', 'convert', '-i', str(h5ad_path), '-o', str(rds_out), '-f', 'seurat']
    )
    if ok:
        print(f'  Success: output {Path(rds_out).stat().st_size:,} bytes')
    else:
        print(f'  FAILED')
        if stdout.strip():
            print(f'  stdout: {stdout[:500]}')
        if stderr.strip():
            print(f'  stderr: {stderr[:500]}')
    return ok, stderr


# ══════════════════════════════════════════════
# Test 1: anndataR
# ══════════════════════════════════════════════
def test_anndataR():
    print('=' * 60)
    print('TEST: anndataR RDS → H5AD → (CrossCell) RDS')
    print('=' * 60)
    
    h5ad_out = str(TMP / 'anndataR_pbmc3k.h5ad')
    rds_out = str(TMP / 'anndataR_pbmc3k_rt.rds')
    
    # Step 1: anndataR converts RDS → H5AD
    print('\n[Step 1] anndataR: RDS → H5AD')
    r_script = f'''
suppressPackageStartupMessages({{
    library(anndataR)
    library(Seurat)
}})
cat("anndataR version:", as.character(packageVersion("anndataR")), "\\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\\n")

obj <- readRDS("{RDS_FILE}")
tryCatch({{obj <- UpdateSeuratObject(obj)}}, error=function(e){{}})
cat("Seurat object:", ncol(obj), "cells x", nrow(obj), "genes\\n")

if (file.exists("{h5ad_out}")) file.remove("{h5ad_out}")
ad <- as_AnnData(obj)
cat("AnnData object created\\n")
cat("X shape:", dim(ad)[1], "x", dim(ad)[2], "\\n")
write_h5ad(ad, "{h5ad_out}")
cat("Written to:", "{h5ad_out}", "\\n")
cat("File size:", file.size("{h5ad_out}"), "bytes\\n")
'''
    ok, stdout, stderr = run_cmd(['Rscript', '-e', r_script])
    print(f'  R exit: {"OK" if ok else "FAILED"}')
    if stdout.strip():
        for line in stdout.strip().split('\n'):
            print(f'  R stdout: {line}')
    if stderr.strip():
        # Filter out R startup messages
        for line in stderr.strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('Loading required') and not line.startswith('Attaching'):
                print(f'  R stderr: {line}')
    
    if not ok or not Path(h5ad_out).exists():
        print('  ❌ anndataR conversion failed, cannot continue')
        return
    
    # Step 2: Inspect H5AD
    inspect_h5ad(h5ad_out, 'anndataR output')
    
    # Step 3: Try anndata read
    try_anndata_read(h5ad_out, 'anndataR output')
    
    # Step 4: Try CrossCell H5AD → RDS
    ok_cc, stderr_cc = try_crosscell_h2r(h5ad_out, rds_out, 'anndataR output')
    
    print(f'\n  CONCLUSION: anndataR → CrossCell H2R: {"✅ OK" if ok_cc else "❌ FAILED"}')
    if not ok_cc:
        print(f'  Error summary: {stderr_cc[:300]}')


# ══════════════════════════════════════════════
# Test 2: easySCF
# ══════════════════════════════════════════════
def test_easySCF():
    print('\n' + '=' * 60)
    print('TEST: easySCF RDS → H5AD → (CrossCell) RDS')
    print('=' * 60)
    
    h5ad_out = str(TMP / 'easySCF_pbmc3k.h5')
    rds_out = str(TMP / 'easySCF_pbmc3k_rt.rds')
    
    # Step 1: easySCF converts RDS → H5
    print('\n[Step 1] easySCF: RDS → H5')
    r_script = f'''
suppressPackageStartupMessages({{
    library(easySCFr)
    library(Seurat)
}})
cat("easySCFr loaded\\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\\n")

obj <- readRDS("{RDS_FILE}")
tryCatch({{obj <- UpdateSeuratObject(obj)}}, error=function(e){{}})
cat("Seurat object:", ncol(obj), "cells x", nrow(obj), "genes\\n")

if (file.exists("{h5ad_out}")) file.remove("{h5ad_out}")
saveH5(obj, "{h5ad_out}")
cat("Written to:", "{h5ad_out}", "\\n")
if (file.exists("{h5ad_out}")) {{
    cat("File size:", file.size("{h5ad_out}"), "bytes\\n")
}} else {{
    cat("ERROR: output file not created\\n")
}}
'''
    ok, stdout, stderr = run_cmd(['Rscript', '-e', r_script])
    print(f'  R exit: {"OK" if ok else "FAILED"}')
    if stdout.strip():
        for line in stdout.strip().split('\n'):
            print(f'  R stdout: {line}')
    if stderr.strip():
        for line in stderr.strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('Loading required') and not line.startswith('Attaching'):
                print(f'  R stderr: {line}')
    
    if not ok or not Path(h5ad_out).exists():
        print('  ❌ easySCF conversion failed, cannot continue')
        return
    
    # Step 2: Inspect H5
    inspect_h5ad(h5ad_out, 'easySCF output')
    
    # Step 3: Try anndata read
    try_anndata_read(h5ad_out, 'easySCF output')
    
    # Step 4: Try CrossCell H5AD → RDS
    ok_cc, stderr_cc = try_crosscell_h2r(h5ad_out, rds_out, 'easySCF output')
    
    print(f'\n  CONCLUSION: easySCF → CrossCell H2R: {"✅ OK" if ok_cc else "❌ FAILED"}')
    if not ok_cc:
        print(f'  Error summary: {stderr_cc[:300]}')


# ══════════════════════════════════════════════
# Test 3: Also test Zellkonverter as control
# ══════════════════════════════════════════════
def test_zellkonverter_control():
    print('\n' + '=' * 60)
    print('CONTROL: Zellkonverter RDS → H5AD → (CrossCell) RDS')
    print('=' * 60)
    
    h5ad_out = str(TMP / 'zellkonverter_pbmc3k.h5ad')
    rds_out = str(TMP / 'zellkonverter_pbmc3k_rt.rds')
    
    print('\n[Step 1] Zellkonverter: RDS → H5AD')
    r_script = f'''
suppressPackageStartupMessages({{
    library(zellkonverter)
    library(Seurat)
    library(SingleCellExperiment)
}})
cat("zellkonverter version:", as.character(packageVersion("zellkonverter")), "\\n")

obj <- readRDS("{RDS_FILE}")
tryCatch({{obj <- UpdateSeuratObject(obj)}}, error=function(e){{}})
cat("Seurat object:", ncol(obj), "cells x", nrow(obj), "genes\\n")

sce <- as.SingleCellExperiment(obj)
cat("SCE:", ncol(sce), "cells x", nrow(sce), "genes\\n")
writeH5AD(sce, "{h5ad_out}")
cat("Written to:", "{h5ad_out}", "\\n")
cat("File size:", file.size("{h5ad_out}"), "bytes\\n")
'''
    ok, stdout, stderr = run_cmd(['Rscript', '-e', r_script])
    print(f'  R exit: {"OK" if ok else "FAILED"}')
    if stdout.strip():
        for line in stdout.strip().split('\n'):
            print(f'  R stdout: {line}')
    
    if not ok or not Path(h5ad_out).exists():
        print('  ❌ Zellkonverter conversion failed')
        return
    
    inspect_h5ad(h5ad_out, 'Zellkonverter output')
    try_anndata_read(h5ad_out, 'Zellkonverter output')
    ok_cc, stderr_cc = try_crosscell_h2r(h5ad_out, rds_out, 'Zellkonverter output')
    print(f'\n  CONCLUSION: Zellkonverter → CrossCell H2R: {"✅ OK" if ok_cc else "❌ FAILED"}')


if __name__ == '__main__':
    if not RDS_FILE.exists():
        print(f'ERROR: {RDS_FILE} not found. Run inside Docker container.')
        sys.exit(1)
    
    test_anndataR()
    test_easySCF()
    test_zellkonverter_control()
    
    print('\n' + '=' * 60)
    print('ALL TESTS COMPLETE')
    print('=' * 60)
