#!/usr/bin/env python3
"""Inspect H5AD files to diagnose format issues."""
import h5py
import sys

files = [
    '/benchmark/data/generated/scanpy_pbmc3k.h5ad',
    '/benchmark/data/generated/scanpy_pbmc3k_processed.h5ad',
    '/benchmark/data/generated/scvelo_pancreas.h5ad',
    '/benchmark/data/generated/scvelo_dentategyrus.h5ad',
    '/benchmark/data/generated/squidpy_imc.h5ad',
]

for fpath in files:
    name = fpath.split("/")[-1]
    print("=" * 60)
    print(f"FILE: {name}")
    print("=" * 60)
    try:
        h = h5py.File(fpath, 'r')
        print(f"  Top-level keys: {list(h.keys())}")
        print(f"  Top-level attrs: {dict(h.attrs)}")
        
        # Check X
        if 'X' in h:
            x = h['X']
            if isinstance(x, h5py.Dataset):
                print(f"  X: Dense Dataset, shape={x.shape}, dtype={x.dtype}")
            elif isinstance(x, h5py.Group):
                print(f"  X: Group (sparse), keys={list(x.keys())}")
                for k, v in x.attrs.items():
                    print(f"    X.attrs[{k}] = {v}")
        
        # Check obs
        if 'obs' in h:
            obs = h['obs']
            if isinstance(obs, h5py.Group):
                print(f"  obs: Group, keys={list(obs.keys())[:10]}")
                for k, v in obs.attrs.items():
                    print(f"    obs.attrs[{k}] = {v}")
            elif isinstance(obs, h5py.Dataset):
                print(f"  obs: Dataset, shape={obs.shape}, dtype={obs.dtype}")
        
        # Check var
        if 'var' in h:
            var = h['var']
            if isinstance(var, h5py.Group):
                print(f"  var: Group, keys={list(var.keys())[:10]}")
            elif isinstance(var, h5py.Dataset):
                print(f"  var: Dataset, shape={var.shape}, dtype={var.dtype}")
        
        # Check layers
        if 'layers' in h:
            layers = h['layers']
            print(f"  layers: keys={list(layers.keys())}")
        
        # Check obsm
        if 'obsm' in h:
            obsm = h['obsm']
            print(f"  obsm: keys={list(obsm.keys())}")
        
        h.close()
    except Exception as e:
        print(f"  ERROR: {e}")
    print()
