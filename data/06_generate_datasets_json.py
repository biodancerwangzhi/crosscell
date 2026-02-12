#!/usr/bin/env python3
"""Generate benchmark/data/datasets.json with metadata for all test datasets.

使用方法 (在 benchmark 容器内):
  python3 /benchmark/data/06_generate_datasets_json.py
"""

import json
import os
import csv
from datetime import datetime

def read_rds_summaries():
    """Read RDS dataset info from CSV summaries."""
    datasets = []
    
    for version in ["v4", "v5"]:
        csv_path = f"/benchmark/data/generated/seurat_{version}_summary.csv"
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found, skipping")
            continue
        
        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                datasets.append({
                    "name": row["dataset"],
                    "file": row["file"],
                    "path": f"data/generated/{row['file']}",
                    "format": "rds",
                    "seurat_version": row["version"],
                    "n_cells": int(row["n_cells"]),
                    "n_genes": int(row["n_genes"]),
                    "assays": row["assays"],
                    "default_assay": row["default_assay"],
                    "assay_class": row["assay_class"],
                    "has_spatial": row["has_spatial"] == "TRUE",
                    "reductions": row["reductions"] if row["reductions"] else None,
                    "file_size_mb": float(row["file_size_mb"]),
                    "source": "SeuratData"
                })
    
    # Also check for processed files and raw files not in CSV
    generated_dir = "/benchmark/data/generated"
    if os.path.exists(generated_dir):
        for f in sorted(os.listdir(generated_dir)):
            if not f.endswith(".rds"):
                continue
            # Check if already in list
            if any(d["file"] == f for d in datasets):
                continue
            # Get size
            fpath = os.path.join(generated_dir, f)
            size_mb = os.path.getsize(fpath) / (1024 * 1024)
            version = "V4" if "_v4_" in f else "V5"
            is_processed = "_processed.rds" in f
            dataset_name = f.replace(f"seurat_{version.lower()}_", "").replace("_processed.rds", "").replace("_raw.rds", "")
            entry = {
                "name": dataset_name,
                "file": f,
                "path": f"data/generated/{f}",
                "format": "rds",
                "seurat_version": version,
                "type": "processed" if is_processed else "raw",
                "n_cells": None,
                "n_genes": None,
                "file_size_mb": round(size_mb, 2),
                "source": "SeuratData"
            }
            datasets.append(entry)
    
    return datasets


def read_h5ad_metadata():
    """Read H5AD file metadata using anndata."""
    import anndata
    
    datasets = []
    h5ad_dir = "/benchmark/data/generated"
    
    if not os.path.exists(h5ad_dir):
        print(f"Warning: {h5ad_dir} not found")
        return datasets
    
    # Source mapping based on filename prefix
    source_map = {
        "scanpy_": "scanpy",
        "scvelo_": "scvelo",
        "squidpy_": "squidpy",
        "cellxgene_": "CELLxGENE"
    }
    
    for fname in sorted(os.listdir(h5ad_dir)):
        if not fname.endswith(".h5ad"):
            continue
        
        fpath = os.path.join(h5ad_dir, fname)
        size_mb = os.path.getsize(fpath) / (1024 * 1024)
        
        # Determine source
        source = "unknown"
        for prefix, src in source_map.items():
            if fname.startswith(prefix):
                source = src
                break
        
        # Read with anndata
        try:
            adata = anndata.read_h5ad(fpath, backed="r")
            n_cells, n_genes = adata.shape
            
            # Get available components
            components = []
            if adata.obs is not None and len(adata.obs.columns) > 0:
                components.append("obs")
            if adata.var is not None and len(adata.var.columns) > 0:
                components.append("var")
            if adata.obsm is not None and len(adata.obsm) > 0:
                components.append("obsm")
            if adata.obsp is not None and len(adata.obsp) > 0:
                components.append("obsp")
            if adata.uns is not None and len(adata.uns) > 0:
                components.append("uns")
            if adata.layers is not None and len(adata.layers) > 0:
                components.append("layers")
            
            obsm_keys = list(adata.obsm.keys()) if adata.obsm is not None else []
            layer_keys = list(adata.layers.keys()) if adata.layers is not None else []
            
            adata.file.close()
            
            dataset_name = fname.replace(".h5ad", "")
            
            datasets.append({
                "name": dataset_name,
                "file": fname,
                "path": f"data/generated/{fname}",
                "format": "h5ad",
                "n_cells": n_cells,
                "n_genes": n_genes,
                "file_size_mb": round(size_mb, 2),
                "source": source,
                "components": components,
                "obsm_keys": obsm_keys,
                "layer_keys": layer_keys
            })
            print(f"  ✓ {fname}: {n_cells} cells × {n_genes} genes ({size_mb:.1f} MB)")
        except Exception as e:
            print(f"  ✗ {fname}: {e}")
            datasets.append({
                "name": fname.replace(".h5ad", ""),
                "file": fname,
                "path": f"data/generated/{fname}",
                "format": "h5ad",
                "n_cells": None,
                "n_genes": None,
                "file_size_mb": round(size_mb, 2),
                "source": source,
                "error": str(e)
            })
    
    return datasets


def classify_scale(n_cells):
    """Classify dataset scale based on cell count."""
    if n_cells is None:
        return "unknown"
    if n_cells < 5000:
        return "small"
    elif n_cells < 20000:
        return "medium"
    else:
        return "large"


def main():
    print("=" * 60)
    print("Generating datasets.json")
    print("=" * 60)
    
    all_datasets = []
    
    # 1. RDS datasets
    print("\n📦 Reading RDS datasets from CSV summaries...")
    rds_datasets = read_rds_summaries()
    print(f"  Found {len(rds_datasets)} RDS files")
    for d in rds_datasets:
        d["scale"] = classify_scale(d.get("n_cells"))
    all_datasets.extend(rds_datasets)
    
    # 2. H5AD datasets
    print("\n📦 Reading H5AD datasets...")
    h5ad_datasets = read_h5ad_metadata()
    print(f"  Found {len(h5ad_datasets)} H5AD files")
    for d in h5ad_datasets:
        d["scale"] = classify_scale(d.get("n_cells"))
    all_datasets.extend(h5ad_datasets)
    
    # Build output
    output = {
        "version": "1.0",
        "generated_at": datetime.now().isoformat(),
        "summary": {
            "total_datasets": len(all_datasets),
            "rds_count": len(rds_datasets),
            "h5ad_count": len(h5ad_datasets),
            "by_source": {},
            "by_scale": {"small": 0, "medium": 0, "large": 0, "unknown": 0}
        },
        "datasets": all_datasets
    }
    
    # Count by source
    for d in all_datasets:
        src = d.get("source", "unknown")
        output["summary"]["by_source"][src] = output["summary"]["by_source"].get(src, 0) + 1
    
    # Count by scale
    for d in all_datasets:
        scale = d.get("scale", "unknown")
        output["summary"]["by_scale"][scale] = output["summary"]["by_scale"].get(scale, 0) + 1
    
    # Write output
    out_path = "/benchmark/results/datasets.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ Written to {out_path}")
    print(f"   Total: {len(all_datasets)} datasets ({len(rds_datasets)} RDS + {len(h5ad_datasets)} H5AD)")


if __name__ == "__main__":
    main()
