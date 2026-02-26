#!/usr/bin/env python3
"""
Package benchmark data and results for Zenodo upload.

按类别分包，方便用户按需下载：
  1. seurat_v4_datasets.tar.gz      — Seurat V4 RDS 文件 (~1.3 GB)
  2. seurat_v5_datasets.tar.gz      — Seurat V5 RDS 文件 (~0.7 GB)
  3. h5ad_small_datasets.tar.gz     — scanpy + scvelo H5AD (~0.1 GB)
  4. h5ad_spatial_datasets.tar.gz   — squidpy 空间数据 H5AD (~0.7 GB)
  5. h5ad_cellxgene_datasets.tar.gz — CELLxGENE 小规模 H5AD (≤40k cells, ~0.7 GB)
  6. benchmark_results.tar.gz       — 预计算结果 + 脚本 (~1 MB)

注意：10 个 CELLxGENE benchmark 大规模数据集（10k-836k cells）不打包到 Zenodo，
因为它们可以从 CELLxGENE Discover 免费下载。脚本会生成 cellxgene_sources.json
记录下载 URL 和 dataset ID，供用户自行下载。

Usage (from project root on host):
  python3 data/07_package_for_zenodo.py

Output:
  data/zenodo/
  ├── seurat_v4_datasets.tar.gz
  ├── seurat_v5_datasets.tar.gz
  ├── h5ad_small_datasets.tar.gz
  ├── h5ad_spatial_datasets.tar.gz
  ├── h5ad_cellxgene_datasets.tar.gz
  ├── benchmark_results.tar.gz
  ├── cellxgene_sources.json
  ├── README.md
  └── manifest.json
"""
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent  # project root
ZENODO_DIR = ROOT / "data" / "zenodo"
DATA_DIR = ROOT / "data" / "generated"
RESULTS_DIR = ROOT / "benchmark" / "results"
SCRIPTS_DIR = ROOT / "benchmark" / "scripts"
BENCHMARK_DIR = ROOT / "benchmark"


def sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def file_size_mb(filepath):
    return os.path.getsize(filepath) / (1024 * 1024)


def classify_data_files():
    """Classify data files into groups for separate archives.
    
    Large CELLxGENE benchmark datasets (>40k cells) are excluded from Zenodo
    packaging — they can be downloaded directly from CELLxGENE Discover.
    Only the 3 small CELLxGENE files (pbmc_15k, heart_23k, brain_40k) are included.
    """
    # 10 large CELLxGENE benchmark files to EXCLUDE from Zenodo tar.gz
    # (users download these directly from CELLxGENE Discover)
    CELLXGENE_BENCHMARK_EXCLUDE = {
        "cellxgene_combat_pbmc_836k.h5ad",
        "cellxgene_hlca_core_585k.h5ad",
        "cellxgene_heart_486k.h5ad",
        "cellxgene_gut_428k.h5ad",
        "cellxgene_retina_244k.h5ad",
        "cellxgene_brain_dlpfc_172k.h5ad",
        "cellxgene_pancreas_122k.h5ad",
        "cellxgene_brain_multiome_102k.h5ad",
        "cellxgene_kidney_atacseq_37k.h5ad",
        "cellxgene_tabula_liver_22k.h5ad",
        "cellxgene_skin_bcc_10k.h5ad",
        # Also exclude the 1.8M scalability dataset
        "cellxgene_immune_1.8M.h5ad",
    }
    
    groups = {
        "seurat_v4_datasets": [],
        "seurat_v5_datasets": [],
        "h5ad_small_datasets": [],      # scanpy + scvelo
        "h5ad_spatial_datasets": [],     # squidpy
        "h5ad_cellxgene_datasets": [],   # cellxgene (small only: pbmc_15k, heart_23k, brain_40k)
    }

    if not DATA_DIR.exists():
        print(f"WARNING: {DATA_DIR} does not exist")
        return groups

    excluded_count = 0
    for f in sorted(DATA_DIR.iterdir()):
        name = f.name
        if name.startswith("seurat_v4_") and name.endswith(".rds"):
            groups["seurat_v4_datasets"].append(f)
        elif name.startswith("seurat_v5_") and name.endswith(".rds"):
            groups["seurat_v5_datasets"].append(f)
        elif name.endswith(".h5ad"):
            if name in CELLXGENE_BENCHMARK_EXCLUDE:
                excluded_count += 1
                continue
            if name.startswith(("scanpy_", "scvelo_")):
                groups["h5ad_small_datasets"].append(f)
            elif name.startswith("squidpy_"):
                groups["h5ad_spatial_datasets"].append(f)
            elif name.startswith("cellxgene_"):
                groups["h5ad_cellxgene_datasets"].append(f)
            else:
                # 未知来源的 h5ad 放到 small
                groups["h5ad_small_datasets"].append(f)

    if excluded_count > 0:
        print(f"  ℹ️  Excluded {excluded_count} large CELLxGENE files (download from CELLxGENE Discover)")

    return groups


def create_archive(name, files):
    """Create a tar.gz archive from a list of files."""
    archive = ZENODO_DIR / f"{name}.tar.gz"
    if not files:
        print(f"  {name}: no files, skipping")
        return None

    rel_paths = [str(f.relative_to(ROOT)) for f in files]
    total_mb = sum(file_size_mb(f) for f in files)

    print(f"\n  {name}.tar.gz")
    print(f"    Files: {len(files)}, uncompressed: {total_mb:.1f} MB")

    cmd = ["tar", "czf", str(archive), "-C", str(ROOT)] + rel_paths
    subprocess.run(cmd, check=True)

    archive_mb = file_size_mb(archive)
    print(f"    Archive: {archive_mb:.1f} MB")
    return archive


def create_results_archive():
    """Create archive of results, scripts, and configs."""
    all_files = []

    # Result JSONs
    if RESULTS_DIR.exists():
        for f in sorted(RESULTS_DIR.iterdir()):
            if f.suffix == ".json":
                all_files.append(f)
        # Include cellxgene_11/ subdirectory results
        cellxgene_dir = RESULTS_DIR / "cellxgene_11"
        if cellxgene_dir.exists():
            for f in sorted(cellxgene_dir.iterdir()):
                if f.suffix == ".json":
                    all_files.append(f)

    # Benchmark scripts
    if SCRIPTS_DIR.exists():
        for f in sorted(SCRIPTS_DIR.iterdir()):
            if f.suffix in (".py", ".sh", ".R") and not f.name.startswith("__"):
                all_files.append(f)

    # Config files
    config_dir = BENCHMARK_DIR / "data"
    if config_dir.exists():
        for f in sorted(config_dir.iterdir()):
            if f.suffix == ".json":
                all_files.append(f)

    # Key benchmark files
    for name in ("README.md", "Dockerfile.benchmark", "docker-compose.benchmark.yml"):
        f = BENCHMARK_DIR / name
        if f.exists():
            all_files.append(f)

    return create_archive("benchmark_results", all_files)


def create_cellxgene_sources_json():
    """Generate cellxgene_sources.json with download URLs for large CELLxGENE datasets.
    
    These datasets are NOT included in Zenodo archives (too large).
    Users download them directly from CELLxGENE Discover using these URLs.
    """
    sources = {
        "description": "CELLxGENE benchmark datasets — download directly from CELLxGENE Discover",
        "note": "These files are NOT included in Zenodo archives. Use the URLs below to download.",
        "excluded_from_benchmark": {
            "cellxgene_retina_244k": {
                "reason": "Abnormal dense matrix structure (18 GB file, inspect requires 62 GB RAM)",
                "url": "https://datasets.cellxgene.cziscience.com/d4f1a799-9616-4089-90e9-f3fe843f01e5.h5ad"
            }
        },
        "datasets": [
            {
                "name": "cellxgene_skin_bcc_10k",
                "filename": "cellxgene_skin_bcc_10k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/a6670528-1e94-4920-a78f-49a76b96ca9e.h5ad",
                "cells": 9841, "genes": 26886, "tissue": "skin", "assay": "multi-assay",
                "disease": "basal cell carcinoma", "size_mb": 89
            },
            {
                "name": "cellxgene_tabula_liver_22k",
                "filename": "cellxgene_tabula_liver_22k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/62e4e805-9a9a-42ca-a6e7-85df258b382c.h5ad",
                "cells": 22214, "genes": 60606, "tissue": "liver", "assay": "Smart-seq2 + 10x",
                "disease": "normal", "size_mb": 796
            },
            {
                "name": "cellxgene_kidney_atacseq_37k",
                "filename": "cellxgene_kidney_atacseq_37k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/fa0f9f5e-e7fb-42a4-a801-3410b0923ccc.h5ad",
                "cells": 37747, "genes": 19276, "tissue": "kidney", "assay": "10x scATAC-seq",
                "disease": "normal", "size_mb": 898
            },
            {
                "name": "cellxgene_brain_multiome_102k",
                "filename": "cellxgene_brain_multiome_102k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/412352dd-a919-4d8e-9f74-e210627328b5.h5ad",
                "cells": 101924, "genes": 35451, "tissue": "brain", "assay": "10x multiome",
                "disease": "normal", "size_mb": 1142
            },
            {
                "name": "cellxgene_pancreas_122k",
                "filename": "cellxgene_pancreas_122k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/72e6a74d-8804-4305-81d7-a12b93f2ac0d.h5ad",
                "cells": 121916, "genes": 32356, "tissue": "pancreas", "assay": "multi-assay",
                "disease": "normal", "size_mb": 915
            },
            {
                "name": "cellxgene_brain_dlpfc_172k",
                "filename": "cellxgene_brain_dlpfc_172k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/9198437d-f781-4215-bf4f-32c3177e57df.h5ad",
                "cells": 172120, "genes": 37490, "tissue": "brain", "assay": "10x 3' v3",
                "disease": "normal", "size_mb": 2234
            },
            {
                "name": "cellxgene_gut_428k",
                "filename": "cellxgene_gut_428k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad",
                "cells": 428469, "genes": 32383, "tissue": "gut", "assay": "10x 3'/5' v2",
                "disease": "normal", "size_mb": 5462
            },
            {
                "name": "cellxgene_heart_486k",
                "filename": "cellxgene_heart_486k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/a7f6822d-0e0e-451e-9858-af81392fcb9b.h5ad",
                "cells": 486134, "genes": 32383, "tissue": "heart", "assay": "10x 3' v2/v3",
                "disease": "normal", "size_mb": 2816
            },
            {
                "name": "cellxgene_hlca_core_585k",
                "filename": "cellxgene_hlca_core_585k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/4cb45d80-499a-48ae-a056-c71ac3552c94.h5ad",
                "cells": 584944, "genes": 27402, "tissue": "lung", "assay": "multi-assay",
                "disease": "normal", "size_mb": 5602
            },
            {
                "name": "cellxgene_combat_pbmc_836k",
                "filename": "cellxgene_combat_pbmc_836k.h5ad",
                "url": "https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad",
                "cells": 836148, "genes": 36306, "tissue": "blood", "assay": "10x 5' v1",
                "disease": "COVID-19/sepsis/flu", "size_mb": 5337
            }
        ],
        "download_script": "python3 data/05_download_h5ad_datasets.py --source cellxgene-benchmark"
    }
    
    out_path = ZENODO_DIR / "cellxgene_sources.json"
    with open(out_path, "w") as f:
        json.dump(sources, f, indent=2, ensure_ascii=False)
    print(f"\nCELLxGENE sources: {out_path}")
    return out_path


def create_manifest(groups, archives):
    """Create manifest.json with file listing and checksums."""
    manifest = {
        "created": datetime.now().isoformat(),
        "project": "CrossCell Benchmark",
        "description": "Benchmark data and results for CrossCell paper",
        "archives": {},
        "datasets": {},
    }

    for name, path in archives.items():
        if path and path.exists():
            manifest["archives"][name] = {
                "sha256": sha256(path),
                "size_mb": round(file_size_mb(path), 2),
            }

    for group_name, files in groups.items():
        manifest["datasets"][group_name] = []
        for f in files:
            manifest["datasets"][group_name].append({
                "file": f.name,
                "size_mb": round(file_size_mb(f), 2),
                "sha256": sha256(f),
            })

    manifest_path = ZENODO_DIR / "manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
    print(f"\nManifest: {manifest_path}")
    return manifest_path


def create_zenodo_readme(groups, archives):
    """Create README for Zenodo deposit."""
    readme = ZENODO_DIR / "README.md"

    # Build archive table
    rows = []
    for name, path in sorted(archives.items()):
        if path and path.exists():
            mb = file_size_mb(path)
            n = len(groups.get(name.replace(".tar.gz", ""), []))
            rows.append(f"| {name}.tar.gz | {n} files | {mb:.0f} MB |")

    archive_table = "\n".join(rows)

    content = f"""# CrossCell Benchmark Data and Results

## Description

Benchmark datasets and pre-computed results for CrossCell,
a high-performance single-cell data format converter (Seurat RDS ↔ AnnData H5AD).

Total: 65 datasets (42 RDS + 13 H5AD + 10 CELLxGENE benchmark H5AD)

## Archives

Data is split into separate archives for convenient downloading:

| Archive | Contents | Size |
|---------|----------|------|
| seurat_v4_datasets.tar.gz | Seurat V4 RDS (raw + processed) | ~1.3 GB |
| seurat_v5_datasets.tar.gz | Seurat V5 RDS (raw + processed) | ~0.7 GB |
| h5ad_small_datasets.tar.gz | scanpy + scvelo H5AD | ~0.1 GB |
| h5ad_spatial_datasets.tar.gz | squidpy spatial H5AD | ~0.7 GB |
| h5ad_cellxgene_datasets.tar.gz | CELLxGENE small H5AD (3 files) | ~0.7 GB |
| benchmark_results.tar.gz | Pre-computed results + scripts | ~1 MB |

All archives extract to `data/generated/` or `benchmark/` relative paths.

## CELLxGENE Benchmark Datasets (Not Included)

10 large-scale CELLxGENE datasets (10k–836k cells) used for scalability benchmarking
are NOT included in this deposit due to their size (total ~25 GB).

Download them directly from CELLxGENE Discover using the URLs in `cellxgene_sources.json`,
or use the provided download script:

```bash
# Download all 10 CELLxGENE benchmark datasets
python3 data/05_download_h5ad_datasets.py --source cellxgene-benchmark
```

Note: `cellxgene_retina_244k` was excluded from the benchmark due to an abnormal
dense matrix structure (18 GB file requiring 62 GB RAM for inspection).

## Quick Start

```bash
# Clone the repository
git clone https://github.com/TODO/crosscell
cd crosscell

# Download and extract only what you need
# Example: just Seurat V4 + small H5AD for a quick test
tar xzf seurat_v4_datasets.tar.gz
tar xzf h5ad_small_datasets.tar.gz

# Download CELLxGENE benchmark datasets (optional, ~25 GB)
python3 data/05_download_h5ad_datasets.py --source cellxgene-benchmark

# Or extract everything
for f in *.tar.gz; do tar xzf "$f"; done

# Build and run benchmark
docker-compose run --rm dev cargo build --release
cp target/release/crosscell benchmark/crosscell
docker-compose -f benchmark/docker-compose.benchmark.yml build
docker-compose -f benchmark/docker-compose.benchmark.yml up notebook
```

## Tools Compared

| Tool | Version | Language |
|------|---------|----------|
| CrossCell | 0.1.0 | Rust |
| Zellkonverter | 1.20.1 | R (Bioconductor) |
| anndataR | 1.0.0 | R (scverse) |
| convert2anndata | 0.2.0 | R |
| easySCF | 0.2.2.7 (R) / 0.2.5 (Python) | R + Python |

## Data Sources

| Source | Datasets | Reference |
|--------|----------|-----------|
| SeuratData | 42 RDS (11 datasets × V4/V5 × raw/processed) | Hao et al. Cell 2024 |
| CELLxGENE | 3 small + 10 benchmark H5AD | CZ CELLxGENE Discover |
| scanpy | 2 H5AD | Wolf et al. Genome Biology 2018 |
| scvelo | 2 H5AD | Bergen et al. Nature Biotechnology 2020 |
| squidpy | 6 H5AD | Palla et al. Nature Methods 2022 |

## License

Benchmark scripts: MIT OR Apache-2.0
Test datasets: Subject to original data licenses (see individual sources)
"""
    with open(readme, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Zenodo README: {readme}")


def main():
    print("=" * 60)
    print("  CrossCell Benchmark — Zenodo Packaging")
    print("=" * 60)

    ZENODO_DIR.mkdir(parents=True, exist_ok=True)

    # Classify data files
    groups = classify_data_files()
    total_files = sum(len(v) for v in groups.values())
    print(f"\nData files: {total_files}")
    for name, files in groups.items():
        mb = sum(file_size_mb(f) for f in files) if files else 0
        print(f"  {name}: {len(files)} files ({mb:.0f} MB)")

    # Create data archives
    print("\nCreating archives...")
    archives = {}
    for name, files in groups.items():
        archive = create_archive(name, files)
        if archive:
            archives[name] = archive

    # Create results archive
    results_archive = create_results_archive()
    if results_archive:
        archives["benchmark_results"] = results_archive

    # Create CELLxGENE sources JSON (download URLs for large datasets)
    create_cellxgene_sources_json()

    # Create manifest and README
    create_manifest(groups, archives)
    create_zenodo_readme(groups, archives)

    # Summary
    print("\n" + "=" * 60)
    print("  Done! Files ready for Zenodo upload:")
    print("=" * 60)
    total_archive_mb = 0
    for f in sorted(ZENODO_DIR.iterdir()):
        mb = file_size_mb(f)
        total_archive_mb += mb
        print(f"  {f.name:45s} {mb:8.1f} MB")
    print(f"  {'TOTAL':45s} {total_archive_mb:8.1f} MB")

    print(f"\nUpload to: https://zenodo.org/deposit/new")
    print("Tips:")
    print("  1. Upload all .tar.gz + manifest.json + README.md")
    print("  2. License: CC-BY-4.0 (data), MIT (code)")
    print("  3. Keywords: single-cell, benchmark, format conversion")
    print("  4. Link to GitHub repository")


if __name__ == "__main__":
    main()
