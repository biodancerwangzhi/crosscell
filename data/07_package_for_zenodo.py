#!/usr/bin/env python3
"""
Package benchmark data and results for Zenodo upload.

按类别分包，方便用户按需下载：
  1. seurat_v4_datasets.tar.gz      — Seurat V4 RDS 文件 (~1.3 GB)
  2. seurat_v5_datasets.tar.gz      — Seurat V5 RDS 文件 (~0.7 GB)
  3. h5ad_small_datasets.tar.gz     — scanpy + scvelo H5AD (~0.1 GB)
  4. h5ad_spatial_datasets.tar.gz   — squidpy 空间数据 H5AD (~0.7 GB)
  5. h5ad_cellxgene_datasets.tar.gz — CELLxGENE 中大规模 H5AD (~0.7 GB)
  6. benchmark_results.tar.gz       — 预计算结果 + 脚本 (~1 MB)

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
    """Classify data files into groups for separate archives."""
    groups = {
        "seurat_v4_datasets": [],
        "seurat_v5_datasets": [],
        "h5ad_small_datasets": [],      # scanpy + scvelo
        "h5ad_spatial_datasets": [],     # squidpy
        "h5ad_cellxgene_datasets": [],   # cellxgene
    }

    if not DATA_DIR.exists():
        print(f"WARNING: {DATA_DIR} does not exist")
        return groups

    for f in sorted(DATA_DIR.iterdir()):
        name = f.name
        if name.startswith("seurat_v4_") and name.endswith(".rds"):
            groups["seurat_v4_datasets"].append(f)
        elif name.startswith("seurat_v5_") and name.endswith(".rds"):
            groups["seurat_v5_datasets"].append(f)
        elif name.endswith(".h5ad"):
            if name.startswith(("scanpy_", "scvelo_")):
                groups["h5ad_small_datasets"].append(f)
            elif name.startswith("squidpy_"):
                groups["h5ad_spatial_datasets"].append(f)
            elif name.startswith("cellxgene_"):
                groups["h5ad_cellxgene_datasets"].append(f)
            else:
                # 未知来源的 h5ad 放到 small
                groups["h5ad_small_datasets"].append(f)

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

## Archives

Data is split into separate archives for convenient downloading:

| Archive | Contents | Size |
|---------|----------|------|
| seurat_v4_datasets.tar.gz | Seurat V4 RDS (raw + processed) | ~1.3 GB |
| seurat_v5_datasets.tar.gz | Seurat V5 RDS (raw + processed) | ~0.7 GB |
| h5ad_small_datasets.tar.gz | scanpy + scvelo H5AD | ~0.1 GB |
| h5ad_spatial_datasets.tar.gz | squidpy spatial H5AD | ~0.7 GB |
| h5ad_cellxgene_datasets.tar.gz | CELLxGENE medium/large H5AD | ~0.7 GB |
| benchmark_results.tar.gz | Pre-computed results + scripts | ~1 MB |

All archives extract to `data/generated/` or `benchmark/` relative paths.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/TODO/crosscell
cd crosscell

# Download and extract only what you need
# Example: just Seurat V4 + small H5AD for a quick test
tar xzf seurat_v4_datasets.tar.gz
tar xzf h5ad_small_datasets.tar.gz

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

| Source | Reference |
|--------|-----------|
| SeuratData | Hao et al. Cell 2024 |
| CELLxGENE | CZ CELLxGENE Discover |
| scanpy | Wolf et al. Genome Biology 2018 |
| scvelo | Bergen et al. Nature Biotechnology 2020 |
| squidpy | Palla et al. Nature Methods 2022 |

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
