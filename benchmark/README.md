# CrossCell Benchmark: Reproducibility Guide

## Overview

This benchmark compares CrossCell against 4 existing single-cell data conversion tools:
- **Zellkonverter** (Bioconductor, R/Python via basilisk)
- **anndataR** (scverse, pure R)
- **convert2anndata** (R/Python)
- **easySCF** (R/Python, custom HDF5 format)

All benchmarks run inside a Docker container for full reproducibility.

## Quick Verification (< 5 minutes)

Pre-computed results are available in `benchmark/results/`:

```bash
# View pass rates for all tools
python3 benchmark/scripts/check_results.py

# View detailed comparison
python3 benchmark/scripts/compare_tools.py
```

## Full Reproduction

### Prerequisites

- Docker and Docker Compose
- ~10 GB disk space (Docker image + test data)
- ~16 GB RAM recommended

### Step 1: Obtain Test Data

Test datasets are archived on Zenodo (DOI: [10.5281/zenodo.18610876](https://doi.org/10.5281/zenodo.18610876)).

```bash
# Option A: Download from Zenodo (split archives, download what you need)
# See data/README.md for full Zenodo upload/download instructions
wget https://zenodo.org/records/18610876/files/seurat_v4_datasets.tar.gz
wget https://zenodo.org/records/18610876/files/seurat_v5_datasets.tar.gz
wget https://zenodo.org/records/18610876/files/h5ad_small_datasets.tar.gz
wget https://zenodo.org/records/18610876/files/h5ad_spatial_datasets.tar.gz
wget https://zenodo.org/records/18610876/files/h5ad_cellxgene_datasets.tar.gz
for f in *.tar.gz; do tar xzf "$f"; done

# Option B: Generate from SeuratData (requires R environment)
# See data/README.md for instructions
```

**Dataset summary**: 42 Seurat RDS files (V4/V5 × raw/processed × 11 datasets from SeuratData) + 13 H5AD files (CELLxGENE, scanpy, scvelo, squidpy). Cell counts range from 1,438 to 73,655. Full details in `benchmark/config/datasets.json`.

### Step 2: Build CrossCell

```bash
# Build CrossCell binary in dev container
docker-compose run --rm dev cargo build --release

# Copy binary to benchmark directory
cp target/release/crosscell benchmark/crosscell
```

### Step 3: Build Benchmark Docker Image

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml build
```

This installs all 5 tools in a single Docker image. Tool versions are recorded in `benchmark/results/tool_versions.json`.

### Step 4: Run Benchmarks

```bash
# Run all tools, all datasets, 3 runs each (~2-3 hours)
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  bash /benchmark/scripts/run_all.sh

# Or run a specific tool
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  bash /benchmark/scripts/run_all.sh --tool crosscell

# Or run a specific direction
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  bash /benchmark/scripts/run_all.sh --tool zellkonverter --direction rds_to_h5ad
```

### Step 5: Reproduce easySCF Round-trip Test

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  python3 /benchmark/scripts/test_easyscf_roundtrip.py --skip-download
```

This tests three conversion pathways:
- **Test A** (R → H5 → Python): Uses `easySCFr::saveH5()` then `easySCFpy.loadH5()`
- **Test B** (Python → H5 → R): Uses `easySCFpy.saveH5()` then `easySCFr::loadH5()`
- **Test C** (R → H5 → Python → H5 → R): Full round-trip

Results are saved to `benchmark/results/easyscf_roundtrip.json`.

## Results Summary

### Conversion Pass Rates (3 runs each)

| Tool | rds_to_h5ad | h5ad_to_rds | Total |
|------|-------------|-------------|-------|
| CrossCell | 42/42 (100%) | 13/13 (100%) | 55/55 (100%) |
| Zellkonverter | 41/42 (98%) | 13/13 (100%) | 54/55 (98%) |
| anndataR | 41/42 (98%) | 3/13 (23%) | 44/55 (80%) |
| convert2anndata | 19/42 (45%) | N/A | 19/42 (45%) |
| easySCF | 41/42 (98%) | N/A | 41/42 (98%) |

### easySCF Round-trip Results

| Pathway | Pass Rate | Root Cause of Failure |
|---------|-----------|----------------------|
| R → H5 → Python | 5/5 (100%) | — |
| Python → H5 → R | 0/3 (0%) | Encoding incompatibility (see below) |
| R → H5 → Py → H5 → R | 0/5 (0%) | Fails at Python-to-R leg |

**Root cause**: `easySCFpy.saveH5()` uses `anndata._io.specs.write_elem()` to serialize var DataFrames with anndata's categorical encoding, but `easySCFr::loadH5()` expects `h5[["var/rawvar/_index"]]` as a simple HDF5 string dataset. Additionally, `loadH5()` has an error-handling bug where `tryCatch` swallows the real error, causing the misleading `object 'sce' not found` message.

### Performance (CrossCell vs others, rds_to_h5ad median)

| Metric | vs Zellkonverter | vs anndataR | vs convert2anndata | vs easySCF |
|--------|-----------------|-------------|-------------------|-----------|
| Speed | 7.2x faster | 5.8x faster | 10x+ faster | 5-24x faster |
| Memory | 2.5x less | 2.1x less | 3-9x less | 2-5x less |

## File Structure

```
benchmark/
├── Dockerfile.benchmark          # Complete environment definition
├── docker-compose.benchmark.yml  # Docker Compose config
├── crosscell                     # Pre-built CrossCell binary (Linux x86_64)
├── scripts/
│   ├── run_benchmark.py          # Main benchmark script (5 tools × 2 directions × 3 runs)
│   ├── run_all.sh                # One-click entry point
│   ├── test_easyscf_roundtrip.py # easySCF round-trip test
│   ├── validate_output.py        # Output validation
│   ├── check_results.py          # Quick results viewer
│   └── compare_tools.py          # Detailed comparison
├── data/
│   ├── benchmark_testcases.json   # Benchmark test cases (tool × dataset × direction)
│   └── datasets.json             # Dataset metadata (cells, genes, format, source)
└── results/
    ├── crosscell.json            # CrossCell results (3 runs)
    ├── zellkonverter.json        # Zellkonverter results (3 runs)
    ├── anndataR.json             # anndataR results (3 runs)
    ├── convert2anndata.json      # convert2anndata results (3 runs)
    ├── easySCF.json              # easySCF results (3 runs)
    ├── easyscf_roundtrip.json    # easySCF round-trip test results
    ├── tool_versions.json        # Software versions
    ├── dataset_summary.json      # Dataset scale summary
    └── validation.json           # Output validation results
```

## Software Versions

| Software | Version |
|----------|---------|
| CrossCell | 0.1.0 |
| Zellkonverter | 1.20.1 |
| anndataR | 1.0.0 |
| convert2anndata | 0.2.0 |
| easySCFr | 0.2.2.7 |
| easySCFpy | 0.2.5 |
| Seurat | 5.4.0 |
| anndata (Python) | 0.12.9 |
| R | 4.5.2 |
| Python | 3.12.12 |

## Test Datasets

All datasets are publicly available:

| Source | Datasets | Reference |
|--------|----------|-----------|
| SeuratData | pbmc3k, cbmc, ifnb, panc8, bmcite, pbmcsca, hcabm40k, ssHippo, stxKidney, celegans.embryo, thp1.eccite | Hao et al. 2024 |
| CELLxGENE | pbmc_15k, heart_23k, brain_40k | CZ CELLxGENE Discover |
| scanpy | pbmc3k, pbmc3k_processed | Wolf et al. 2018 |
| scvelo | pancreas, dentategyrus | Bergen et al. 2020 |
| squidpy | imc, seqfish, visium_hne, merfish, slideseqv2, mibitof | Palla et al. 2022 |

## License

Benchmark scripts and configurations are released under the same license as CrossCell (MIT OR Apache-2.0).
