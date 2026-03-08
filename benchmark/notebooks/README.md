# Benchmark Notebooks: Reproducibility Guide

This directory contains all Jupyter notebooks used to generate the figures and compute the metrics in the manuscript. Follow the steps below to reproduce every result from scratch.

## Prerequisites

- Docker and Docker Compose installed
- ~10 GB disk space (Docker image + test data)
- ~16 GB RAM for standard datasets; 64 GB for CELLxGENE large-scale datasets

## Step 1: Pull the CrossCell benchmark docker image


```bash
docker pull biodancer/crosscell-benchmark:latest
# rename the image
docker tag biodancer/crosscell-benchmark crosscell-benchmark:latest

```
This image `crosscell-benchmark:latest`, which contains all five tools (CrossCell, Zellkonverter, anndataR, convert2anndata, easySCF), R 4.5, Python 3.12, and Jupyter.

## Step 2: Clone the CrossCell repo from github 
```bash
git clone https://github.com/biodancerwangzhi/crosscell
# go to the repository root:
cd crosscell
```

## Step 3: Download test datasets

Datasets are archived on Zenodo: [10.5281/zenodo.18610876](https://doi.org/10.5281/zenodo.18610876)

```bash
# Download split archives
wget https://zenodo.org/records/18610877/files/seurat_v4_datasets.tar.gz
wget https://zenodo.org/records/18610877/files/seurat_v5_datasets.tar.gz
wget https://zenodo.org/records/18610877/files/h5ad_small_datasets.tar.gz
wget https://zenodo.org/records/18610877/files/h5ad_spatial_datasets.tar.gz
wget https://zenodo.org/records/18610877/files/h5ad_cellxgene_datasets.tar.gz

# Extract into data/ directory
for f in *.tar.gz; do tar xzf "$f"; done
```
Alternatively, generate RDS datasets from SeuratData (see `data/README.md`).

## Step 4: Launch Jupyter

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml up notebook
```
Open your browser at `http://localhost:8888`. No token or password is required.


## Step 5: Run notebooks

### Recommended execution order

Notebooks are organized into two groups: **compute** notebooks that generate result JSON files, and **figure** notebooks that read those results and produce manuscript figures.

#### 1. Environment check (optional)

| Notebook | Description |
|----------|-------------|
| `environment_check.ipynb` | Verify all tools are installed and datasets are accessible |

#### 2. Compute notebooks (generate results)

Run these first. They read raw data and write JSON files to `benchmark/results/`.

| Notebook | Description | Output |
|----------|-------------|--------|
| `compute_robustness.ipynb` | Run all 5 tools × 55 datasets × 3 runs | `robustness_benchmark.json` |
| `compute_fidelity.ipynb` | Expression fidelity (Pearson r, MSE) | `tools_fidelity.json` |
| `compute_type_fidelity.ipynb` | Metadata type preservation (Int32, Bool, etc.) | `type_fidelity.json` |
| `compute_roundtrip_fidelity.ipynb` | Bidirectional roundtrip exact match | `roundtrip_fidelity.json` |
| `compute_scalability.ipynb` | CELLxGENE large-scale datasets (CrossCell only) | `cellxgene_scalability.json` |

#### 3. Figure notebooks (generate figures)

These read from `benchmark/results/` and write figures to `benchmark/figures/main/`.

| Notebook | Figure | Description |
|----------|--------|-------------|
| `fig1_design.ipynb` | Figure 1 | Benchmark design and evaluation framework |
| `fig2_robustness.ipynb` | Figure 2 | Conversion robustness and failure modes |
| `fig3_performance.ipynb` | Figure 3 | Runtime performance and memory efficiency |
| `fig4_scalability.ipynb` | Figure 4 | Large-scale scalability (CELLxGENE) |
| `fig5_data_fidelity.ipynb` | Figure 5 | Expression and type preservation |
| `fig6_split_layers.ipynb` | Figure 6 | Seurat V5 Split Layers compatibility |
| `fig7_ecosystem.ipynb` | Figure 7 | CrossCell as ecosystem infrastructure |

#### 4. Supplementary figure notebooks

| Notebook | Description |
|----------|-------------|
| `supp_fig1_performance_detail.ipynb` | Detailed per-dataset performance breakdown |
| `supp_fig2_memory_detail.ipynb` | Detailed memory scaling analysis |
| `supp_fig3_cellxgene.ipynb` | CELLxGENE dataset details |
| `supp_fig4_aiready.ipynb` | AI-Ready workflow validation |

## Quick start: figures only (using pre-computed results)

If you just want to regenerate figures from the pre-computed results already in `benchmark/results/`, skip the compute notebooks and run only the figure notebooks (fig1–fig7, supp_fig1–supp_fig4).

```bash
# Launch Jupyter
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm -p 8888:8888 notebook

# Then open fig1_design.ipynb, fig2_robustness.ipynb, ... in the browser
```

## Utility

| File | Description |
|------|-------------|
| `benchmark_tutorial.ipynb` | Interactive tutorial for using CrossCell |
| `_gen_notebooks.py` | Script used to generate notebook templates |

## Troubleshooting

- **"No such file" errors in compute notebooks**: Ensure datasets are extracted into `data/` (Step 3).
- **Out of memory**: CELLxGENE large-scale tests require 64 GB RAM. Standard tests work with 16 GB.
- **Port 8888 already in use**: Change the port mapping, e.g., `-p 9999:8888`, then open `http://localhost:9999`.
- **Figures look different**: Ensure you are using the same Docker image version. Matplotlib rendering may vary slightly across platforms.
