# CrossCell

High-performance bidirectional converter between AnnData (`.h5ad`) and Seurat (`.rds`) for single-cell transcriptomics data, written in Rust.

## Features

- **Bidirectional**: H5AD ↔ RDS full roundtrip support
- **Runtime-independent**: no R or Python runtime required (CLI)
- **Seurat V5**: Split Layers detection, merging, and `--keep-layers` option
- **AI-Ready**: built-in normalize, HVG selection, gene ID mapping
- **Scalable**: streaming mode for 1M+ cell datasets
- **Three interfaces**: CLI, Python (PyO3), R (extendr)

## Installation

### CLI (prebuilt binaries)

Download from [GitHub Releases](https://github.com/biodancerwangzhi/crosscell/releases):

```bash
# Linux x86_64
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-linux-x64.tar.gz | tar xz
sudo mv crosscell-linux-x64/crosscell /usr/local/bin/

# macOS Apple Silicon
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-macos-arm64.tar.gz | tar xz
sudo mv crosscell-macos-arm64/crosscell /usr/local/bin/
```

### Python

```bash
pip install crosscell
```

### R

```r
install.packages("crosscell", repos = "https://biodancerwangzhi.r-universe.dev")
```

## Quick Start

### CLI

```bash
# AnnData → Seurat
crosscell convert -i input.h5ad -o output.rds -f seurat

# Seurat → AnnData
crosscell convert -i input.rds -o output.h5ad -f anndata

# AI-Ready export
crosscell convert -i seurat.rds -o ai_ready.h5ad -f anndata \
  --normalize --top-genes 2000 --gene-id-column ensembl_id

# Inspect file
crosscell inspect -i input.h5ad
```

### Python

```python
import crosscell

adata = crosscell.read_rds("seurat.rds")
crosscell.write_h5ad(adata, "output.h5ad")
```

### R

```r
library(crosscell)

seurat_obj <- read_h5ad_as_seurat("data.h5ad")
write_as_h5ad(seurat_obj, "output.h5ad")
```

## Architecture

```
H5AD ──→ [HDF5 Reader] ──→ Apache Arrow IR ──→ [RDS Writer] ──→ RDS
RDS  ──→ [RDS Reader]  ──→ Apache Arrow IR ──→ [H5AD Writer] ──→ H5AD
```

## Documentation

See [`docs/`](docs/) for full documentation:
[Installation](docs/installation.md) · [Quick Start](docs/quick_start.md) · [CLI Reference](docs/cli_reference.md) · [Python API](docs/python_api.md) · [R API](docs/r_api.md) · [Advanced Usage](docs/advanced_usage.md)

## Benchmark

A containerized benchmark suite covering 66 datasets is included in [`benchmark/`](benchmark/).

## Citation

```
[To be added upon publication]
```

## License

[MIT](LICENSE-MIT) OR [Apache-2.0](LICENSE-APACHE)
