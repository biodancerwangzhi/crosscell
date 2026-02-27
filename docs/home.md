# CrossCell

High-performance, zero-dependency single-cell data format converter.

## Overview

CrossCell converts between **AnnData (.h5ad)** and **Seurat (.rds)** formats with three interfaces:

- **CLI** — Standalone binary, no runtime dependencies
- **Python** — Native AnnData integration via PyO3
- **R** — Native Seurat/SCE integration via extendr

All three interfaces share the same Rust core engine.

## Key Features

| Feature | Description |
|---------|-------------|
| Zero dependency | No R or Python runtime required |
| High fidelity | Pearson R = 1.0, MSE < 10⁻¹⁴ in roundtrip |
| Seurat V5 | Auto-merge Split Layers, optional `--keep-layers` |
| AI-Ready | Built-in normalize, HVG selection, gene ID mapping |
| CELLxGENE | Schema 5.x validation with precise type mapping |
| Streaming | Constant-memory mode for 1M+ cell datasets |
| Validation | Roundtrip fidelity audit (9 metrics + ARI/NMI) |

## Quick Example

```bash
# CLI
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata --normalize --top-genes 2000
```

```python
# Python
import crosscell
adata = crosscell.read_rds("seurat.rds", normalize=True, top_genes=2000)
```

```r
# R
library(crosscell)
seurat_obj <- read_h5ad_as_seurat("data.h5ad", normalize = TRUE, top_genes = 2000L)
```

## Supported Data Components

| Component | AnnData | Seurat | Status |
|-----------|---------|--------|--------|
| Expression Matrix | X | counts/data | ✅ |
| Cell Metadata | obs | meta.data | ✅ |
| Gene Metadata | var | meta.features | ✅ |
| PCA/UMAP/tSNE | obsm | reductions | ✅ |
| Layers | layers | layers | ✅ |
| Categorical Variables | categorical | factors | ✅ |
| Spatial Coordinates | obsm['spatial'] | images | ✅ |
| Seurat V5 Assay5 | — | Assay5 | ✅ |
| Split Layers | — | counts.1/2/… | ✅ |

## Documentation

| Page | Description |
|------|-------------|
| [Installation](installation.md) | CLI, Python, R installation |
| [Quick Start](quick_start.md) | First conversion in 2 minutes |
| [CLI Reference](cli_reference.md) | All commands and parameters |
| [Python API](python_api.md) | Python package reference |
| [R API](r_api.md) | R package reference |
| [Advanced Usage](advanced_usage.md) | Streaming, batch, AI-Ready, CELLxGENE |
| [FAQ](faq.md) | Common questions |
| [Contributing](contributing.md) | Development guide |

## Installation

```bash
# CLI
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-linux-x64.tar.gz | tar xz

# Python
pip install crosscell

# R
install.packages("crosscell", repos = "https://biodancerwangzhi.r-universe.dev")
```

## License

MIT OR Apache-2.0
