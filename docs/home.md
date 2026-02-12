# CrossCell

High-performance single-cell transcriptomics data format converter.

## Overview

CrossCell enables bidirectional conversion between **AnnData (.h5ad)** and **Seurat (.rds)** formats with:

- **Pure Rust implementation** - No Python/R runtime dependencies
- **High performance** - Parallel processing, lazy loading, memory optimization
- **Streaming support** - Handle ultra-large datasets (10M+ cells) with constant memory
- **Data integrity** - Preserves expression matrices, metadata, embeddings, spatial data
- **CLI tool** - Standalone executable, easy deployment

## Supported Data

| Component | AnnData | Seurat | Status |
|-----------|---------|--------|--------|
| Expression Matrix | X | counts/data | ✅ |
| Cell Metadata | obs | meta.data | ✅ |
| Gene Metadata | var | meta.features | ✅ |
| PCA/UMAP/tSNE | obsm | reductions | ✅ |
| Layers | layers | layers | ✅ |
| Categorical | categorical | factors | ✅ |
| Spatial Coords | obsm['spatial'] | images | ✅ |
| Seurat V5 Assay5 | - | Assay5 | ✅ |

## Quick Example

```bash
# AnnData → Seurat
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# Seurat → AnnData  
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata

# Streaming mode for large files
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming

# Inspect file
crosscell inspect -i pbmc3k.h5ad
```

## Processing Modes

| Mode | Use Case | Memory Usage |
|------|----------|--------------|
| Standard | < 50k cells | Full dataset in memory |
| Lazy Loading | 50k - 1M cells | Reduced memory |
| Streaming | 1M+ cells | Constant (~2-5 GB) |
| Disk-backed | 10M+ cells | Minimal memory |

## Documentation

- [Installation](installation.md)
- [Quick Start](quick_start.md)
- [Advanced Usage](advanced_usage.md)
- [FAQ](faq.md)
- [Contributing](contributing.md)

## License

MIT OR Apache-2.0
