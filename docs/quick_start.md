# Quick Start

## CLI

```bash
# AnnData → Seurat
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# Seurat → AnnData
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata

# With validation
crosscell convert -i input.h5ad -o output.rds -f seurat --validate

# AI-Ready (normalize + HVG + gene ID mapping)
crosscell convert -i seurat.rds -o ai_ready.h5ad -f anndata \
  --normalize --top-genes 2000 --gene-id-column ensembl_id

# Inspect
crosscell inspect -i data.h5ad --detailed

# Validate roundtrip
crosscell validate --original original.h5ad --converted roundtrip.h5ad
```

## Python

```python
import crosscell

# Read Seurat RDS → AnnData
adata = crosscell.read_rds("seurat.rds")

# Read with AI-Ready transforms
adata = crosscell.read_rds("seurat.rds", normalize=True, top_genes=2000)

# Read H5AD
adata = crosscell.read_h5ad("data.h5ad")

# Write
crosscell.write_h5ad(adata, "output.h5ad")
crosscell.write_rds(adata, "output.rds")

# Inspect
info = crosscell.inspect("data.h5ad")
print(f"{info['n_cells']} cells × {info['n_genes']} genes")

# Validate
result = crosscell.validate("original.h5ad", "roundtrip.h5ad")
print(result["passed"])  # True
```

## R

```r
library(crosscell)

# Read H5AD → Seurat
seurat_obj <- read_h5ad_as_seurat("data.h5ad")

# Read with AI-Ready transforms
seurat_obj <- read_h5ad_as_seurat("data.h5ad", normalize = TRUE, top_genes = 2000L)

# Fast read Seurat RDS
seurat_obj <- read_rds_fast("seurat.rds")

# Read with split layers preserved
seurat_obj <- read_rds_fast("integrated.rds", keep_layers = TRUE)

# Write
write_as_h5ad(seurat_obj, "output.h5ad")
write_rds_fast(seurat_obj, "output.rds")

# Inspect
info <- inspect_file("data.h5ad")

# Validate
result <- validate_conversion("original.h5ad", "roundtrip.h5ad")
```

## Cross-Ecosystem Workflow

### R → Python

```r
# R: preprocess
library(Seurat); library(crosscell)
seurat_obj <- readRDS("raw.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
write_as_h5ad(seurat_obj, "processed.h5ad")
```

```python
# Python: continue
import crosscell, scanpy as sc
adata = crosscell.read_h5ad("processed.h5ad")
sc.tl.leiden(adata)
sc.pl.umap(adata, color="leiden")
```

### Python → R

```python
# Python: preprocess
import crosscell, scanpy as sc
adata = crosscell.read_rds("seurat.rds")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
crosscell.write_h5ad(adata, "processed.h5ad")
```

```r
# R: continue
library(crosscell); library(Seurat)
seurat_obj <- read_h5ad_as_seurat("processed.h5ad")
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj)
```

## Next Steps

- [CLI Reference](cli_reference.md) — All commands and parameters
- [Python API](python_api.md) — Full Python reference
- [R API](r_api.md) — Full R reference
- [Advanced Usage](advanced_usage.md) — Streaming, batch, AI-Ready, CELLxGENE
