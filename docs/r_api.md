# R API Reference

## Installation

```r
install.packages("path/to/crosscell-r", repos = NULL, type = "source")
library(crosscell)
```

## Functions

### read_h5ad_as_seurat

Read an H5AD file and return a Seurat object.

```r
read_h5ad_as_seurat(path, normalize = FALSE, top_genes = NULL, gene_id_column = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | character | required | Path to .h5ad file |
| `normalize` | logical | FALSE | Library-size normalization + log1p |
| `top_genes` | integer or NULL | NULL | Select top N variable genes by variance |
| `gene_id_column` | character or NULL | NULL | Use this var column as gene identifiers |

Returns: Seurat object

```r
# Basic
seurat_obj <- read_h5ad_as_seurat("data.h5ad")

# AI-Ready
seurat_obj <- read_h5ad_as_seurat(
  "data.h5ad",
  normalize = TRUE,
  top_genes = 2000L,
  gene_id_column = "ensembl_id"
)
```

---

### read_h5ad_as_sce

Read an H5AD file and return a SingleCellExperiment object.

```r
read_h5ad_as_sce(path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | character | Path to .h5ad file |

Returns: SingleCellExperiment object

```r
sce_obj <- read_h5ad_as_sce("data.h5ad")
```

---

### read_rds_fast

Fast read of a Seurat RDS file using Rust-native parser.

```r
read_rds_fast(path, normalize = FALSE, top_genes = NULL,
              gene_id_column = NULL, keep_layers = FALSE)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | character | required | Path to .rds file |
| `normalize` | logical | FALSE | Library-size normalization + log1p |
| `top_genes` | integer or NULL | NULL | Select top N variable genes by variance |
| `gene_id_column` | character or NULL | NULL | Use this var column as gene identifiers |
| `keep_layers` | logical | FALSE | Preserve Seurat V5 split layers |

Returns: Seurat object

```r
# Basic
seurat_obj <- read_rds_fast("seurat.rds")

# AI-Ready
seurat_obj <- read_rds_fast(
  "seurat.rds",
  normalize = TRUE,
  top_genes = 2000L,
  gene_id_column = "ensembl_id"
)

# Keep split layers
seurat_obj <- read_rds_fast("integrated.rds", keep_layers = TRUE)
```

---

### write_as_h5ad

Write a Seurat or SingleCellExperiment object to H5AD format.

```r
write_as_h5ad(obj, path, normalize = FALSE, top_genes = NULL, gene_id_column = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | Seurat or SCE | required | Object to write |
| `path` | character | required | Output .h5ad file path |
| `normalize` | logical | FALSE | Library-size normalization + log1p |
| `top_genes` | integer or NULL | NULL | Select top N variable genes by variance |
| `gene_id_column` | character or NULL | NULL | Use this var column as gene identifiers |

```r
# Basic
write_as_h5ad(seurat_obj, "output.h5ad")

# With transforms
write_as_h5ad(seurat_obj, "ai_ready.h5ad",
              normalize = TRUE, top_genes = 2000L)
```

---

### write_rds_fast

Fast write of a Seurat object to RDS format.

```r
write_rds_fast(obj, path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | Seurat | Seurat object to write |
| `path` | character | Output .rds file path |

```r
write_rds_fast(seurat_obj, "output.rds")
```

---

### inspect_file

Return file metadata as a list.

```r
inspect_file(path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | character | Path to .h5ad or .rds file |

Returns: list with `format`, `n_cells`, `n_genes`, `obs_columns`, `var_columns`, `has_embeddings`, `embedding_names`, `has_layers`, `layer_names`.

```r
info <- inspect_file("data.h5ad")
cat(sprintf("%d cells × %d genes\n", info$n_cells, info$n_genes))
```

---

### validate_conversion

Compare two files and return a fidelity report.

```r
validate_conversion(original, converted, tolerance = 1e-7, cluster_column = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `original` | character | required | Path to original file (.h5ad or .rds) |
| `converted` | character | required | Path to converted file (.h5ad or .rds) |
| `tolerance` | numeric | 1e-7 | Numerical tolerance |
| `cluster_column` | character or NULL | NULL | obs column for cluster labels (enables ARI/NMI) |

Returns: list with `passed` (logical), `summary` (character), `detailed_report` (character), and optionally `ari` (numeric), `nmi` (numeric).

```r
result <- validate_conversion("original.h5ad", "roundtrip.h5ad")
if (result$passed) message("Validation passed!")

# With cluster metrics
result <- validate_conversion(
  "original.h5ad", "roundtrip.h5ad",
  cluster_column = "cell_type"
)
cat(sprintf("ARI: %.4f, NMI: %.4f\n", result$ari, result$nmi))
```

## Typical Workflows

### H5AD → Seurat Analysis

```r
library(crosscell)
library(Seurat)

seurat_obj <- read_h5ad_as_seurat("scanpy_processed.h5ad")
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj, reduction = "umap")
write_as_h5ad(seurat_obj, "seurat_analyzed.h5ad")
```

### AI-Ready Export

```r
library(crosscell)

# One-step: Seurat → normalized H5AD with HVG
write_as_h5ad(
  seurat_obj, "ai_ready.h5ad",
  normalize = TRUE,
  top_genes = 2000L,
  gene_id_column = "ensembl_id"
)
```

### Roundtrip Validation

```r
library(crosscell)

# Convert and validate
write_as_h5ad(seurat_obj, "converted.h5ad")
result <- validate_conversion("original.rds", "converted.h5ad")
stopifnot(result$passed)
```
