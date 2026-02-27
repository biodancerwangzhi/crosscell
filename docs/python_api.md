# Python API Reference

## Installation

```bash
pip install maturin
cd crosscell-py
maturin develop --release
```

```python
import crosscell
print(crosscell.__version__)
```

## Functions

### read_rds

Read a Seurat RDS file and return an AnnData object.

```python
crosscell.read_rds(
    path,
    *,
    normalize=False,
    top_genes=None,
    gene_id_column=None,
    keep_layers=False,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | str | required | Path to .rds file |
| `normalize` | bool | False | Library-size normalization + log1p |
| `top_genes` | int or None | None | Select top N variable genes by variance |
| `gene_id_column` | str or None | None | Use this var column as gene identifiers |
| `keep_layers` | bool | False | Preserve Seurat V5 split layers as separate AnnData layers |

Returns: `anndata.AnnData`

```python
# Basic
adata = crosscell.read_rds("seurat.rds")

# AI-Ready pipeline
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# Keep split layers
adata = crosscell.read_rds("integrated.rds", keep_layers=True)
print(adata.layers.keys())  # dict_keys(['counts.1', 'counts.2'])
```

---

### read_h5ad

Read an H5AD file and return an AnnData object.

```python
crosscell.read_h5ad(
    path,
    *,
    normalize=False,
    top_genes=None,
    gene_id_column=None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | str | required | Path to .h5ad file |
| `normalize` | bool | False | Library-size normalization + log1p |
| `top_genes` | int or None | None | Select top N variable genes by variance |
| `gene_id_column` | str or None | None | Use this var column as gene identifiers |

Returns: `anndata.AnnData`

```python
adata = crosscell.read_h5ad("data.h5ad", normalize=True, top_genes=2000)
```

---

### write_h5ad

Write an AnnData object to H5AD format.

```python
crosscell.write_h5ad(adata, path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `adata` | anndata.AnnData | AnnData object to write |
| `path` | str | Output .h5ad file path |

```python
crosscell.write_h5ad(adata, "output.h5ad")
```

---

### write_rds

Write an AnnData object to Seurat RDS format.

```python
crosscell.write_rds(adata, path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `adata` | anndata.AnnData | AnnData object (must be non-empty) |
| `path` | str | Output .rds file path |

```python
crosscell.write_rds(adata, "output.rds")
```

---

### inspect

Return file metadata as a dictionary.

```python
crosscell.inspect(path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | str | Path to .h5ad or .rds file |

Returns: `dict` with keys like `format`, `n_cells`, `n_genes`, `is_sparse`, `embedding_names`, `layer_names`, etc.

```python
info = crosscell.inspect("data.h5ad")
print(f"{info['n_cells']} cells × {info['n_genes']} genes")
```

---

### validate

Compare two files and return a fidelity report.

```python
crosscell.validate(
    original,
    converted,
    *,
    tolerance=1e-7,
    cluster_column=None,
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `original` | str | required | Path to original file (.h5ad or .rds) |
| `converted` | str | required | Path to converted file (.h5ad or .rds) |
| `tolerance` | float | 1e-7 | Numerical tolerance |
| `cluster_column` | str or None | None | obs column for cluster labels (enables ARI/NMI) |

Returns: `dict` with keys `passed` (bool), `summary` (str), `detailed_report` (str), and optionally `ari` (float), `nmi` (float).

```python
result = crosscell.validate("original.h5ad", "roundtrip.h5ad")
print(result["passed"])   # True
print(result["summary"])

# With cluster metrics
result = crosscell.validate(
    "original.h5ad",
    "roundtrip.h5ad",
    cluster_column="cell_type",
)
print(f"ARI: {result['ari']:.4f}, NMI: {result['nmi']:.4f}")
```

## Typical Workflows

### Seurat → Scanpy

```python
import crosscell
import scanpy as sc

adata = crosscell.read_rds("seurat.rds")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
crosscell.write_h5ad(adata, "analyzed.h5ad")
```

### AI-Ready Pipeline (Geneformer)

```python
import crosscell
import numpy as np

# One-step: RDS → normalized H5AD with top 2000 genes
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# Per-cell rank encoding (~15 lines)
for i in range(adata.n_obs):
    row = adata.X[i].toarray().flatten() if hasattr(adata.X[i], 'toarray') else adata.X[i]
    nonzero = np.nonzero(row)[0]
    ranked = nonzero[np.argsort(-row[nonzero])][:2048]
    # ranked contains gene indices sorted by expression (descending)
```

### Roundtrip Validation

```python
import crosscell

crosscell.write_h5ad(crosscell.read_rds("input.rds"), "converted.h5ad")
adata2 = crosscell.read_h5ad("converted.h5ad")
crosscell.write_rds(adata2, "roundtrip.rds")

result = crosscell.validate("input.rds", "roundtrip.rds", cluster_column="seurat_clusters")
assert result["passed"]
```
