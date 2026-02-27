# Python API 参考

## 安装

```bash
pip install maturin
cd crosscell-py
maturin develop --release
```

```python
import crosscell
print(crosscell.__version__)
```

## 函数

### read_rds

读取 Seurat RDS 文件，返回 AnnData 对象。

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

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `path` | str | 必需 | .rds 文件路径 |
| `normalize` | bool | False | 库大小归一化 + log1p |
| `top_genes` | int 或 None | None | 按方差选择前 N 个高变基因 |
| `gene_id_column` | str 或 None | None | 用指定 var 列作为基因标识符 |
| `keep_layers` | bool | False | 保留 Seurat V5 split layers 为独立的 AnnData layers |

返回：`anndata.AnnData`

```python
# 基础用法
adata = crosscell.read_rds("seurat.rds")

# AI-Ready 管道
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# 保留 split layers
adata = crosscell.read_rds("integrated.rds", keep_layers=True)
print(adata.layers.keys())  # dict_keys(['counts.1', 'counts.2'])
```

---

### read_h5ad

读取 H5AD 文件，返回 AnnData 对象。

```python
crosscell.read_h5ad(
    path,
    *,
    normalize=False,
    top_genes=None,
    gene_id_column=None,
)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `path` | str | 必需 | .h5ad 文件路径 |
| `normalize` | bool | False | 库大小归一化 + log1p |
| `top_genes` | int 或 None | None | 按方差选择前 N 个高变基因 |
| `gene_id_column` | str 或 None | None | 用指定 var 列作为基因标识符 |

返回：`anndata.AnnData`

```python
adata = crosscell.read_h5ad("data.h5ad", normalize=True, top_genes=2000)
```

---

### write_h5ad

将 AnnData 对象写入 H5AD 格式。

```python
crosscell.write_h5ad(adata, path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `adata` | anndata.AnnData | 要写入的 AnnData 对象 |
| `path` | str | 输出 .h5ad 文件路径 |

---

### write_rds

将 AnnData 对象写入 Seurat RDS 格式。

```python
crosscell.write_rds(adata, path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `adata` | anndata.AnnData | AnnData 对象（不能为空） |
| `path` | str | 输出 .rds 文件路径 |

---

### inspect

返回文件元数据字典。

```python
crosscell.inspect(path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `path` | str | .h5ad 或 .rds 文件路径 |

返回：`dict`，包含 `format`、`n_cells`、`n_genes`、`is_sparse`、`embedding_names`、`layer_names` 等。

```python
info = crosscell.inspect("data.h5ad")
print(f"{info['n_cells']} cells × {info['n_genes']} genes")
```

---

### validate

比较两个文件，返回保真度报告。

```python
crosscell.validate(
    original,
    converted,
    *,
    tolerance=1e-7,
    cluster_column=None,
)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `original` | str | 必需 | 原始文件路径 |
| `converted` | str | 必需 | 转换后文件路径 |
| `tolerance` | float | 1e-7 | 数值容差 |
| `cluster_column` | str 或 None | None | obs 中的聚类标签列（启用 ARI/NMI） |

返回：`dict`，包含 `passed`（bool）、`summary`（str）、`detailed_report`（str），以及可选的 `ari`（float）、`nmi`（float）。

```python
result = crosscell.validate("original.h5ad", "roundtrip.h5ad")
print(result["passed"])  # True

# 带聚类指标
result = crosscell.validate(
    "original.h5ad", "roundtrip.h5ad",
    cluster_column="cell_type",
)
print(f"ARI: {result['ari']:.4f}, NMI: {result['nmi']:.4f}")
```

## 典型工作流

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

### AI-Ready 管道（Geneformer）

```python
import crosscell
import numpy as np

# 一步完成：RDS → 归一化 H5AD + 前 2000 高变基因
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# Per-cell rank encoding（约 15 行代码）
for i in range(adata.n_obs):
    row = adata.X[i].toarray().flatten() if hasattr(adata.X[i], 'toarray') else adata.X[i]
    nonzero = np.nonzero(row)[0]
    ranked = nonzero[np.argsort(-row[nonzero])][:2048]
```

### 往返验证

```python
import crosscell

crosscell.write_h5ad(crosscell.read_rds("input.rds"), "converted.h5ad")
adata2 = crosscell.read_h5ad("converted.h5ad")
crosscell.write_rds(adata2, "roundtrip.rds")

result = crosscell.validate("input.rds", "roundtrip.rds", cluster_column="seurat_clusters")
assert result["passed"]
```
