# R API 参考

## 安装

```r
install.packages("path/to/crosscell-r", repos = NULL, type = "source")
library(crosscell)
```

## 函数

### read_h5ad_as_seurat

读取 H5AD 文件，返回 Seurat 对象。

```r
read_h5ad_as_seurat(path, normalize = FALSE, top_genes = NULL, gene_id_column = NULL)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `path` | character | 必需 | .h5ad 文件路径 |
| `normalize` | logical | FALSE | 库大小归一化 + log1p |
| `top_genes` | integer 或 NULL | NULL | 按方差选择前 N 个高变基因 |
| `gene_id_column` | character 或 NULL | NULL | 用指定 var 列作为基因标识符 |

返回：Seurat 对象

```r
# 基础用法
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

读取 H5AD 文件，返回 SingleCellExperiment 对象。

```r
read_h5ad_as_sce(path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `path` | character | .h5ad 文件路径 |

返回：SingleCellExperiment 对象

---

### read_rds_fast

使用 Rust 原生解析器快速读取 Seurat RDS 文件。

```r
read_rds_fast(path, normalize = FALSE, top_genes = NULL,
              gene_id_column = NULL, keep_layers = FALSE)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `path` | character | 必需 | .rds 文件路径 |
| `normalize` | logical | FALSE | 库大小归一化 + log1p |
| `top_genes` | integer 或 NULL | NULL | 按方差选择前 N 个高变基因 |
| `gene_id_column` | character 或 NULL | NULL | 用指定 var 列作为基因标识符 |
| `keep_layers` | logical | FALSE | 保留 Seurat V5 split layers |

返回：Seurat 对象

```r
# 基础用法
seurat_obj <- read_rds_fast("seurat.rds")

# AI-Ready
seurat_obj <- read_rds_fast(
  "seurat.rds",
  normalize = TRUE,
  top_genes = 2000L,
  gene_id_column = "ensembl_id"
)

# 保留 split layers
seurat_obj <- read_rds_fast("integrated.rds", keep_layers = TRUE)
```

---

### write_as_h5ad

将 Seurat 或 SingleCellExperiment 对象写入 H5AD 格式。

```r
write_as_h5ad(obj, path, normalize = FALSE, top_genes = NULL, gene_id_column = NULL)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `obj` | Seurat 或 SCE | 必需 | 要写入的对象 |
| `path` | character | 必需 | 输出 .h5ad 文件路径 |
| `normalize` | logical | FALSE | 库大小归一化 + log1p |
| `top_genes` | integer 或 NULL | NULL | 按方差选择前 N 个高变基因 |
| `gene_id_column` | character 或 NULL | NULL | 用指定 var 列作为基因标识符 |

```r
# 基础用法
write_as_h5ad(seurat_obj, "output.h5ad")

# 带变换
write_as_h5ad(seurat_obj, "ai_ready.h5ad",
              normalize = TRUE, top_genes = 2000L)
```

---

### write_rds_fast

快速将 Seurat 对象写入 RDS 格式。

```r
write_rds_fast(obj, path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `obj` | Seurat | 要写入的 Seurat 对象 |
| `path` | character | 输出 .rds 文件路径 |

---

### inspect_file

返回文件元数据列表。

```r
inspect_file(path)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `path` | character | .h5ad 或 .rds 文件路径 |

返回：list，包含 `format`、`n_cells`、`n_genes`、`obs_columns`、`var_columns`、`has_embeddings`、`embedding_names`、`has_layers`、`layer_names`。

```r
info <- inspect_file("data.h5ad")
cat(sprintf("%d cells × %d genes\n", info$n_cells, info$n_genes))
```

---

### validate_conversion

比较两个文件，返回保真度报告。

```r
validate_conversion(original, converted, tolerance = 1e-7, cluster_column = NULL)
```

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `original` | character | 必需 | 原始文件路径 |
| `converted` | character | 必需 | 转换后文件路径 |
| `tolerance` | numeric | 1e-7 | 数值容差 |
| `cluster_column` | character 或 NULL | NULL | obs 中的聚类标签列（启用 ARI/NMI） |

返回：list，包含 `passed`（logical）、`summary`（character）、`detailed_report`（character），以及可选的 `ari`（numeric）、`nmi`（numeric）。

```r
result <- validate_conversion("original.h5ad", "roundtrip.h5ad")
if (result$passed) message("验证通过！")

# 带聚类指标
result <- validate_conversion(
  "original.h5ad", "roundtrip.h5ad",
  cluster_column = "cell_type"
)
cat(sprintf("ARI: %.4f, NMI: %.4f\n", result$ari, result$nmi))
```

## 典型工作流

### H5AD → Seurat 分析

```r
library(crosscell)
library(Seurat)

seurat_obj <- read_h5ad_as_seurat("scanpy_processed.h5ad")
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj, reduction = "umap")
write_as_h5ad(seurat_obj, "seurat_analyzed.h5ad")
```

### AI-Ready 导出

```r
library(crosscell)

# 一步完成：Seurat → 归一化 H5AD + 高变基因
write_as_h5ad(
  seurat_obj, "ai_ready.h5ad",
  normalize = TRUE,
  top_genes = 2000L,
  gene_id_column = "ensembl_id"
)
```

### 往返验证

```r
library(crosscell)

write_as_h5ad(seurat_obj, "converted.h5ad")
result <- validate_conversion("original.rds", "converted.h5ad")
stopifnot(result$passed)
```
