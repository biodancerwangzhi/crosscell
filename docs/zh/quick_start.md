# 快速入门

## CLI

```bash
# AnnData → Seurat
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# Seurat → AnnData
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata

# 带验证
crosscell convert -i input.h5ad -o output.rds -f seurat --validate

# AI-Ready（归一化 + 高变基因 + 基因 ID 映射）
crosscell convert -i seurat.rds -o ai_ready.h5ad -f anndata \
  --normalize --top-genes 2000 --gene-id-column ensembl_id

# 检查文件
crosscell inspect -i data.h5ad --detailed

# 验证往返
crosscell validate --original original.h5ad --converted roundtrip.h5ad
```

## Python

```python
import crosscell

# 读取 Seurat RDS → AnnData
adata = crosscell.read_rds("seurat.rds")

# 带 AI-Ready 变换读取
adata = crosscell.read_rds("seurat.rds", normalize=True, top_genes=2000)

# 读取 H5AD
adata = crosscell.read_h5ad("data.h5ad")

# 写入
crosscell.write_h5ad(adata, "output.h5ad")
crosscell.write_rds(adata, "output.rds")

# 检查
info = crosscell.inspect("data.h5ad")
print(f"{info['n_cells']} cells × {info['n_genes']} genes")

# 验证
result = crosscell.validate("original.h5ad", "roundtrip.h5ad")
print(result["passed"])  # True
```

## R

```r
library(crosscell)

# 读取 H5AD → Seurat
seurat_obj <- read_h5ad_as_seurat("data.h5ad")

# 带 AI-Ready 变换读取
seurat_obj <- read_h5ad_as_seurat("data.h5ad", normalize = TRUE, top_genes = 2000L)

# 快速读取 Seurat RDS
seurat_obj <- read_rds_fast("seurat.rds")

# 保留 split layers 读取
seurat_obj <- read_rds_fast("integrated.rds", keep_layers = TRUE)

# 写入
write_as_h5ad(seurat_obj, "output.h5ad")
write_rds_fast(seurat_obj, "output.rds")

# 检查
info <- inspect_file("data.h5ad")

# 验证
result <- validate_conversion("original.h5ad", "roundtrip.h5ad")
```

## 跨生态系统工作流

### R → Python

```r
# R: 预处理
library(Seurat); library(crosscell)
seurat_obj <- readRDS("raw.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
write_as_h5ad(seurat_obj, "processed.h5ad")
```

```python
# Python: 继续分析
import crosscell, scanpy as sc
adata = crosscell.read_h5ad("processed.h5ad")
sc.tl.leiden(adata)
sc.pl.umap(adata, color="leiden")
```

### Python → R

```python
# Python: 预处理
import crosscell, scanpy as sc
adata = crosscell.read_rds("seurat.rds")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
crosscell.write_h5ad(adata, "processed.h5ad")
```

```r
# R: 继续分析
library(crosscell); library(Seurat)
seurat_obj <- read_h5ad_as_seurat("processed.h5ad")
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj)
```

## 下一步

- [CLI 参考](cli_reference.md) — 所有命令和参数
- [Python API](python_api.md) — Python 包完整参考
- [R API](r_api.md) — R 包完整参考
- [高级用法](advanced_usage.md) — 流式处理、批量转换、AI-Ready、CELLxGENE
