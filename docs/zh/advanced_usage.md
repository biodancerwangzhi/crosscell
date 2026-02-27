# 高级用法

## 大数据集处理

CrossCell 提供多种模式处理大数据集：

| 模式 | 细胞数 | 内存 | 速度 | 适用场景 |
|------|--------|------|------|----------|
| 标准 | < 50k | 高 | 最快 | 小数据集 |
| 延迟加载 | 50k - 1M | 中 | 快 | 中等数据集 |
| 流式处理 | 1M+ | 恒定 | 良好 | 大数据集 |
| 磁盘缓存 | 10M+ | 最小 | 较慢 | 超大数据集 |

### 内存估算

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --estimate-memory
```

输出：
```
📊 Memory Estimation:
   File size: 5.2 GB
   Full load: ~15.6 GB
   Lazy load: ~2.6 GB (peak)
   Streaming: ~2.4 GB (constant)

💡 Tip: Use --streaming flag for large datasets
```

### 流式模式

对于 1M+ 细胞的数据集，流式模式以恒定内存分块处理数据：

```bash
# 基本流式处理
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming

# 自定义块大小（默认：10000）
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --chunk-size 20000

# 预览模式
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --dry-run
```

流式模式支持：
- H5AD → H5AD（重新分块）
- H5AD → RDS
- RDS → H5AD

### 延迟加载

适用于 50k-1M 细胞的数据集：

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --lazy --chunk-size 5000
```

### 磁盘缓存模式

适用于超大数据集（10M+ 细胞）：

```bash
# 使用磁盘缓存
crosscell convert -i huge.h5ad -o huge.rds -f seurat --disk-backed

# 自定义临时目录
crosscell convert -i huge.h5ad -o huge.rds -f seurat --disk-backed --temp-dir /fast/ssd
```

### 压缩选项

```bash
# 快速 (LZ4)
crosscell convert -i input.h5ad -o output.rds -f seurat -c lz4

# 最大压缩
crosscell convert -i input.h5ad -o output.rds -f seurat -c gzip --compression-level 9

# 无压缩（最快）
crosscell convert -i input.h5ad -o output.rds -f seurat -c none
```

## 批量处理

```bash
#!/bin/bash
for file in data/*.h5ad; do
    name=$(basename "$file" .h5ad)
    crosscell convert -i "$file" -o "output/${name}.rds" -f seurat --validate
done
```

## Python 集成

Python API 现在支持与 CLI 相同的 AI-Ready 变换：

```python
import crosscell

# AI-Ready 管道（无需 subprocess）
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# 验证往返
result = crosscell.validate("original.h5ad", "roundtrip.h5ad", cluster_column="cell_type")
print(f"Passed: {result['passed']}, ARI: {result.get('ari', 'N/A')}")
```

详见 [Python API 参考](python_api.md)。

## R 集成

R API 同样支持 AI-Ready 变换：

```r
library(crosscell)

# AI-Ready 导出
write_as_h5ad(seurat_obj, "ai_ready.h5ad",
              normalize = TRUE, top_genes = 2000L, gene_id_column = "ensembl_id")

# 验证往返
result <- validate_conversion("original.rds", "converted.h5ad", cluster_column = "cell_type")
cat(sprintf("Passed: %s, ARI: %.4f\n", result$passed, result$ari))
```

详见 [R API 参考](r_api.md)。

## 空间转录组学

```bash
# 转换 Visium 数据
crosscell convert -i visium.h5ad -o visium.rds -f seurat

# 验证空间坐标
crosscell validate --original visium.h5ad --converted visium_roundtrip.h5ad
```

## AI-Ready 导出

CrossCell 可以在一次转换中为 AI 基础模型（Geneformer、scGPT、scBERT 等）准备数据。三个变换参数可自由组合，执行顺序固定为：normalize → top-genes → gene-id-column。

### 库大小归一化

对每个细胞执行库大小归一化（缩放到 10,000）+ log1p 变换：

```bash
crosscell convert -i raw.h5ad -o normalized.h5ad -f anndata --normalize
```

全零细胞（空液滴）保持为零，避免除零错误。

### 高变基因筛选

按方差降序保留前 N 个基因：

```bash
crosscell convert -i raw.h5ad -o hvg.h5ad -f anndata --top-genes 2000
```

如果 N 超过总基因数，保留所有基因并输出警告。

### 基因标识符替换

用 `var` 中指定列的值替换基因索引：

```bash
crosscell convert -i raw.h5ad -o ensembl.h5ad -f anndata --gene-id-column ensembl_id
```

### 组合示例

一条命令完成典型的 AI-Ready 管道：

```bash
crosscell convert \
  -i seurat_obj.rds \
  -o ai_ready.h5ad \
  -f anndata \
  --normalize \
  --top-genes 2000 \
  --gene-id-column ensembl_id
```

## CELLxGENE Schema 验证

在提交数据到 CELLxGENE Discover 之前，验证数据集是否符合 Schema 5.0 要求：

```bash
crosscell validate \
  --original data.h5ad \
  --converted data.h5ad \
  --validate-schema cellxgene-v5
```

检查 obs 的 12 个必需字段（如 `assay_ontology_term_id`、`cell_type_ontology_term_id`、`donor_id`、`organism_ontology_term_id` 等）和 var 的 1 个必需字段（`feature_is_filtered`）。

## 聚类感知验证（ARI / NMI）

验证往返转换时，可以通过指定聚类标签列来衡量对聚类结果的影响：

```bash
crosscell validate \
  --original original.h5ad \
  --converted roundtrip.h5ad \
  --cluster-column cell_type
```

这会在标准指标（Pearson 相关系数、MSE、精确匹配率等）之外，额外计算调整兰德指数（ARI）和归一化互信息（NMI）。

## convert 命令选项

| 选项 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | 必需 | 输入文件 |
| `-o, --output` | 必需 | 输出文件 |
| `-f, --format` | 必需 | 目标格式 (anndata/seurat/loom) |
| `--validate` | false | 转换后验证 |
| `--sparse-format` | auto | 稀疏格式 (csr/csc/auto) |
| `-c, --compression` | gzip | 压缩 (none/gzip/lz4) |
| `--compression-level` | 6 | 压缩级别 (0-9) |
| `--lazy` | false | 启用延迟加载 |
| `--streaming` | false | 启用流式模式 |
| `--chunk-size` | 10000 | 延迟加载/流式块大小 |
| `--disk-backed` | false | 启用磁盘缓存模式 |
| `--temp-dir` | 系统默认 | 磁盘缓存临时目录 |
| `--estimate-memory` | false | 估算内存使用 |
| `--dry-run` | false | 预览不转换 |
| `--direct` | true | 直接读取 Seurat（默认） |
| `--simplify-first` | false | 旧模式，需要 R 预处理 |
| `--normalize` | false | 库大小归一化 + log1p |
| `--top-genes` | | 按方差选择前 N 个高变基因 |
| `--gene-id-column` | | 用指定 var 列作为基因标识符 |
| `--keep-layers` | false | 保留 Seurat V5 split layers 为独立的 AnnData layers |

## validate 命令选项

| 选项 | 默认值 | 描述 |
|------|--------|------|
| `--original` | 必需 | 原始文件 |
| `--converted` | 必需 | 转换后文件 |
| `--tolerance` | 1e-7 | 数值容差 |
| `-o, --output` | | 保存报告为 JSON |
| `--strict` | false | 任何差异都失败 |
| `--cluster-column` | | obs 中的聚类标签列（启用 ARI/NMI 计算） |
| `--validate-schema` | | 验证的 Schema 名称（如 `cellxgene-v5`） |

## 环境变量

| 变量 | 描述 |
|------|------|
| `CROSSCELL_THREADS` | 默认线程数 |
| `CROSSCELL_LOG_LEVEL` | 日志级别 (error/warn/info/debug) |
| `RAYON_NUM_THREADS` | Rayon 线程池大小 |

## 退出码

| 代码 | 描述 |
|------|------|
| 0 | 成功 |
| 1 | 一般错误 |
| 2 | 参数错误 |
| 3 | 文件不存在 |
| 4 | 不支持的格式 |
| 5 | 内存不足 |
| 6 | 验证失败 |
| 7 | 转换失败 |

## 性能优化建议

1. **大文件使用流式模式**：> 5GB 的文件会自动建议流式模式
2. **调整块大小**：更大的块 = 更快但更多内存
3. **磁盘缓存使用 SSD**：比 HDD 快很多
4. **选择性验证**：对于可信转换跳过 `--validate`
5. **使用 LZ4 压缩**：大多数情况下速度/大小最佳平衡
