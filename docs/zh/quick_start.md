# 快速入门

## 基本转换

### AnnData 转 Seurat

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat
```

### Seurat 转 AnnData

```bash
crosscell convert -i input.rds -o output.h5ad -f anndata
```

### 带验证的转换

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

## 检查文件

```bash
# 基本信息
crosscell inspect -i data.h5ad

# 详细信息
crosscell inspect -i data.h5ad --detailed
```

输出：
```
📊 文件信息
格式: AnnData (h5ad)
细胞数: 2,700
基因数: 32,738
稀疏度: 93.4%

✅ 数据组件
  ├─ 表达矩阵 (X): 稀疏 CSR
  ├─ 细胞元数据 (obs): 5 列
  ├─ 降维结果 (obsm): PCA, UMAP
  └─ Layers: counts, log1p
```

## 大数据集处理

### 流式模式（推荐用于 1M+ 细胞）

```bash
# 启用流式处理，自定义块大小
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --chunk-size 10000
```

输出：
```
🔄 CrossCell Streaming Conversion
   huge.h5ad → huge.rds
   Format: ANNDATA → SEURAT
   Mode: Streaming (constant memory)
   Chunk size: 10000 cells

📊 File Information:
   File size: 15.2 GB
   Estimated cells: ~5,000,000
   Estimated memory: ~2.4 GB

   Progress: 100% (5,000,000/5,000,000 cells)

✅ Streaming conversion completed!
   Time: 245.32s
   Output size: 12.8 GB
```

### 自动检测

对于 > 5GB 的文件，CrossCell 会自动建议使用流式模式：

```bash
crosscell convert -i huge.h5ad -o huge.rds -f seurat
# ⚠ Large file detected (15 GB, ~5M cells)
# 💡 Using streaming mode with chunk_size=10000
```

## 验证往返转换

```bash
# 双向转换
crosscell convert -i original.h5ad -o converted.rds -f seurat
crosscell convert -i converted.rds -o roundtrip.h5ad -f anndata

# 验证
crosscell validate --original original.h5ad --converted roundtrip.h5ad
```

## 简化 Seurat 对象

复杂的 Seurat 对象可能需要先简化再转换：

```bash
# 简化
crosscell simplify-seurat -i complex.rds -o simplified.rds

# 然后转换
crosscell convert -i simplified.rds -o output.h5ad -f anndata
```

## 查看支持的格式

```bash
crosscell formats
```

输出：
```
📋 Supported Formats
─────────────────────────────────────────

  anndata      (.h5ad)  [read/write]
               AnnData format for Python/Scanpy
  seurat       (.rds)   [read/write]
               Seurat format for R
  loom         (.loom)  [read/write]
               Loom format
```

## 命令参考

| 命令 | 描述 |
|------|------|
| `convert` | 格式转换 |
| `inspect` | 显示文件信息 |
| `validate` | 验证转换准确性 |
| `simplify-seurat` | 简化 Seurat 对象 |
| `formats` | 列出支持的格式 |

## 全局选项

| 选项 | 描述 |
|------|------|
| `-v, --verbose` | 详细输出 |
| `-q, --quiet` | 仅显示错误 |
| `-t, --threads` | 线程数 |
| `-h, --help` | 显示帮助 |
| `-V, --version` | 显示版本 |
