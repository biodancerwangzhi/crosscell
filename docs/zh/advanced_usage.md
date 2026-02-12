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

```python
import subprocess
import scanpy as sc

# 在 Python 中处理
adata = sc.read_h5ad('raw.h5ad')
sc.pp.normalize_total(adata)
sc.tl.pca(adata)
adata.write('processed.h5ad')

# 转换为 Seurat
subprocess.run(['crosscell', 'convert', '-i', 'processed.h5ad', '-o', 'for_r.rds', '-f', 'seurat'])
```

## R 集成

```r
# 从 Python 转换后
seurat_obj <- readRDS('from_python.rds')
seurat_obj <- FindClusters(seurat_obj)
saveRDS(seurat_obj, 'analyzed.rds')

# 转换回 Python（在 shell 中运行）
# crosscell simplify-seurat -i analyzed.rds -o simplified.rds
# crosscell convert -i simplified.rds -o back_to_python.h5ad -f anndata
```

## 空间转录组学

```bash
# 转换 Visium 数据
crosscell convert -i visium.h5ad -o visium.rds -f seurat

# 验证空间坐标
crosscell validate --original visium.h5ad --converted visium_roundtrip.h5ad
```

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

## validate 命令选项

| 选项 | 默认值 | 描述 |
|------|--------|------|
| `--original` | 必需 | 原始文件 |
| `--converted` | 必需 | 转换后文件 |
| `--tolerance` | 1e-7 | 数值容差 |
| `-o, --output` | | 保存报告为 JSON |
| `--strict` | false | 任何差异都失败 |

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
