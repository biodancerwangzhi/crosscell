# 常见问题

## 通用问题

### 支持哪些格式？

- **AnnData** (.h5ad) - Python/Scanpy 格式
- **Seurat** (.rds) - R/Seurat 格式 (V3, V4, V5)

### 需要安装 Python/R 吗？

不需要。CrossCell 是独立的 Rust 二进制文件，无运行时依赖。

### 哪些数据会被保留？

表达矩阵、细胞/基因元数据、降维结果（PCA、UMAP、tSNE）、layers、类别变量和空间坐标都会完整保留。

---

## 转换问题

### "内存不足" 错误

对大数据集使用延迟加载：

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --lazy --chunk-size 5000
```

或限制线程数：

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat -t 4
```

### "不支持的格式" 错误

确保文件扩展名与内容匹配：
- `.h5ad` 用于 AnnData
- `.rds` 用于 Seurat

### Seurat 转换失败

复杂的 Seurat 对象可能包含运行时对象（闭包、环境）。先简化：

```bash
crosscell simplify-seurat -i complex.rds -o simplified.rds
crosscell convert -i simplified.rds -o output.h5ad -f anndata
```

### 验证失败

尝试放宽容差：

```bash
crosscell validate --original a.h5ad --converted b.h5ad --tolerance 1e-5
```

---

## 性能

### 如何加速转换？

1. 使用 LZ4 压缩：`-c lz4`
2. 禁用压缩：`-c none`
3. 使用所有 CPU 核心（默认）

### 如何减少内存使用？

1. 启用延迟加载：`--lazy`
2. 减小块大小：`--chunk-size 2000`
3. 限制线程数：`-t 4`

### 预计转换时间？

| 数据规模 | 时间 |
|----------|------|
| 3k 细胞 | ~5 秒 |
| 20k 细胞 | ~30 秒 |
| 100k 细胞 | ~3 分钟 |

---

## 数据完整性

### 会丢失数据吗？

部分 `.uns`（非结构化）数据可能不会保留。核心数据（表达矩阵、元数据、降维结果）会完整保留。

### 如何验证转换？

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

或手动验证：

```bash
crosscell validate --original input.h5ad --converted roundtrip.h5ad
```

### 可接受的容差是多少？

默认容差是 `1e-7`。对于大多数用例，`1e-5` 是可接受的。

---

## Docker

### 为什么使用 Docker？

Docker 确保跨平台的一致环境，避免依赖问题。

### 如何运行命令？

```bash
docker-compose run --rm dev crosscell <命令>
```

### 如何访问文件？

工作区挂载在 `/workspace`。使用相对路径：

```bash
docker-compose run --rm dev crosscell convert -i data/input.h5ad -o data/output.rds -f seurat
```
