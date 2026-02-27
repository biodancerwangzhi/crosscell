# CLI 参考

## 全局选项

| 选项 | 简写 | 说明 |
|------|------|------|
| `--verbose` | `-v` | 启用 debug 级别输出 |
| `--quiet` | `-q` | 仅显示错误（与 `-v` 冲突） |
| `--threads` | `-t` | 线程数（默认：所有核心） |
| `--log-file` | | 日志写入文件 |
| `--version` | `-V` | 显示版本 |
| `--help` | `-h` | 显示帮助 |

## 命令

### convert

H5AD 和 RDS 格式之间的转换。

```bash
crosscell convert -i <INPUT> -o <OUTPUT> -f <FORMAT> [OPTIONS]
```

#### 必需参数

| 选项 | 简写 | 说明 |
|------|------|------|
| `--input` | `-i` | 输入文件路径（.h5ad 或 .rds） |
| `--output` | `-o` | 输出文件路径 |
| `--format` | `-f` | 目标格式：`anndata`、`seurat` |

#### AI-Ready 变换

执行顺序固定为：normalize → top-genes → gene-id-column。

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--normalize` | false | 库大小归一化（缩放到 10k）+ log1p |
| `--top-genes` | — | 按方差选择前 N 个基因 |
| `--gene-id-column` | — | 用指定 var 列替换基因索引 |

#### Seurat V5

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--keep-layers` | false | 保留 split layers 为独立的 AnnData layers |
| `--direct` | true | 直接 Rust 原生 RDS 解析（默认） |
| `--simplify-first` | false | 旧模式：需要 R 预处理的简化 RDS |

#### 质量控制

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--validate` | false | 转换后运行往返验证 |
| `--dry-run` | false | 预览转换，不写入输出 |
| `--estimate-memory` | false | 打印内存使用估算 |

#### 性能

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--streaming` | false | 恒定内存流式模式（1M+ 细胞） |
| `--lazy` | false | 延迟加载模式（50k–1M 细胞） |
| `--chunk-size` | 10000 | 每块行数 |
| `--disk-backed` | false | 磁盘缓存模式（10M+ 细胞） |
| `--temp-dir` | 系统默认 | 磁盘缓存临时目录 |

#### 压缩

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--compression` / `-c` | gzip | 算法：`none`、`gzip`、`lz4` |
| `--compression-level` | 6 | 级别 0–9（仅 gzip） |

#### 稀疏矩阵

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--sparse-format` | auto | 输出稀疏格式：`csr`、`csc`、`auto` |

#### 示例

```bash
# 基础转换
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# AI-Ready 管道
crosscell convert -i seurat.rds -o ai_ready.h5ad -f anndata \
  --normalize --top-genes 2000 --gene-id-column ensembl_id

# 保留 Seurat V5 split layers
crosscell convert -i integrated.rds -o integrated.h5ad -f anndata --keep-layers

# 大数据集流式处理
crosscell convert -i big.h5ad -o big.rds -f seurat --streaming --chunk-size 20000

# 带验证的转换
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

---

### inspect

显示文件元数据。

```bash
crosscell inspect -i <INPUT> [OPTIONS]
```

| 选项 | 简写 | 说明 |
|------|------|------|
| `--input` | `-i` | 输入文件路径 |
| `--detailed` | `-d` | 显示详细信息 |
| `--output` | `-o` | 保存报告为 JSON |

---

### validate

比较原始文件和转换后文件的保真度。

```bash
crosscell validate --original <FILE> --converted <FILE> [OPTIONS]
```

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--original` | 必需 | 原始文件路径 |
| `--converted` | 必需 | 转换后文件路径 |
| `--tolerance` | 1e-7 | 浮点数值容差 |
| `--output` / `-o` | — | 保存报告为 JSON |
| `--strict` | false | 任何差异都返回非零退出码 |
| `--cluster-column` | — | obs 中的聚类标签列（启用 ARI/NMI） |
| `--validate-schema` | — | 验证的 Schema（如 `cellxgene-v5`） |

报告指标：Pearson R、Spearman ρ、MSE、RMSE、Max Diff、Mean Diff、Exact Match %，以及可选的 ARI + NMI。

---

### simplify-seurat

移除 Seurat RDS 文件中的运行时对象（闭包、环境、字节码）。

```bash
crosscell simplify-seurat -i <INPUT> -o <OUTPUT> [OPTIONS]
```

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--input` / `-i` | 必需 | 输入 Seurat RDS 文件 |
| `--output` / `-o` | 必需 | 输出简化 RDS 文件 |
| `--keep-graphs` | false | 保留图对象 |
| `--keep-commands` | false | 保留命令历史 |
| `--keep-tools` | false | 保留工具函数 |

---

### formats

列出所有支持的格式。

```bash
crosscell formats
```

## 环境变量

| 变量 | 说明 |
|------|------|
| `CROSSCELL_THREADS` | 默认线程数 |
| `CROSSCELL_LOG_LEVEL` | 日志级别：error、warn、info、debug、trace |
| `RAYON_NUM_THREADS` | Rayon 线程池大小 |

## 退出码

| 代码 | 含义 |
|------|------|
| 0 | 成功 |
| 1 | 一般错误 |
| 2 | 参数错误 |
| 3 | 文件不存在 |
| 4 | 不支持的格式 |
| 5 | 内存不足 |
| 6 | 验证失败 |
| 7 | 转换失败 |
| 8 | 用户中断 (Ctrl+C) |
