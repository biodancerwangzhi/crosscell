# CrossCell

高性能、零依赖的单细胞数据格式转换器。

## 概述

CrossCell 实现 **AnnData (.h5ad)** 和 **Seurat (.rds)** 格式之间的双向转换，提供三种接口：

- **CLI** — 独立二进制文件，无运行时依赖
- **Python 包** — 通过 PyO3 原生集成 AnnData
- **R 包** — 通过 extendr 原生集成 Seurat/SCE

三种接口共享同一个 Rust 核心引擎。

## 核心特性

| 特性 | 说明 |
|------|------|
| 零依赖 | 无需 R 或 Python 运行时 |
| 高保真度 | 往返转换 Pearson R = 1.0，MSE < 10⁻¹⁴ |
| Seurat V5 | 自动合并 Split Layers，可选 `--keep-layers` |
| AI-Ready | 内置归一化、高变基因筛选、基因 ID 映射 |
| CELLxGENE | Schema 5.x 验证 + 精确类型映射 |
| 流式处理 | 恒定内存模式，支持 1M+ 细胞 |
| 验证 | 往返保真度审计（9 项指标 + ARI/NMI） |

## 快速示例

```bash
# CLI
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata --normalize --top-genes 2000
```

```python
# Python
import crosscell
adata = crosscell.read_rds("seurat.rds", normalize=True, top_genes=2000)
```

```r
# R
library(crosscell)
seurat_obj <- read_h5ad_as_seurat("data.h5ad", normalize = TRUE, top_genes = 2000L)
```

## 支持的数据组件

| 组件 | AnnData | Seurat | 状态 |
|------|---------|--------|------|
| 表达矩阵 | X | counts/data | ✅ |
| 细胞元数据 | obs | meta.data | ✅ |
| 基因元数据 | var | meta.features | ✅ |
| PCA/UMAP/tSNE | obsm | reductions | ✅ |
| Layers | layers | layers | ✅ |
| 类别变量 | categorical | factors | ✅ |
| 空间坐标 | obsm['spatial'] | images | ✅ |
| Seurat V5 Assay5 | — | Assay5 | ✅ |
| Split Layers | — | counts.1/2/… | ✅ |

## 文档

| 页面 | 说明 |
|------|------|
| [安装](installation.md) | CLI、Python、R 安装 |
| [快速入门](quick_start.md) | 2 分钟完成首次转换 |
| [CLI 参考](cli_reference.md) | 所有命令和参数 |
| [Python API](python_api.md) | Python 包参考 |
| [R API](r_api.md) | R 包参考 |
| [高级用法](advanced_usage.md) | 流式处理、批量转换、AI-Ready、CELLxGENE |
| [常见问题](faq.md) | 常见问题解答 |
| [贡献指南](contributing.md) | 开发指南 |

## 安装

```bash
# CLI
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-linux-x64.tar.gz | tar xz

# Python
pip install crosscell

# R
install.packages("crosscell", repos = "https://biodancerwangzhi.r-universe.dev")
```

## 许可证

MIT OR Apache-2.0
