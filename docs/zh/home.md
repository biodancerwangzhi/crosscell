# CrossCell

高性能单细胞转录组数据格式转换器。

## 概述

CrossCell 实现 **AnnData (.h5ad)** 和 **Seurat (.rds)** 格式之间的双向转换：

- **纯 Rust 实现** - 无需 Python/R 运行时依赖
- **高性能** - 并行处理、延迟加载、内存优化
- **流式处理** - 支持超大数据集（10M+ 细胞），内存占用恒定
- **数据完整性** - 保留表达矩阵、元数据、降维结果、空间数据
- **命令行工具** - 独立可执行文件，易于部署

## 支持的数据

| 组件 | AnnData | Seurat | 状态 |
|------|---------|--------|------|
| 表达矩阵 | X | counts/data | ✅ |
| 细胞元数据 | obs | meta.data | ✅ |
| 基因元数据 | var | meta.features | ✅ |
| PCA/UMAP/tSNE | obsm | reductions | ✅ |
| Layers | layers | layers | ✅ |
| 类别变量 | categorical | factors | ✅ |
| 空间坐标 | obsm['spatial'] | images | ✅ |
| Seurat V5 Assay5 | - | Assay5 | ✅ |

## 快速示例

```bash
# AnnData → Seurat
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# Seurat → AnnData  
crosscell convert -i seurat.rds -o anndata.h5ad -f anndata

# 流式模式处理大文件
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming

# 检查文件
crosscell inspect -i pbmc3k.h5ad
```

## 处理模式

| 模式 | 适用场景 | 内存占用 |
|------|----------|----------|
| 标准模式 | < 50k 细胞 | 完整数据集 |
| 延迟加载 | 50k - 1M 细胞 | 较低 |
| 流式处理 | 1M+ 细胞 | 恒定（~2-5 GB） |
| 磁盘缓存 | 10M+ 细胞 | 最小 |

## 文档

- [安装](installation.md)
- [快速入门](quick_start.md)
- [高级用法](advanced_usage.md)
- [常见问题](faq.md)
- [贡献指南](contributing.md)

## 许可证

MIT OR Apache-2.0
