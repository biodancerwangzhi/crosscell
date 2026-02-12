# CrossCell

基于 Rust 的单细胞转录组数据格式转换器。

## 概述

CrossCell 提供 AnnData (.h5ad) 和 Seurat (.rds) 格式之间的高性能双向转换。

### 核心特性

- **纯 Rust 实现**: 无需 Python/R 运行时依赖
- **高性能**: 并行处理、Lazy Loading、内存优化
- **流式处理**: 支持超大数据集（10M+ 细胞），内存占用恒定
- **命令行工具**: 独立可执行文件，易于部署
- **双向转换**: H5AD ↔ RDS 完整支持
- **数据完整性**: 保留表达矩阵、元数据、降维结果

## 安装

### 使用 Docker（推荐）

```bash
# 构建 Release 版本（静态链接 HDF5）
docker-compose run --rm dev cargo build --release

# 二进制文件直接生成在本地 target/release/crosscell
# 在 Linux/WSL 中运行
chmod +x target/release/crosscell
./target/release/crosscell --version
```

### 从源码编译

```bash
cargo build --release
```

## 使用方法

### 格式转换

```bash
# AnnData → Seurat
crosscell convert -i input.h5ad -o output.rds -f seurat

# Seurat → AnnData
crosscell convert -i input.rds -o output.h5ad -f anndata

# 带验证的转换
crosscell convert -i input.h5ad -o output.rds -f seurat --validate

# 内存估算
crosscell convert -i input.h5ad -o output.rds -f seurat --estimate-memory
```

### 大数据集处理

```bash
# 流式转换（推荐用于 10M+ 细胞）
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --chunk-size 10000

# Lazy Loading（适用于 50k-1M 细胞）
crosscell convert -i large.h5ad -o large.rds -f seurat --lazy --chunk-size 10000

# 磁盘缓存模式（超大数据集）
crosscell convert -i huge.h5ad -o huge.rds -f seurat --disk-backed --temp-dir /tmp
```

### 文件检查

```bash
# 检查文件信息
crosscell inspect -i input.h5ad

# 详细信息
crosscell inspect -i input.h5ad --detailed

# 验证转换结果
crosscell validate --original input.h5ad --converted output_roundtrip.h5ad
```

### 查看支持的格式

```bash
crosscell formats
```

## 项目结构

```
crosscell/
├── src/                    # 源代码
│   ├── anndata/            # AnnData (.h5ad) 读写 + 流式处理
│   ├── seurat/             # Seurat 对象处理
│   ├── rds/               # RDS 格式读写
│   ├── ir/                 # 中间表示层
│   ├── sparse/             # 稀疏矩阵操作
│   ├── backend/            # 存储后端抽象
│   ├── storage/            # 磁盘缓存存储
│   ├── validation/         # 验证引擎
│   ├── diagnostics/        # 数据诊断
│   ├── formats/            # 格式插件系统
│   └── bin/                # CLI 入口
├── tests/                  # 测试文件
│   ├── data/               # 测试数据
│   └── test_*.rs           # 集成测试
├── benches/                # 性能基准测试
├── scripts/                # 工具脚本
├── docs/                   # 文档
├── .devcontainer/          # Docker 开发环境
└── .kiro/specs/            # 项目规格
```

## 开发

### 环境准备

```bash
# 验证开发环境
docker-compose run --rm dev bash scripts/verify_environment.sh
```

### 运行测试

```bash
# 运行所有测试
docker-compose run --rm dev cargo test

# 运行特定测试
docker-compose run --rm dev cargo test test_name

# 运行流式处理测试
docker-compose run --rm dev cargo test streaming

# 运行基准测试
docker-compose run --rm dev cargo bench
```

### 工具脚本

```bash
# 生成配对测试数据
docker-compose run --rm dev python3 scripts/generate_paired_data_easyscf.py

# 准确性验证
docker-compose run --rm dev python3 scripts/verify_accuracy.py

# 性能测试
docker-compose run --rm dev python3 scripts/test_performance.py

# 发布前检查
docker-compose run --rm dev python3 scripts/check_no_mock_data.py
```

## 技术栈

- **hdf5**: HDF5 文件读写
- **nalgebra-sparse**: 稀疏矩阵操作
- **rayon**: 并行计算
- **arrow**: 列式数据处理
- **clap**: 命令行界面
- **thiserror**: 错误处理

## 文档

**English**
- [Home](docs/home.md) - Overview
- [Installation](docs/installation.md) - Setup guide
- [Quick Start](docs/quick_start.md) - Get started
- [Advanced Usage](docs/advanced_usage.md) - Full reference
- [FAQ](docs/faq.md) - Common questions
- [Contributing](docs/contributing.md) - Development guide

**中文**
- [首页](docs/zh/home.md) - 概述
- [安装](docs/zh/installation.md) - 安装指南
- [快速入门](docs/zh/quick_start.md) - 快速上手
- [高级用法](docs/zh/advanced_usage.md) - 完整参考
- [常见问题](docs/zh/faq.md) - 常见问题
- [贡献指南](docs/zh/contributing.md) - 开发指南

**Technical**
- [Performance](docs/PERFORMANCE.md) - Performance optimization
- [RDS Spec](docs/RDS_SPEC.md) - RDS format specification

## 许可证

MIT OR Apache-2.0
