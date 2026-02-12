# CrossCell 脚本工具

本目录包含 CrossCell 项目的工具脚本。

## �️ 运行环境

脚本支持两种运行方式：

### 方式 1：WSL（推荐用于测试）

直接在 WSL 中运行，需要先安装依赖：

```bash
# 安装 Python 依赖
pip install numpy scipy scanpy anndata

# 安装 R 依赖
Rscript -e "install.packages(c('Matrix', 'Seurat'), repos='https://cloud.r-project.org')"
```

### 方式 2：Docker（推荐用于开发）

使用 Docker 容器，环境已预配置：

```bash
docker-compose run --rm dev <命令>
```

---

## 📋 脚本列表

| 脚本 | 用途 | 阶段 |
|------|------|------|
| `verify_environment.sh` | 验证开发环境 | 1. 环境准备 |
| `download_test_datasets.py` | 下载 AnnData 测试数据 | 2. 数据准备 |
| `download_test_datasets.R` | 下载 Seurat 测试数据 | 2. 数据准备 |
| `generate_paired_data_easyscf.py` | 生成配对测试数据 | 2. 数据准备 |
| `generate_seurat_v5_test.R` | 生成 Seurat V5 测试数据 | 2. 数据准备 |
| `verify_accuracy.py` | 转换准确性验证 | 3. 功能测试 |
| `test_performance.py` | 性能优化测试 | 3. 功能测试 |
| `check_no_mock_data.py` | Mock 数据检查 | 4. 发布前检查 |

---

## 1️⃣ 环境验证

### verify_environment.sh

验证开发环境是否正确配置。

```bash
# WSL
bash scripts/verify_environment.sh

# Docker
docker-compose run --rm dev bash scripts/verify_environment.sh
```

**检查内容**：系统工具、Rust 工具链、Python 环境、R 环境、系统库

---

## 2️⃣ 数据准备

### download_test_datasets.py

下载真实测试数据集（AnnData 格式）。

```bash
# WSL
python3 scripts/download_test_datasets.py

# Docker
docker-compose run --rm dev python3 scripts/download_test_datasets.py
```

**输出文件**：
- `tests/data/real_datasets/pbmc3k.h5ad`
- `tests/data/real_datasets/pancreas.h5ad`
- `tests/data/real_datasets/visium_heart.h5ad`

---

### download_test_datasets.R

下载真实测试数据集（Seurat 格式）。

```bash
# WSL
Rscript scripts/download_test_datasets.R

# Docker
docker-compose run --rm dev Rscript scripts/download_test_datasets.R
```

---

### generate_seurat_v5_test.R

生成 Seurat V5 Assay5 测试数据（SimplifiedSeurat 格式）。

```bash
# WSL
Rscript scripts/generate_seurat_v5_test.R

# Docker
docker-compose run --rm dev Rscript scripts/generate_seurat_v5_test.R
```

**输出文件**：
- `tests/data/real_datasets/seurat_v5_small.rds`
- `tests/data/real_datasets/seurat_v5_medium.rds`
- `tests/data/real_datasets/seurat_v5_multi_assay.rds`

---

### generate_paired_data_easyscf.py

使用 easySCF 生成配对的 h5ad 和 rds 数据。

```bash
# WSL（需要安装 easySCFpy）
pip install easySCFpy
python3 scripts/generate_paired_data_easyscf.py

# Docker
docker-compose run --rm dev python3 scripts/generate_paired_data_easyscf.py
```

---

## 3️⃣ 功能测试

### verify_accuracy.py

验证 CrossCell 的格式转换准确性。

```bash
# WSL
python3 scripts/verify_accuracy.py

# Docker
docker-compose run --rm dev python3 scripts/verify_accuracy.py
```

---

### test_performance.py

测试性能优化功能（Lazy Loading、流式处理、并行转换、内存估算）。

```bash
# WSL
python3 scripts/test_performance.py

# Docker
docker-compose run --rm dev python3 scripts/test_performance.py
```

---

### Rust 基准测试

运行 Rust 性能基准测试：

```bash
# 稀疏矩阵转换基准测试
docker-compose run --rm dev cargo bench --bench sparse_conversion

# 流式处理基准测试
docker-compose run --rm dev cargo bench --bench streaming_conversion

# 运行所有基准测试
docker-compose run --rm dev cargo bench
```

**基准测试内容**：
- `sparse_conversion.rs` - CSR/CSC 转换、并行处理、内存估算
- `streaming_conversion.rs` - 流式读取、分块处理、大文件转换

---

## 4️⃣ 发布前检查

### check_no_mock_data.py

检查源码和测试中是否存在 mock 数据。

```bash
# WSL
python3 scripts/check_no_mock_data.py
python3 scripts/check_no_mock_data.py --strict  # 严格模式

# Docker
docker-compose run --rm dev python3 scripts/check_no_mock_data.py
```

---

## 🔄 完整工作流

### WSL 方式

```bash
# 1. 验证环境
bash scripts/verify_environment.sh

# 2. 下载测试数据
python3 scripts/download_test_datasets.py
Rscript scripts/download_test_datasets.R
Rscript scripts/generate_seurat_v5_test.R

# 3. 运行 Rust 测试
cargo test

# 4. 运行准确性验证
python3 scripts/verify_accuracy.py

# 5. 发布前检查
python3 scripts/check_no_mock_data.py --strict
```

### Docker 方式

```bash
# 1. 验证环境
docker-compose run --rm dev bash scripts/verify_environment.sh

# 2. 下载测试数据
docker-compose run --rm dev python3 scripts/download_test_datasets.py
docker-compose run --rm dev Rscript scripts/download_test_datasets.R
docker-compose run --rm dev Rscript scripts/generate_seurat_v5_test.R

# 3. 运行 Rust 测试
docker-compose run --rm dev cargo test

# 4. 运行准确性验证
docker-compose run --rm dev python3 scripts/verify_accuracy.py

# 5. 运行性能测试
docker-compose run --rm dev python3 scripts/test_performance.py

# 6. 运行基准测试（可选）
docker-compose run --rm dev cargo bench

# 7. 发布前检查
docker-compose run --rm dev python3 scripts/check_no_mock_data.py --strict
```

---

## ⚠️ 注意事项

1. **所有脚本都应从项目根目录运行**
2. **WSL 需要先安装依赖**（Python: numpy, scipy, scanpy; R: Matrix, Seurat）
3. **Docker 方式环境已预配置**，无需额外安装

---

**最后更新**: 2026-01-31
