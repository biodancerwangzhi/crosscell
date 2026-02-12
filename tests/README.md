# Tests 目录

CrossCell 项目的测试文件组织结构。

## 目录结构

```
tests/
├── data/                          # 测试数据
│   ├── real_datasets/            # 真实数据集（PBMC 3k 等）
│   └── *.rds, *.h5ad             # 测试用的小型数据文件
├── test_*.rs                      # Rust 集成测试
├── *_test.py                      # Python 验证脚本
├── *.R                            # R 验证脚本
└── *.sh                           # Shell 测试脚本
```

## 测试文件分类

### 1. Rust 集成测试 (test_*.rs)

核心功能测试，使用 Rust 的测试框架：

- **RDS 格式测试**
  - `test_rds_roundtrip.rs` - RDS 往返测试
  - `test_rds_validation.rs` - RDS 验证测试
  - `test_read_r_files.rs` - 读取 R 文件测试
  - `rds_basic.rs` - RDS 基础功能
  - `rds_r_compat.rs` - R 兼容性测试

- **H5AD 格式测试**
  - `test_h5ad_read.rs` - H5AD 读取测试
  - `test_h5ad_write.rs` - H5AD 写入测试
  - `test_h5ad_roundtrip.rs` - H5AD 往返测试

- **Seurat 对象测试**
  - `test_seurat_to_ir.rs` - Seurat → IR 转换
  - `test_ir_to_seurat.rs` - IR → Seurat 转换
  - `test_seurat_rds_read.rs` - Seurat RDS 读取
  - `test_seurat_spatial.rs` - 空间转录组测试

- **空间转录组测试**
  - `test_spatial_h5ad.rs` - 空间 H5AD 测试
  - `test_spatial_roundtrip.rs` - 空间数据往返
  - `test_spatial_comprehensive.rs` - 综合空间测试

- **数据结构测试**
  - `test_dataframe_read.rs` - DataFrame 读取
  - `test_dgcmatrix_write.rs` - dgCMatrix 写入
  - `test_list_with_*.rs` - 列表结构测试
  - `test_named_list_write.rs` - 命名列表写入

### 2. Python 验证脚本 (*_test.py, verify_*.py)

用于验证 Rust 代码生成的数据：

- `verify_h5ad_write.py` - 验证 H5AD 写入
- `verify_rust_written.py` - 验证 Rust 写入的数据
- `verify_embeddings_layers.py` - 验证降维和层数据
- `verify_test_datasets.py` - 验证测试数据集

### 3. R 验证脚本 (*.R)

用于验证 RDS 文件的 R 兼容性：

- `check_rds.R` - 检查 RDS 文件结构

### 4. Shell 脚本 (*.sh)

自动化测试流程：

- `run_h5ad_tests.sh` - 运行 H5AD 测试套件
- `setup_test_datasets.sh` - 设置测试数据集

## 测试数据

### 小型测试文件 (tests/data/)

用于单元测试的小型数据文件：

- `*.rds` - R 数据文件（各种数据结构）
- `*.h5ad` - AnnData 文件（小型测试数据）
- `empty*.rds` - 空数据测试
- `simple*.rds` - 简单数据结构测试
- `seurat*.rds` - Seurat 对象测试

### 真实数据集 (tests/data/real_datasets/)

用于集成测试的真实数据：

- `pbmc3k.h5ad` - PBMC 3k 数据集（2,700 细胞）
- `pbmc3k_from_h5ad.rds` - 从 h5ad 转换的 RDS
- 其他真实数据集...

详见 `tests/data/real_datasets/README.md`

## 运行测试

### 运行所有 Rust 测试

```bash
docker-compose run --rm dev cargo test
```

### 运行特定测试

```bash
# RDS 测试
docker-compose run --rm dev cargo test rds

# H5AD 测试
docker-compose run --rm dev cargo test h5ad

# 空间转录组测试
docker-compose run --rm dev cargo test spatial
```

### 运行 Python 验证

```bash
docker-compose run --rm dev python3 tests/verify_h5ad_write.py
```

### 运行 R 验证

```bash
docker-compose run --rm dev Rscript tests/check_rds.R
```

## 添加新测试

### 1. Rust 集成测试

在 `tests/` 目录创建 `test_*.rs` 文件：

```rust
#[test]
fn test_my_feature() {
    // 测试代码
}
```

### 2. Python 验证脚本

在 `tests/` 目录创建 `verify_*.py` 文件：

```python
#!/usr/bin/env python3
import scanpy as sc

def verify_my_feature():
    # 验证代码
    pass

if __name__ == "__main__":
    verify_my_feature()
```

### 3. R 验证脚本

在 `tests/` 目录创建 `check_*.R` 文件：

```r
#!/usr/bin/env Rscript

# 验证代码
```

## 测试数据管理

### 下载测试数据

```bash
# 下载真实数据集
python3 scripts/download_test_datasets.py

# 或使用 R
Rscript scripts/download_test_datasets.R
```

### 生成配对数据

```bash
# 使用 easySCF 生成金标准配对数据
python3 scripts/generate_paired_data_easyscf.py
```

## 注意事项

1. **测试隔离**：每个测试应该独立运行，不依赖其他测试
2. **数据清理**：测试后清理临时文件
3. **错误处理**：测试应该有清晰的错误信息
4. **性能**：避免在测试中使用过大的数据集
5. **文档**：为复杂测试添加注释说明

## 相关文档

- [测试指南](../TESTING_GUIDE.md) - 详细的测试指南
- [Docker 测试工作流](.kiro/steering/docker-testing-workflow.md) - Docker 测试规范
- [数据集说明](data/real_datasets/README.md) - 测试数据集详情
