# CrossCell Benchmark 数据生成

## 概述

使用两个独立的 Docker 镜像生成 Seurat V4/V5 格式的测试数据。

基于 **satijalab/seurat** 官方镜像，已预装 Seurat，构建速度快。

| 服务 | 基础镜像 | Seurat 版本 |
|------|----------|-------------|
| datagen-v4 | satijalab/seurat:4.3.0 | 4.3.0 |
| datagen-v5 | satijalab/seurat:5.0.0 | 5.0.0 |

> 注：官方没有 Seurat V3 镜像，暂不支持 V3 数据生成。

## ⚠️ 网络问题

SeuratData 数据托管在美国服务器，国内下载不稳定。需要用 IDM 手动下载。

---

## 完整流程

### 步骤 1：构建镜像

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml build datagen-v4 datagen-v5
```

### 步骤 2：导出下载链接

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm datagen-v5 \
  Rscript /workspace/scripts/export_download_links.R
```

链接文件生成到 `benchmark/data/` 目录。

### 步骤 3：用 IDM 下载数据包

1. 打开 `benchmark/data/seuratdata_download_links.txt`
2. 导入 IDM 下载所有 `.tar.gz` 文件
3. 保存到 `benchmark/data/` 目录

推荐下载：pbmc3k, panc8, cbmc, ifnb, bmcite

### 步骤 4：安装数据包（每个版本都要）

```bash
# V4
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm datagen-v4 bash
# 在容器内执行:
for pkg in /workspace/data/*.tar.gz; do R CMD INSTALL "$pkg"; done
exit

# V5
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm datagen-v5 bash
# 在容器内执行:
for pkg in /workspace/data/*.tar.gz; do R CMD INSTALL "$pkg"; done
exit
```

数据包安装到 Docker volume，下次运行时自动保留。

### 步骤 5：生成测试数据

```bash
# 生成 V4 格式
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm datagen-v4 \
  Rscript /workspace/scripts/generate_data.R

# 生成 V5 格式
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm datagen-v5 \
  Rscript /workspace/scripts/generate_data.R
```

生成的文件保存在 `data/generated/` 目录：
- `seurat_v4_pbmc3k_raw.rds`
- `seurat_v4_pbmc3k_processed.rds`
- `seurat_v5_pbmc3k_raw.rds`
- ...

---

## 脚本列表

| 脚本 | 用途 |
|------|------|
| `export_download_links.R` | 导出 SeuratData 下载链接 |
| `generate_data.R` | 生成测试数据（自动检测版本） |
| `run_benchmark.py` | 运行性能测试 |
| `verify_tools.sh` | 验证工具安装 |

---

## 清理 Docker Volume

如需重新安装数据包：

```bash
docker volume rm crosscell-datagen-v4-libs crosscell-datagen-v5-libs
```
