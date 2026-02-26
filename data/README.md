# CrossCell 测试数据

## 概述

本目录包含 CrossCell benchmark 论文使用的 67 个数据集的生成、下载和打包脚本。

| 类别 | 数据集数 | 格式 | 来源 |
|------|---------|------|------|
| Seurat V4 | 21 (11 raw + 10 processed) | RDS | SeuratData |
| Seurat V5 | 21 (11 raw + 10 processed) | RDS | SeuratData |
| scanpy | 2 | H5AD | scanpy built-in |
| scvelo | 2 | H5AD | scvelo built-in |
| squidpy | 6 | H5AD | squidpy/figshare |
| CELLxGENE (小) | 3 | H5AD | CELLxGENE Discover |
| CELLxGENE (benchmark) | 10 | H5AD | CELLxGENE Discover |
| CELLxGENE (百万级) | 2 | H5AD | CELLxGENE Discover |
| **合计** | **67** | | |

## 数据分类依据

### 为什么选择这些数据集？

1. **SeuratData (42 RDS)**：Satija Lab 官方数据包，覆盖 Seurat V4 (Assay) 和 V5 (Assay5) 两种格式，包含 CITE-seq (bmcite/cbmc)、空间数据 (stxKidney)、RNA velocity 等多种数据类型。

2. **Python 生态 H5AD (13)**：覆盖 scanpy/scvelo/squidpy 三大 Python 单细胞工具的内置数据集，包含标准 scRNA-seq、RNA velocity、6 种空间转录组技术。

3. **CELLxGENE benchmark (10)**：分层采样策略，用于大规模 scalability 测试：
   - 组织多样性：blood, brain, gut, heart, kidney, liver, lung, pancreas, skin（9 种组织）
   - 技术多样性：10x 3'/5' v2/v3, 10x multiome, scATAC-seq, Smart-seq2
   - 规模梯度：10k → 37k → 102k → 172k → 428k → 836k cells
   - 疾病状态：normal + disease (COVID-19, BCC tumor)

4. **CELLxGENE 百万级 (2)**：极限 scalability 测试：
   - `cellxgene_eqtl_autoimmune_1.2M`：1,248,980 cells, blood, 10x 3' v2, ~4 GB
   - `cellxgene_immune_1.8M`：1,821,725 cells, blood, ~9.5 GB（64GB RAM OOM，作为 limitation 记录）

### 排除说明

- `cellxgene_retina_244k`：18 GB 文件含异常 dense 矩阵结构（正常 244k 细胞文件应 ~2 GB），inspect 需 62 GB 内存。非工具问题，属于数据异常。

## 目录结构

```
data/
├── generated/                      # 输出：生成的 RDS + 下载的 H5AD
├── seurat-packages/                # SeuratData 包 (IDM 手动下载)
├── zenodo/                         # Zenodo 上传包（打包生成）
├── docker-compose.datagen.yml      # Docker 编排（Seurat V4/V5）
├── Dockerfile.seurat-v4            # Seurat 4.3.0 镜像
├── Dockerfile.seurat-v5            # Seurat 5.0.0 镜像
│
│  # === 数据生成脚本（按执行顺序） ===
├── 01_export_download_links.R      # Step 1: 导出 SeuratData 下载链接
├── 02_generate_seurat_data.R       # Step 2: 生成 RDS 数据
├── 03_inspect_seurat_data.R        # Step 3: 检查 RDS 并生成 CSV
├── 04_convert_stxKidney_v5.R       # Step 4: stxKidney 空间数据转换
├── 05_download_h5ad_datasets.py    # Step 5: 下载 H5AD 数据集
├── 05_download_h5ad.sh             # Step 5: shell 入口
├── 06_generate_datasets_json.py    # Step 6: 生成 datasets.json 元数据
├── 07_package_for_zenodo.py        # Step 7: 打包上传到 Zenodo
├── 08_generate_split_layers_test.R # Step 8: 生成 split layers 测试数据
│
│  # === 辅助脚本 ===
├── 09_query_cellxgene_datasets.py  # CELLxGENE 数据集选择文档和来源信息
├── check_h5ad_files.py             # 检查 H5AD 文件完整性
└── table.csv                       # CELLxGENE 浏览器查询原始导出（参考用）
```

## 完整复现流程

### 前置条件

- Docker + Docker Compose
- 磁盘空间：~30 GB（RDS + 小 H5AD），额外 ~25 GB（CELLxGENE benchmark）

### Step 0: 构建 Docker 镜像

```bash
# Seurat 数据生成镜像
docker-compose -f data/docker-compose.datagen.yml build

# 开发/测试镜像（含 Rust + Python + R）
docker-compose build dev
```

构建两个 Seurat 镜像：

| 服务 | 基础镜像 | Assay 类型 |
|------|----------|------------|
| seurat-v4 | satijalab/seurat:4.3.0 | Assay |
| seurat-v5 | satijalab/seurat:5.0.0 | Assay5 |

### Step 1: 导出 SeuratData 下载链接

```bash
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 \
  Rscript /workspace/export_download_links.R
```

输出 `seurat-packages/seuratdata_download_links.txt`，将链接导入 IDM 批量下载 `*.tar.gz` 到 `seurat-packages/`。

### Step 2: 安装数据包 + 生成 RDS

先在两个容器中分别安装数据包：

```bash
# V4
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v4 bash
# 容器内:
for pkg in /workspace/packages/*.tar.gz; do R CMD INSTALL "$pkg"; done
exit

# V5
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 bash
# 容器内:
for pkg in /workspace/packages/*.tar.gz; do R CMD INSTALL "$pkg"; done
exit
```

然后生成数据：

```bash
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v4 \
  Rscript /workspace/generate_seurat_data.R
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 \
  Rscript /workspace/generate_seurat_data.R
```

每个数据集生成 `*_raw.rds` 和 `*_processed.rds`。

### Step 3: 检查 RDS 数据

```bash
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v4 \
  Rscript /workspace/inspect_seurat_data.R
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 \
  Rscript /workspace/inspect_seurat_data.R
```

生成 `generated/seurat_v4_summary.csv` 和 `generated/seurat_v5_summary.csv`。

### Step 4: 转换 stxKidney 空间数据

stxKidney 使用 SliceImage 类，与 Seurat V5 不兼容，需要单独转换：

```bash
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 \
  Rscript /workspace/convert_stxKidney_v5.R
```

### Step 5: 下载 H5AD 数据集

```bash
# 下载 scanpy + scvelo + squidpy + CELLxGENE 小数据集 (13 files, ~1.5 GB)
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  python3 /benchmark/data/05_download_h5ad_datasets.py --source all

# 下载 CELLxGENE benchmark 数据集 (10 files, ~25 GB)
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  python3 /benchmark/data/05_download_h5ad_datasets.py --source cellxgene-benchmark
```

### Step 6: 生成 datasets.json

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  python3 /benchmark/data/06_generate_datasets_json.py
```

扫描 `generated/` 中所有 RDS 和 H5AD 文件，生成 `benchmark/results/datasets.json`。
大文件（>500 MB）使用 h5py 读取元数据以避免 OOM。

### Step 7: 打包上传到 Zenodo

```bash
python3 data/07_package_for_zenodo.py
```

详见下方 Zenodo 部分。

### Step 8: 生成 Split Layers 测试数据

```bash
docker-compose -f data/docker-compose.datagen.yml run --rm seurat-v5 \
  Rscript /workspace/08_generate_split_layers_test.R
```

## Docker 镜像说明

### 镜像构建和推送

项目使用 3 类 Docker 镜像：

| 镜像 | Dockerfile | 用途 |
|------|-----------|------|
| crosscell-seurat-v4 | `data/Dockerfile.seurat-v4` | 生成 Seurat V4 RDS |
| crosscell-seurat-v5 | `data/Dockerfile.seurat-v5` | 生成 Seurat V5 RDS |
| crosscell-benchmark | `benchmark/Dockerfile.benchmark` | 运行 benchmark |

Seurat 镜像基于 Satija Lab 官方镜像（`satijalab/seurat:4.3.0` / `satijalab/seurat:5.0.0`），
benchmark 镜像包含 Rust 编译的 CrossCell 二进制 + Python + R 环境。

### 在云服务器上使用

```bash
# 方法 1: 本地构建后导出
docker save -o crosscell-benchmark.tar crosscell-benchmark:latest
# 拷贝到服务器
docker load -i crosscell-benchmark.tar

# 方法 2: 从 Docker Hub 拉取（如已推送）
docker pull username/crosscell-benchmark:latest
```

## 上传到 Zenodo

### 打包策略

| 文件 | 内容 | 大小 | 包含在 Zenodo |
|------|------|------|:---:|
| `seurat_v4_datasets.tar.gz` | 21 个 Seurat V4 RDS | ~3.3 GB | ✅ |
| `seurat_v5_datasets.tar.gz` | 21 个 Seurat V5 RDS | ~2.3 GB | ✅ |
| `h5ad_small_datasets.tar.gz` | scanpy + scvelo H5AD (4 个) | ~101 MB | ✅ |
| `h5ad_spatial_datasets.tar.gz` | squidpy 空间数据 H5AD (6 个) | ~190 MB | ✅ |
| `h5ad_cellxgene_datasets.tar.gz` | CELLxGENE 小 H5AD (3 个) | ~613 MB | ✅ |
| `benchmark_results.tar.gz` | 预计算结果 + 脚本 | < 1 MB | ✅ |
| `cellxgene_sources.json` | 10 个大 CELLxGENE 下载 URL | < 1 KB | ✅ |
| `manifest.json` | 文件清单和 SHA256 校验 | < 1 KB | ✅ |
| CELLxGENE benchmark H5AD (10 个) | 10k–836k cells | ~25 GB | ❌ |

10 个 CELLxGENE benchmark 数据集不包含在 Zenodo 中（总计 ~25 GB），因为：
1. 它们可以从 CELLxGENE Discover 免费下载（永久 URL）
2. `cellxgene_sources.json` 提供了所有下载 URL 和 dataset ID
3. `05_download_h5ad_datasets.py --source cellxgene-benchmark` 可一键下载

### Zenodo 元数据

| 字段 | 值 |
|------|-----|
| Title | CrossCell Benchmark: Single-Cell Data Format Conversion |
| License | Creative Commons Attribution 4.0 (数据), MIT (代码) |
| Keywords | single-cell, scRNA-seq, benchmark, format conversion, Seurat, AnnData |

### Colab 下载验证

可以在 Google Colab 中验证 Zenodo 数据下载：

```python
# 在 Colab 中运行
import urllib.request
import json

ZENODO_RECORD_ID = "XXXXXXX"  # 替换为实际 Record ID
BASE_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"

# 获取文件列表
resp = urllib.request.urlopen(BASE_URL)
record = json.loads(resp.read())

print(f"Title: {record['metadata']['title']}")
print(f"DOI: {record['doi']}")
print(f"\nFiles:")
for f in record['files']:
    size_mb = f['size'] / 1024 / 1024
    print(f"  {f['key']:45s} {size_mb:8.1f} MB")
    # 下载: urllib.request.urlretrieve(f['links']['self'], f['key'])
```

下载 CELLxGENE benchmark 数据集：

```python
# 从 cellxgene_sources.json 下载
import json, urllib.request

with open("cellxgene_sources.json") as f:
    sources = json.load(f)

for ds in sources["datasets"]:
    print(f"Downloading {ds['filename']} ({ds['size_mb']} MB)...")
    urllib.request.urlretrieve(ds["url"], ds["filename"])
```

## 数据集汇总

### RDS 数据集 (42 files)

| 数据集 | 细胞数 | 基因数 | V4 大小 | V5 大小 | 特殊特性 |
|--------|--------|--------|---------|---------|----------|
| bmcite | 30,672 | 17,009 | 167 MB | 60 MB | CITE-seq (ADT) |
| cbmc | 8,617 | 20,501 | 34 MB | 17 MB | CITE-seq (ADT) |
| celegans.embryo | 6,188 | 20,222 | 21 MB | 11 MB | |
| hcabm40k | 40,000 | 17,369 | 180 MB | 90 MB | |
| ifnb | 13,999 | 14,053 | 40 MB | 20 MB | |
| panc8 | 14,890 | 34,363 | 275 MB | 138 MB | |
| pbmc3k | 2,700 | 13,714 | 10 MB | 5 MB | |
| pbmcsca | 31,021 | 33,694 | 121 MB | 61 MB | |
| ssHippo | 53,173 | 23,264 | 94 MB | 48 MB | |
| stxKidney* | 1,438 | 31,053 | 38 MB | 18 MB | 10x Visium 空间 |
| thp1.eccite | 20,729 | 18,649 | 309 MB | 154 MB | ECCITE-seq |

*stxKidney V5 版本通过 `04_convert_stxKidney_v5.R` 单独转换。

### H5AD 数据集 — Python 生态 (13 files)

| 数据集 | 来源 | 细胞数 | 基因数 | 大小 |
|--------|------|--------|--------|------|
| scanpy_pbmc3k | scanpy | 2,700 | 32,738 | 6 MB |
| scanpy_pbmc3k_processed | scanpy | 2,638 | 1,838 | 24 MB |
| scvelo_pancreas | scvelo | 3,696 | 27,998 | 50 MB |
| scvelo_dentategyrus | scvelo | 2,930 | 13,913 | 24 MB |
| squidpy_imc | squidpy | 4,668 | 34 | 2 MB |
| squidpy_seqfish | squidpy | 19,416 | 351 | 31 MB |
| squidpy_visium_hne | squidpy | 2,688 | 18,078 | 314 MB |
| squidpy_merfish | squidpy | 73,655 | 161 | 49 MB |
| squidpy_slideseqv2 | squidpy | 41,786 | 4,000 | 251 MB |
| squidpy_mibitof | squidpy | 3,309 | 36 | 19 MB |
| cellxgene_pbmc_15k | CELLxGENE | 14,783 | 17,222 | 71 MB |
| cellxgene_heart_23k | CELLxGENE | 23,483 | 32,330 | 164 MB |
| cellxgene_brain_40k | CELLxGENE | 40,165 | 58,232 | 464 MB |

### H5AD 数据集 — CELLxGENE Benchmark (10 files)

| 数据集 | 细胞数 | 基因数 | 大小 | 组织 | 技术 | 测试环境 |
|--------|--------|--------|------|------|------|----------|
| skin_bcc_10k | 9,841 | 26,886 | 89 MB | skin | multi-assay | 16 GB |
| tabula_liver_22k | 22,214 | 60,606 | 796 MB | liver | Smart-seq2+10x | 16 GB |
| kidney_atacseq_37k | 37,747 | 19,276 | 898 MB | kidney | scATAC-seq | 16 GB |
| brain_multiome_102k | 101,924 | 35,451 | 1.1 GB | brain | 10x multiome | 16 GB |
| pancreas_122k | 121,916 | 32,356 | 915 MB | pancreas | multi-assay | 16 GB |
| brain_dlpfc_172k | 172,120 | 37,490 | 2.2 GB | brain | 10x 3' v3 | 64 GB |
| gut_428k | 428,469 | 32,383 | 5.3 GB | gut | 10x 3'/5' v2 | 64 GB |
| heart_486k | 486,134 | 32,383 | 2.8 GB | heart | 10x 3' v2/v3 | 64 GB |
| hlca_core_585k | 584,944 | 27,402 | 5.5 GB | lung | multi-assay | 64 GB |
| combat_pbmc_836k | 836,148 | 36,306 | 5.2 GB | blood | 10x 5' v1 | 64 GB |

选择标准详见 `09_query_cellxgene_datasets.py`。

### H5AD 数据集 — CELLxGENE 百万级 (2 files)

| 数据集 | 细胞数 | 基因数 | 大小 | 组织 | 技术 | 测试环境 | 结果 |
|--------|--------|--------|------|------|------|----------|------|
| eqtl_autoimmune_1.2M | 1,248,980 | TBD | ~4 GB | blood | 10x 3' v2 | 64 GB | 🔄 测试中 |
| immune_1.8M | 1,821,725 | 32,357 | 9.5 GB | blood | 10x 3' v2 | 64 GB | ❌ OOM (61.9 GB peak) |

参考文献：
- eqtl_autoimmune_1.2M: Yazar et al. "Single-cell eQTL mapping identifies cell type specific genetic control of autoimmune disease"
- immune_1.8M: Domínguez Conde et al. "Cross-tissue immune cell analysis"

## References

1. Hao, Y., et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293–304.
2. Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
3. Bergen, V., et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology*, 38, 1408–1414.
4. Palla, G., et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.
5. CZI Single-Cell Biology Program, et al. (2024). CZ CELLxGENE Discover. *Nucleic Acids Research*, 52(D1), D223–D230.
