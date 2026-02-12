# CrossCell 测试数据

## 数据来源

### RDS 数据
- **SeuratData**: Satija Lab 提供的单细胞数据集 (http://seurat.nygenome.org)
- 22 个数据包中 11 个是 ref 类型（用于 Azimuth），不是 Seurat 对象，脚本自动跳过

### H5AD 数据
- scanpy, scvelo, squidpy, CELLxGENE 等 Python 生态数据集

## 目录结构

```
data/
├── generated/                      # 输出：生成的 RDS + 下载的 H5AD
├── seurat-packages/                # SeuratData 包 (IDM 手动下载)
├── zenodo/                         # Zenodo 上传包（打包生成）
├── docker-compose.datagen.yml      # Docker 编排
├── Dockerfile.seurat-v4            # Seurat 4.3.0 镜像
├── Dockerfile.seurat-v5            # Seurat 5.0.0 镜像
├── 01_export_download_links.R      # Step 1: 导出下载链接
├── 02_generate_seurat_data.R       # Step 2: 生成 RDS 数据
├── 03_inspect_seurat_data.R        # Step 3: 检查 RDS 并生成 CSV
├── 04_convert_stxKidney_v5.R       # Step 4: stxKidney 空间数据转换
├── 05_download_h5ad_datasets.py    # Step 5: 下载 H5AD 数据集
├── 05_download_h5ad.sh             # Step 5: shell 入口
├── 06_generate_datasets_json.py    # Step 6: 生成 datasets.json 元数据
└── 07_package_for_zenodo.py        # Step 7: 打包上传到 Zenodo
```

## 数据生成流程

### Step 0: 构建 Docker 镜像

```bash
docker-compose -f data/docker-compose.datagen.yml build
```

构建两个镜像：

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

空间坐标保存在 `meta.data` 的 `spatial_x`/`spatial_y` 列中。

### Step 5: 下载 H5AD 数据集

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  bash /benchmark/data/05_download_h5ad.sh
```

从 scanpy, scvelo, squidpy, CELLxGENE 下载 13 个 H5AD 文件到 `generated/`。

### Step 6: 生成 datasets.json

```bash
docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
  python3 /benchmark/data/06_generate_datasets_json.py
```

扫描 `generated/` 中所有 RDS 和 H5AD 文件，生成 `benchmark/results/datasets.json`。

## 上传到 Zenodo

数据生成完成后，打包上传到 Zenodo 供他人下载复现。

### 1. 打包

在项目根目录运行打包脚本：

```bash
python3 data/07_package_for_zenodo.py
```

输出到 `data/zenodo/`，按类别分包方便下载：

| 文件 | 内容 | 大小 |
|------|------|------|
| `seurat_v4_datasets.tar.gz` | 23 个 Seurat V4 RDS | 3.3 GB |
| `seurat_v5_datasets.tar.gz` | 21 个 Seurat V5 RDS | 2.3 GB |
| `h5ad_small_datasets.tar.gz` | scanpy + scvelo H5AD (4 个) | 101 MB |
| `h5ad_spatial_datasets.tar.gz` | squidpy 空间数据 H5AD (6 个) | 190 MB |
| `h5ad_cellxgene_datasets.tar.gz` | CELLxGENE 中大规模 H5AD (3 个) | 613 MB |
| `benchmark_results.tar.gz` | 预计算结果 + 脚本 | < 1 MB |
| `manifest.json` | 文件清单和 SHA256 校验 | - |
| `README.md` | Zenodo 页面描述 | - |

### 2. 创建 Zenodo Deposit

1. 注册 [Zenodo](https://zenodo.org) 账号（可用 GitHub/ORCID 登录）
2. 打开 https://zenodo.org/deposit/new
3. Upload type 选 **Dataset**

### 3. 上传文件

将 `benchmark/zenodo/` 下的所有文件拖入上传区域（6 个 tar.gz + manifest.json + README.md）。

> Zenodo 单文件限制 50 GB，总容量不限。分包后每个文件都在 1.5 GB 以内，下载更稳定。

### 4. 填写元数据

| 字段 | 值 |
|------|-----|
| Title | CrossCell Benchmark: Single-Cell Data Format Conversion |
| Authors | 你的姓名和单位 |
| Description | 从 `data/zenodo/README.md` 复制 |
| License | Creative Commons Attribution 4.0 (数据) |
| Keywords | single-cell, scRNA-seq, benchmark, format conversion, Seurat, AnnData |

### 5. 填写 Related Works

每条填 4 个字段：Relation、Identifier、Scheme、Resource type。

**GitHub 仓库**

| 字段 | 值 |
|------|-----|
| Relation | Is supplemented by |
| Identifier | https://github.com/你的用户名/crosscell |
| Scheme | URL |
| Resource type | Software |

**SeuratData（RDS 数据来源）**

| 字段 | 值 |
|------|-----|
| Relation | Is derived from |
| Identifier | 10.1016/j.cell.2024.04.038 |
| Scheme | DOI |
| Resource type | Dataset |

**CELLxGENE（H5AD 数据来源）**

| 字段 | 值 |
|------|-----|
| Relation | Is derived from |
| Identifier | https://cellxgene.cziscience.com |
| Scheme | URL |
| Resource type | Dataset |

**scanpy（H5AD 数据来源）**

| 字段 | 值 |
|------|-----|
| Relation | Is derived from |
| Identifier | 10.1186/s13059-017-1382-0 |
| Scheme | DOI |
| Resource type | Software |

**scvelo（H5AD 数据来源）**

| 字段 | 值 |
|------|-----|
| Relation | Is derived from |
| Identifier | 10.1038/s41587-020-0591-3 |
| Scheme | DOI |
| Resource type | Software |

**squidpy（H5AD 数据来源）**

| 字段 | 值 |
|------|-----|
| Relation | Is derived from |
| Identifier | 10.1038/s41592-021-01358-2 |
| Scheme | DOI |
| Resource type | Software |

### 6. 发布

1. 点击 **Preview** 检查
2. 点击 **Publish** 发布
3. 获得 DOI，格式如 `10.5281/zenodo.XXXXXXX`
4. 记录 Record ID（URL 中的数字），填入 notebook 的 `ZENODO_RECORD_ID`

> 建议先在 [Zenodo Sandbox](https://sandbox.zenodo.org) 测试一遍。Sandbox 和正式环境完全独立，账号需要分别注册。

## References

1. Hao, Y., Stuart, T., Kowalski, M.H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Mez, G., Papalexi, E., Chang, W.M., et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293–304. https://doi.org/10.1038/s41587-023-01767-y

2. Wolf, F.A., Angerer, P., & Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0

3. Bergen, V., Lange, M., Peidli, S., Wolf, F.A., & Theis, F.J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology*, 38, 1408–1414. https://doi.org/10.1038/s41587-020-0591-3

4. Palla, G., Spitzer, H., Klein, M., Fischer, D., Schaar, A.C., Kuemmerle, L.B., Rybakov, S., Ibarra, I.L., Polański, K., Rber, M., et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178. https://doi.org/10.1038/s41592-021-01358-2

5. CZI Single-Cell Biology Program, Abdill, R.J., Ahlmann-Eltze, C., et al. (2024). CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Research*, 52(D1), D223–D230. https://doi.org/10.1093/nar/gkad884

## 数据集汇总

### RDS 数据集

| 数据集 | 细胞数 | 基因数 | V4 大小 | V5 大小 |
|--------|--------|--------|---------|---------|
| bmcite | 30,672 | 17,009 | 167 MB | 60 MB |
| cbmc | 8,617 | 20,501 | 34 MB | 17 MB |
| celegans.embryo | 6,188 | 20,222 | 21 MB | 11 MB |
| hcabm40k | 40,000 | 17,369 | 180 MB | 90 MB |
| ifnb | 13,999 | 14,053 | 40 MB | 20 MB |
| panc8 | 14,890 | 34,363 | 275 MB | 138 MB |
| pbmc3k | 2,700 | 13,714 | 10 MB | 5 MB |
| pbmcsca | 31,021 | 33,694 | 121 MB | 61 MB |
| ssHippo | 53,173 | 23,264 | 94 MB | 48 MB |
| stxKidney* | 1,438 | 31,053 | 38 MB | 18 MB |
| thp1.eccite | 20,729 | 18,649 | 309 MB | 154 MB |

*stxKidney 是 10x Visium 空间数据，V5 版本通过 `04_convert_stxKidney_v5.R` 单独转换。

### H5AD 数据集

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
