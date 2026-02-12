#!/usr/bin/env Rscript
# 测试 zellkonverter readH5AD 是否能用 reticulate 直接调用

Sys.setenv(RETICULATE_PYTHON = "/opt/venv/bin/python3")
library(reticulate)
use_python("/opt/venv/bin/python3", required = TRUE)

cat("Python config:\n")
print(py_config())

# 直接用 reticulate 读取 H5AD
ad <- import("anndata")
adata <- ad$read_h5ad("/benchmark/data/generated/scanpy_pbmc3k.h5ad")
cat("Shape:", adata$shape[[1]], "x", adata$shape[[2]], "\n")

# 尝试用 zellkonverter 的内部函数转换
library(zellkonverter)
library(SingleCellExperiment)

# 用 AnnData2SCE 直接转换 (绕过 basilisk)
sce <- AnnData2SCE(adata)
cat("SCE dim:", dim(sce), "\n")

# 转换为 Seurat
library(Seurat)
obj <- as.Seurat(sce, counts = "X", data = NULL)
cat("Seurat OK\n")
saveRDS(obj, "/tmp/test_h5ad_to_rds.rds")
cat("Saved OK\n")
