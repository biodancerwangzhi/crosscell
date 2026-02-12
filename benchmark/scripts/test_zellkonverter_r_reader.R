#!/usr/bin/env Rscript
# 测试 zellkonverter readH5AD(reader="R") 是否能跑通

# 先确保依赖
for (pkg in c("rhdf5", "HDF5Array")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(zellkonverter)
library(Seurat)

# 测试 H5AD → RDS (reader="R", 纯 R 不需要 basilisk)
test_files <- c(
  "/benchmark/data/generated/scanpy_pbmc3k.h5ad",
  "/benchmark/data/generated/cellxgene_pbmc_15k.h5ad"
)

for (f in test_files) {
  if (!file.exists(f)) {
    cat("SKIP:", f, "(not found)\n")
    next
  }
  cat("\n=== Testing:", basename(f), "===\n")
  tryCatch({
    sce <- readH5AD(f, reader = "R")
    cat("  SCE dim:", dim(sce), "\n")
    
    # 尝试转换为 Seurat
    obj <- as.Seurat(sce, counts = "X", data = NULL)
    cat("  Seurat OK, cells:", ncol(obj), "genes:", nrow(obj), "\n")
    
    # 保存
    out <- paste0("/tmp/test_", basename(f), ".rds")
    saveRDS(obj, out)
    cat("  Saved:", out, "\n")
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n=== All tests done ===\n")
