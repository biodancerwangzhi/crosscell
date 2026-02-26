#!/usr/bin/env Rscript
# Test anndataR as_AnnData + write_h5ad pipeline

library(anndataR)
library(Seurat)

test_files <- c(
  "seurat_v5_pbmc3k_raw.rds",
  "seurat_v4_pbmc3k_raw.rds",
  "seurat_v5_ifnb_raw.rds",
  "seurat_v4_cbmc_raw.rds"
)

for (fname in test_files) {
  cat(sprintf("\n=== %s ===\n", fname))
  path <- file.path("/benchmark/data/generated", fname)
  out <- file.path("/tmp", sub(".rds", "_anndataR.h5ad", fname))
  
  tryCatch({
    obj <- readRDS(path)
    tryCatch({obj <- UpdateSeuratObject(obj)}, error=function(e){})
    cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(obj), nrow(obj)))
    
    t0 <- proc.time()["elapsed"]
    ad <- as_AnnData(obj)
    cat(sprintf("  as_AnnData: OK (%.1fs)\n", proc.time()["elapsed"] - t0))
    cat(sprintf("  AnnData class: %s\n", class(ad)[1]))
    
    t1 <- proc.time()["elapsed"]
    write_h5ad(ad, out)
    cat(sprintf("  write_h5ad: OK (%.1fs)\n", proc.time()["elapsed"] - t1))
    cat(sprintf("  Output: %.1f MB\n", file.info(out)$size / 1024 / 1024))
    
  }, error=function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}
