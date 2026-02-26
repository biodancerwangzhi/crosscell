#!/usr/bin/env Rscript
# Quick test: anndataR RDS→H5AD API check

library(anndataR)
library(Seurat)

cat("anndataR version:", as.character(packageVersion("anndataR")), "\n")

# List exported functions
exports <- ls("package:anndataR")
cat("\nAll anndataR exports:\n")
cat(paste(exports, collapse="\n"), "\n")

# Check for conversion functions
cat("\nConversion-related functions:\n")
conv_funcs <- grep("AnnData|anndata|seurat|Seurat|convert|from_|to_|as_", exports, value=TRUE)
cat(paste(conv_funcs, collapse="\n"), "\n")

# Test with smallest dataset
cat("\n--- Testing with pbmc3k ---\n")
obj <- readRDS("/benchmark/data/generated/seurat_v5_pbmc3k_raw.rds")
tryCatch({obj <- UpdateSeuratObject(obj)}, error=function(e){})

# Try as_AnnData
cat("\nTrying as_AnnData(obj)...\n")
tryCatch({
  ad <- as_AnnData(obj)
  cat("SUCCESS: as_AnnData\n")
}, error=function(e) {
  cat("ERROR as_AnnData:", e$message, "\n")
})

# Try from_Seurat
cat("\nTrying from_Seurat(obj)...\n")
tryCatch({
  ad <- from_Seurat(obj)
  cat("SUCCESS: from_Seurat\n")
}, error=function(e) {
  cat("ERROR from_Seurat:", e$message, "\n")
})

# Try AnnData() constructor
cat("\nTrying AnnData() constructor...\n")
tryCatch({
  ad <- AnnData(obj)
  cat("SUCCESS: AnnData()\n")
}, error=function(e) {
  cat("ERROR AnnData():", e$message, "\n")
})

# Try write_h5ad if we have an ad object
cat("\nTrying to create AnnData manually from Seurat layers...\n")
tryCatch({
  counts <- GetAssayData(obj, layer="counts")
  ad <- AnnData(X=counts, obs=obj@meta.data)
  cat("SUCCESS: manual AnnData creation\n")
  cat("Class:", class(ad), "\n")
  
  out <- "/tmp/test_anndataR_output.h5ad"
  write_h5ad(ad, out)
  cat("SUCCESS: write_h5ad\n")
}, error=function(e) {
  cat("ERROR manual:", e$message, "\n")
})
