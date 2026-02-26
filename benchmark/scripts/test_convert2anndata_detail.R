#!/usr/bin/env Rscript
# Detailed convert2anndata failure analysis
# Tests representative V4, V5, and spatial datasets with full error output
# Usage: docker-compose -f docker-compose.benchmark.yml run --rm benchmark-limited Rscript /benchmark/scripts/test_convert2anndata_detail.R

suppressPackageStartupMessages({
  library(convert2anndata)
  library(Seurat)
  library(anndata)
})

data_dir <- "/benchmark/data/generated"

cat("=== convert2anndata Detailed Failure Analysis ===\n")
cat(sprintf("convert2anndata version: %s\n", packageVersion("convert2anndata")))
cat(sprintf("Seurat version: %s\n", packageVersion("Seurat")))
cat(sprintf("anndata (R) version: %s\n", packageVersion("anndata")))
cat("\n")

# Representative test cases
tests <- list(
  # V4 (all fail)
  list(file = "seurat_v4_pbmc3k_raw.rds", id = "v4_pbmc3k_raw"),
  list(file = "seurat_v4_cbmc_raw.rds", id = "v4_cbmc_raw"),
  # V5 (most succeed)
  list(file = "seurat_v5_pbmc3k_raw.rds", id = "v5_pbmc3k_raw"),
  list(file = "seurat_v5_cbmc_raw.rds", id = "v5_cbmc_raw"),
  # V5 spatial (fails)
  list(file = "seurat_v5_ssHippo_raw.rds", id = "v5_ssHippo_raw"),
  list(file = "seurat_v5_stxKidney_raw.rds", id = "v5_stxKidney_raw")
)

for (tc in tests) {
  path <- file.path(data_dir, tc$file)
  if (!file.exists(path)) {
    cat(sprintf("SKIP %-35s NOT FOUND\n", tc$id))
    next
  }

  cat(sprintf("\n=== %s ===\n", tc$id))

  # Step 1: Load RDS
  obj <- tryCatch({
    o <- readRDS(path)
    tryCatch(UpdateSeuratObject(o), error = function(e) o)
  }, error = function(e) {
    cat(sprintf("  ❌ readRDS failed: %s\n", e$message))
    NULL
  })
  if (is.null(obj)) next

  cat(sprintf("  Seurat version: %s\n", Version(obj)))
  cat(sprintf("  Assays: %s\n", paste(Assays(obj), collapse = ", ")))
  cat(sprintf("  Cells: %d, Genes: %d\n", ncol(obj), nrow(obj)))
  cat(sprintf("  Default assay class: %s\n", class(obj[["RNA"]])[1]))

  # Step 2: convert_seurat_to_sce
  sce <- tryCatch({
    convert_seurat_to_sce(obj)
  }, error = function(e) {
    cat(sprintf("  ❌ convert_seurat_to_sce failed: %s\n", e$message))
    NULL
  })
  if (is.null(sce)) next
  cat(sprintf("  ✅ convert_seurat_to_sce OK\n"))

  # Step 3: convert_to_anndata
  ad <- tryCatch({
    convert_to_anndata(sce)
  }, error = function(e) {
    cat(sprintf("  ❌ convert_to_anndata failed: %s\n", e$message))
    NULL
  })
  if (is.null(ad)) next
  cat(sprintf("  ✅ convert_to_anndata OK\n"))

  # Step 4: write_h5ad
  out <- file.path("/tmp", paste0(tc$id, "_c2a.h5ad"))
  tryCatch({
    anndata::write_h5ad(ad, out)
    cat(sprintf("  ✅ write_h5ad OK → %s\n", out))
  }, error = function(e) {
    cat(sprintf("  ❌ write_h5ad failed: %s\n", e$message))
  })
}

cat("\n=== Summary ===\n")
cat("V4 datasets: All fail at convert_seurat_to_sce (Assay vs Assay5 class)\n")
cat("V5 datasets: Most succeed\n")
cat("V5 spatial: May fail at convert_to_anndata (spatial assay handling)\n")
