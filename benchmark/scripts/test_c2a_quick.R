#!/usr/bin/env Rscript
# Quick convert2anndata test: only test convert_seurat_to_sce step (no Python needed)
suppressPackageStartupMessages({
  library(convert2anndata)
  library(Seurat)
})

data_dir <- "/benchmark/data/generated"

tests <- list(
  list(file = "seurat_v4_pbmc3k_raw.rds", id = "v4_pbmc3k_raw"),
  list(file = "seurat_v4_ifnb_raw.rds", id = "v4_ifnb_raw"),
  list(file = "seurat_v5_pbmc3k_raw.rds", id = "v5_pbmc3k_raw"),
  list(file = "seurat_v5_ssHippo_raw.rds", id = "v5_ssHippo_raw"),
  list(file = "seurat_v5_stxKidney_raw.rds", id = "v5_stxKidney_raw"),
  list(file = "seurat_v5_thp1.eccite_raw.rds", id = "v5_thp1.eccite_raw")
)

cat("=== convert2anndata: convert_seurat_to_sce step only ===\n")
cat(sprintf("convert2anndata: %s, Seurat: %s\n\n", packageVersion("convert2anndata"), packageVersion("Seurat")))

for (tc in tests) {
  path <- file.path(data_dir, tc$file)
  if (!file.exists(path)) { cat(sprintf("SKIP %s\n", tc$id)); next }

  obj <- tryCatch({
    o <- readRDS(path)
    tryCatch(UpdateSeuratObject(o), error = function(e) o)
  }, error = function(e) NULL)
  if (is.null(obj)) { cat(sprintf("❌ %s: readRDS failed\n", tc$id)); next }

  assay_class <- tryCatch(class(obj[[DefaultAssay(obj)]])[1], error = function(e) "unknown")
  assays <- paste(Assays(obj), collapse=",")

  sce <- tryCatch({
    convert_seurat_to_sce(obj)
  }, error = function(e) {
    cat(sprintf("❌ %-30s [%s] assays=%s\n   Error: %s\n", tc$id, assay_class, assays, e$message))
    NULL
  })

  if (!is.null(sce)) {
    cat(sprintf("✅ %-30s [%s] assays=%s → SCE OK (%d×%d)\n",
                tc$id, assay_class, assays, ncol(sce), nrow(sce)))
  }
}
