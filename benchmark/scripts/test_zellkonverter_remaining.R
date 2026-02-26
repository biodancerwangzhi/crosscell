#!/usr/bin/env Rscript
# Test Zellkonverter on remaining 7 H5AD files (the first 6 already confirmed OK)
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Seurat)
})

data_dir <- "/benchmark/data/generated"
files <- c(
  "cellxgene_brain_40k.h5ad",
  "squidpy_visium_hne.h5ad",
  "squidpy_mibitof.h5ad",
  "squidpy_slideseqv2.h5ad",
  "squidpy_seqfish.h5ad",
  "squidpy_merfish.h5ad",
  "squidpy_imc.h5ad"
)

cat("=== Zellkonverter remaining 7 files ===\n\n")

for (f in files) {
  path <- file.path(data_dir, f)
  name <- sub("\\.h5ad$", "", f)

  if (!file.exists(path)) {
    cat(sprintf("SKIP %-35s NOT FOUND\n", name))
    next
  }

  # Step 1: readH5AD
  sce <- NULL
  err1 <- ""
  tryCatch({
    sce <- readH5AD(path)
  }, error = function(e) {
    err1 <<- e$message
  })

  if (is.null(sce)) {
    cat(sprintf("❌  %-35s readH5AD FAILED: %s\n", name, substr(err1, 1, 100)))
    next
  }

  n_obs <- ncol(sce)
  n_var <- nrow(sce)
  assays <- assayNames(sce)

  # Step 2: as.Seurat with explicit counts
  err2 <- ""
  tryCatch({
    obj <- as.Seurat(sce, counts = assays[1], data = NULL)
  }, error = function(e) {
    err2 <<- e$message
  })

  if (nchar(err2) == 0) {
    cat(sprintf("✅  %-35s OK (%d×%d, assays:%s)\n", name, n_obs, n_var, paste(assays, collapse=",")))
  } else {
    cat(sprintf("🔶  %-35s readH5AD OK (%d×%d) → as.Seurat FAILED: %s\n",
                name, n_obs, n_var, substr(err2, 1, 100)))
  }
}
