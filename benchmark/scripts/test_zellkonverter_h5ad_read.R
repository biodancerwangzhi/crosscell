#!/usr/bin/env Rscript
# Test Zellkonverter readH5AD() + as.Seurat() on all 13 H5AD datasets
# Purpose: Find which datasets Zellkonverter can read and convert to Seurat
# Usage: docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark-limited Rscript /benchmark/scripts/test_zellkonverter_h5ad_read.R

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Seurat)
})

data_dir <- "/benchmark/data/generated"

datasets <- c(
  "scanpy_pbmc3k.h5ad",
  "scanpy_pbmc3k_processed.h5ad",
  "scvelo_dentategyrus.h5ad",
  "scvelo_pancreas.h5ad",
  "cellxgene_pbmc_15k.h5ad",
  "cellxgene_heart_23k.h5ad",
  "cellxgene_brain_40k.h5ad",
  "squidpy_visium_hne.h5ad",
  "squidpy_mibitof.h5ad",
  "squidpy_slideseqv2.h5ad",
  "squidpy_seqfish.h5ad",
  "squidpy_merfish.h5ad",
  "squidpy_imc.h5ad"
)

cat("=== Zellkonverter readH5AD() + as.Seurat() Test ===\n")
cat(sprintf("Zellkonverter version: %s\n", packageVersion("zellkonverter")))
cat(sprintf("Seurat version: %s\n", packageVersion("Seurat")))
cat(sprintf("Testing %d H5AD files\n\n", length(datasets)))

results <- list()

for (f in datasets) {
  path <- file.path(data_dir, f)
  name <- sub("\\.h5ad$", "", f)

  if (!file.exists(path)) {
    cat(sprintf("⚠️  %-35s NOT FOUND\n", name))
    results[[name]] <- list(read = "not_found", to_seurat = NA, error = "file not found")
    next
  }

  # Step 1: readH5AD → SCE
  read_ok <- FALSE
  read_err <- ""
  read_warn <- ""
  sce <- NULL

  tryCatch({
    withCallingHandlers(
      { sce <- readH5AD(path) },
      warning = function(w) {
        read_warn <<- paste0(read_warn, conditionMessage(w), "; ")
        invokeRestart("muffleWarning")
      }
    )
    read_ok <- TRUE
  }, error = function(e) {
    read_err <<- conditionMessage(e)
  })

  if (!read_ok) {
    cat(sprintf("❌  %-35s readH5AD FAILED: %s\n", name, substr(read_err, 1, 100)))
    results[[name]] <- list(read = "failed", to_seurat = NA, error = read_err)
    next
  }

  n_obs <- ncol(sce)
  n_var <- nrow(sce)
  assay_names <- assayNames(sce)
  warn_tag <- if (nchar(read_warn) > 0) " ⚠️WARN" else ""

  # Step 2: as.Seurat
  seurat_ok <- FALSE
  seurat_err <- ""

  tryCatch({
    obj <- as.Seurat(sce)
    seurat_ok <- TRUE
  }, error = function(e) {
    seurat_err <<- conditionMessage(e)
  })

  # Step 2b: If as.Seurat failed, try with explicit assay
  if (!seurat_ok && length(assay_names) > 0) {
    tryCatch({
      obj <- as.Seurat(sce, counts = assay_names[1], data = NULL)
      seurat_ok <- TRUE
      seurat_err <- paste0("(fixed with counts=", assay_names[1], ")")
    }, error = function(e) {
      # keep original error
    })
  }

  if (seurat_ok) {
    cat(sprintf("✅  %-35s readH5AD OK (%d×%d, assays:%s)%s → as.Seurat OK %s\n",
                name, n_obs, n_var, paste(assay_names, collapse=","), warn_tag, seurat_err))
    results[[name]] <- list(read = "success", to_seurat = "success",
                            n_obs = n_obs, n_var = n_var,
                            assays = assay_names,
                            read_warning = read_warn, error = seurat_err)
  } else {
    cat(sprintf("🔶  %-35s readH5AD OK (%d×%d, assays:%s)%s → as.Seurat FAILED: %s\n",
                name, n_obs, n_var, paste(assay_names, collapse=","), warn_tag,
                substr(seurat_err, 1, 100)))
    results[[name]] <- list(read = "success", to_seurat = "failed",
                            n_obs = n_obs, n_var = n_var,
                            assays = assay_names,
                            read_warning = read_warn, error = seurat_err)
  }
}

# Summary
cat("\n=== Summary ===\n")
read_ok_count <- sum(sapply(results, function(r) r$read == "success"))
seurat_ok_count <- sum(sapply(results, function(r) identical(r$to_seurat, "success")))
cat(sprintf("readH5AD:   %d/%d success\n", read_ok_count, length(datasets)))
cat(sprintf("as.Seurat:  %d/%d success\n", seurat_ok_count, length(datasets)))

# Save JSON
output_file <- "/benchmark/results/zellkonverter_h5ad_read_test.json"
tryCatch({
  jsonlite::write_json(results, output_file, pretty = TRUE, auto_unbox = TRUE)
  cat(sprintf("\nResults saved to %s\n", output_file))
}, error = function(e) {
  cat(sprintf("\nCould not save JSON: %s\n", e$message))
})
