#!/usr/bin/env Rscript
# Test anndataR::read_h5ad() on all 13 H5AD datasets
# Purpose: Find which datasets anndataR can/cannot read, with error details
# Usage: docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark-limited Rscript /benchmark/scripts/test_anndataR_h5ad_read.R

suppressPackageStartupMessages({
  library(anndataR)
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

cat("=== anndataR read_h5ad() Test ===\n")
cat(sprintf("anndataR version: %s\n", packageVersion("anndataR")))
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

  # Step 1: read_h5ad
  read_ok <- FALSE
  read_err <- ""
  read_warn <- ""
  adata <- NULL

  tryCatch({
    withCallingHandlers(
      { adata <- read_h5ad(path) },
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
    cat(sprintf("❌  %-35s read FAILED: %s\n", name, substr(read_err, 1, 100)))
    results[[name]] <- list(read = "failed", to_seurat = NA, error = read_err)
    next
  }

  # Report read success (with or without warnings)
  n_obs <- tryCatch(adata$n_obs(), error = function(e) NA)
  n_var <- tryCatch(adata$n_vars(), error = function(e) NA)
  warn_tag <- if (nchar(read_warn) > 0) " ⚠️WARN" else ""

  # Step 2: to_Seurat
  seurat_ok <- FALSE
  seurat_err <- ""

  tryCatch({
    obj <- adata$to_Seurat()
    seurat_ok <- TRUE
  }, error = function(e) {
    seurat_err <<- conditionMessage(e)
  })

  if (seurat_ok) {
    cat(sprintf("✅  %-35s read OK (%s×%s)%s → to_Seurat OK\n",
                name, n_obs, n_var, warn_tag))
    results[[name]] <- list(read = "success", to_seurat = "success",
                            n_obs = n_obs, n_var = n_var,
                            read_warning = read_warn, error = "")
  } else {
    cat(sprintf("🔶  %-35s read OK (%s×%s)%s → to_Seurat FAILED: %s\n",
                name, n_obs, n_var, warn_tag, substr(seurat_err, 1, 100)))
    results[[name]] <- list(read = "success", to_seurat = "failed",
                            n_obs = n_obs, n_var = n_var,
                            read_warning = read_warn, error = seurat_err)
  }
}

# Summary
cat("\n=== Summary ===\n")
read_ok_count <- sum(sapply(results, function(r) r$read == "success"))
seurat_ok_count <- sum(sapply(results, function(r) identical(r$to_seurat, "success")))
cat(sprintf("read_h5ad:  %d/%d success\n", read_ok_count, length(datasets)))
cat(sprintf("to_Seurat:  %d/%d success\n", seurat_ok_count, length(datasets)))

# Save JSON
output_file <- "/benchmark/results/anndataR_h5ad_read_test.json"
tryCatch({
  jsonlite::write_json(results, output_file, pretty = TRUE, auto_unbox = TRUE)
  cat(sprintf("\nResults saved to %s\n", output_file))
}, error = function(e) {
  cat(sprintf("\nCould not save JSON: %s\n", e$message))
})
