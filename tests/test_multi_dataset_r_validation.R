#!/usr/bin/env Rscript
# Multi-dataset R validation for CrossCell-generated Seurat RDS files
#
# Usage: docker-compose run --rm dev Rscript tests/test_multi_dataset_r_validation.R
#
# Prerequisites: Run the Rust test first to generate the RDS files:
#   docker-compose run --rm dev cargo test --test test_multi_dataset_r_validation -- --nocapture

library(SeuratObject)

rds_dir <- "tests/data/multi_val"
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  cat("ERROR: No RDS files found in", rds_dir, "\n")
  cat("Run the Rust test first:\n")
  cat("  cargo test --test test_multi_dataset_r_validation -- --nocapture\n")
  quit(status = 1)
}

cat("Found", length(rds_files), "RDS files to validate\n\n")

total <- 0
passed <- 0
failed <- 0
errors <- list()

validate_seurat <- function(path) {
  name <- basename(path)
  cat("--- ", name, " ---\n")

  # 1. readRDS
  obj <- tryCatch(readRDS(path), error = function(e) {
    cat("  FAIL readRDS:", e$message, "\n")
    return(NULL)
  })
  if (is.null(obj)) return(FALSE)

  cat("  class:", paste(class(obj), collapse = ", "), "\n")
  cat("  isS4:", isS4(obj), "\n")

  if (!isS4(obj) || !is(obj, "Seurat")) {
    cat("  FAIL: not a valid Seurat S4 object\n")
    return(FALSE)
  }

  ok <- TRUE

  # 2. ncol / nrow
  nc <- tryCatch(ncol(obj), error = function(e) { cat("  FAIL ncol:", e$message, "\n"); ok <<- FALSE; NA })
  nr <- tryCatch(nrow(obj), error = function(e) { cat("  FAIL nrow:", e$message, "\n"); ok <<- FALSE; NA })
  cat("  dim:", nr, "genes x", nc, "cells\n")

  # 3. active.assay
  aa <- tryCatch(obj@active.assay, error = function(e) { cat("  FAIL active.assay:", e$message, "\n"); ok <<- FALSE; NA })
  cat("  active.assay:", aa, "\n")

  # 4. GetAssayData
  tryCatch({
    counts <- GetAssayData(obj, layer = "counts")
    cat("  GetAssayData: OK (", class(counts), paste(dim(counts), collapse = "x"), ")\n")
  }, error = function(e) {
    cat("  FAIL GetAssayData:", e$message, "\n")
    ok <<- FALSE
  })

  # 5. meta.data
  tryCatch({
    md <- obj@meta.data
    cat("  meta.data: OK (", nrow(md), "x", ncol(md), ")\n")
  }, error = function(e) {
    cat("  FAIL meta.data:", e$message, "\n")
    ok <<- FALSE
  })

  # 6. Cells()
  tryCatch({
    cn <- Cells(obj)
    cat("  Cells(): OK (", length(cn), ")\n")
  }, error = function(e) {
    cat("  FAIL Cells():", e$message, "\n")
    ok <<- FALSE
  })

  # 7. Features()
  tryCatch({
    fn <- Features(obj)
    cat("  Features(): OK (", length(fn), ")\n")
  }, error = function(e) {
    cat("  FAIL Features():", e$message, "\n")
    ok <<- FALSE
  })

  # 8. Idents()
  tryCatch({
    id <- Idents(obj)
    cat("  Idents(): OK (", length(id), "idents)\n")
  }, error = function(e) {
    cat("  FAIL Idents():", e$message, "\n")
    ok <<- FALSE
  })

  # 9. Subset
  if (!is.na(nc) && nc >= 3) {
    tryCatch({
      sub <- obj[, 1:min(3, nc)]
      cat("  Subset: OK (", ncol(sub), "cells)\n")
    }, error = function(e) {
      cat("  FAIL Subset:", e$message, "\n")
      ok <<- FALSE
    })
  }

  # 10. Reductions
  tryCatch({
    reds <- Reductions(obj)
    if (length(reds) > 0) {
      cat("  Reductions:", paste(reds, collapse = ", "), "\n")
      for (r in reds) {
        emb <- Embeddings(obj, reduction = r)
        cat("    ", r, ":", paste(dim(emb), collapse = "x"), "\n")
      }
    } else {
      cat("  Reductions: none\n")
    }
  }, error = function(e) {
    cat("  WARN Reductions:", e$message, "\n")
  })

  # 11. Expression values sanity check
  tryCatch({
    counts <- GetAssayData(obj, layer = "counts")
    nnz <- sum(counts != 0)
    total_el <- prod(dim(counts))
    sparsity <- 1 - nnz / total_el
    cat("  Expression: nnz=", nnz, " sparsity=", round(sparsity * 100, 1), "%\n")
    if (any(is.nan(counts@x)) || any(is.infinite(counts@x))) {
      cat("  WARN: NaN or Inf values in expression matrix\n")
    }
  }, error = function(e) {
    # Dense matrix case
    tryCatch({
      counts <- GetAssayData(obj, layer = "counts")
      cat("  Expression: dense matrix", paste(dim(counts), collapse = "x"), "\n")
    }, error = function(e2) {
      cat("  WARN expression check:", e2$message, "\n")
    })
  })

  cat("  Result:", ifelse(ok, "PASS", "FAIL"), "\n\n")
  return(ok)
}

for (f in rds_files) {
  total <- total + 1
  result <- tryCatch(
    validate_seurat(f),
    error = function(e) {
      cat("  UNEXPECTED ERROR:", e$message, "\n\n")
      errors[[basename(f)]] <<- e$message
      return(FALSE)
    }
  )
  if (isTRUE(result)) {
    passed <- passed + 1
  } else {
    failed <- failed + 1
    if (is.null(errors[[basename(f)]])) {
      errors[[basename(f)]] <- "validation failed"
    }
  }
}

cat("========================================\n")
cat("TOTAL:", total, " PASSED:", passed, " FAILED:", failed, "\n")
if (failed > 0) {
  cat("\nFailed datasets:\n")
  for (name in names(errors)) {
    cat("  -", name, ":", errors[[name]], "\n")
  }
  quit(status = 1)
} else {
  cat("All datasets passed R validation!\n")
  quit(status = 0)
}
