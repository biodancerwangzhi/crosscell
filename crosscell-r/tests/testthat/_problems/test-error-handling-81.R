# Extracted from test-error-handling.R:81

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "crosscell", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(crosscell)

# test -------------------------------------------------------------------------
skip_if_not_installed("Seurat")
library(Seurat)
library(Matrix)
counts <- matrix(1:100, nrow = 10, ncol = 10)
rownames(counts) <- paste0("Gene", 1:10)
colnames(counts) <- paste0("Cell", 1:10)
obj <- CreateSeuratObject(counts = as(counts, "dgCMatrix"))
tmp <- tempfile(fileext = ".h5ad")
on.exit(unlink(tmp))
expect_error(
    write_rds_fast(obj, tmp),
    "Expected .rds"
  )
