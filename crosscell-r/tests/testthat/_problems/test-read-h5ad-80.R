# Extracted from test-read-h5ad.R:80

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "crosscell", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(crosscell)
skip_if_no_test_data <- function() {
  h5ad_path <- "../../tests/data/small_sparse.h5ad"
  if (!file.exists(h5ad_path)) {
    skip("Test H5AD file not found")
  }
}

# test -------------------------------------------------------------------------
tmp <- tempfile(fileext = ".rds")
file.create(tmp)
on.exit(unlink(tmp))
expect_error(
    read_h5ad_as_seurat(tmp),
    "Expected .h5ad"
  )
