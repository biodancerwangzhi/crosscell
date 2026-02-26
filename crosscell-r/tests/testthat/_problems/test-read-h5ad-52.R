# Extracted from test-read-h5ad.R:52

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "crosscell", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(crosscell)
skip_if_no_test_data <- function() {
  # Try multiple possible paths
  paths <- c(
    "../../tests/data/small_sparse.h5ad",
    "../../../tests/data/small_sparse.h5ad",
    "tests/data/small_sparse.h5ad"
  )
  for (p in paths) {
    if (file.exists(p)) {
      return(p)
    }
  }
  skip("Test H5AD file not found")
}
get_test_h5ad_path <- function() {
  paths <- c(
    "../../tests/data/small_sparse.h5ad",
    "../../../tests/data/small_sparse.h5ad",
    "tests/data/small_sparse.h5ad"
  )
  for (p in paths) {
    if (file.exists(p)) {
      return(p)
    }
  }
  NULL
}

# test -------------------------------------------------------------------------
h5ad_path <- get_test_h5ad_path()
skip_if(is.null(h5ad_path), "Test H5AD file not found")
result <- read_h5ad_as_sce(h5ad_path)
