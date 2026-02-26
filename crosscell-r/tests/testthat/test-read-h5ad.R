# Test read_h5ad_as_seurat and read_h5ad_as_sce functions
# Task 13.1: Integration tests for H5AD reading

library(testthat)
library(crosscell)

# Skip if test data not available
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

test_that("read_h5ad_as_seurat returns Seurat object", {
  h5ad_path <- get_test_h5ad_path()
  skip_if(is.null(h5ad_path), "Test H5AD file not found")
  
  result <- read_h5ad_as_seurat(h5ad_path)
  
  expect_s4_class(result, "Seurat")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})

test_that("read_h5ad_as_sce returns SingleCellExperiment object", {
  h5ad_path <- get_test_h5ad_path()
  skip_if(is.null(h5ad_path), "Test H5AD file not found")
  
  result <- read_h5ad_as_sce(h5ad_path)
  
  expect_s4_class(result, "SingleCellExperiment")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})

test_that("read_h5ad_as_seurat preserves dimensions", {
  h5ad_path <- get_test_h5ad_path()
  skip_if(is.null(h5ad_path), "Test H5AD file not found")
  
  info <- inspect_file(h5ad_path)
  result <- read_h5ad_as_seurat(h5ad_path)
  
  expect_equal(ncol(result), info$n_cells)
  expect_equal(nrow(result), info$n_genes)
})

test_that("read_h5ad_as_sce preserves dimensions", {
  h5ad_path <- get_test_h5ad_path()
  skip_if(is.null(h5ad_path), "Test H5AD file not found")
  
  info <- inspect_file(h5ad_path)
  result <- read_h5ad_as_sce(h5ad_path)
  
  expect_equal(ncol(result), info$n_cells)
  expect_equal(nrow(result), info$n_genes)
})

test_that("read_h5ad handles file not found", {
  expect_error(
    read_h5ad_as_seurat("nonexistent.h5ad"),
    "File not found"
  )
  
  expect_error(
    read_h5ad_as_sce("nonexistent.h5ad"),
    "File not found"
  )
})

test_that("read_h5ad handles wrong extension", {
  # Create a temp file with wrong extension
  tmp <- tempfile(fileext = ".rds")
  file.create(tmp)
  on.exit(unlink(tmp))
  
  # extendr wraps errors in panic messages
  expect_error(
    read_h5ad_as_seurat(tmp),
    "Expected .h5ad|panicked"
  )
})
