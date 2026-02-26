# Test inspect_file function
# Task 13.5: Error handling tests

library(testthat)
library(crosscell)

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

get_test_rds_path <- function() {
  paths <- c(
    "../../data/generated/seurat_v4_pbmc3k_raw.rds",
    "../../../data/generated/seurat_v4_pbmc3k_raw.rds",
    "data/generated/seurat_v4_pbmc3k_raw.rds"
  )
  for (p in paths) {
    if (file.exists(p)) {
      return(p)
    }
  }
  NULL
}

test_that("inspect_file returns correct structure for H5AD", {
  h5ad_path <- get_test_h5ad_path()
  skip_if(is.null(h5ad_path), "Test H5AD file not found")
  
  info <- inspect_file(h5ad_path)
  
  expect_type(info, "list")
  expect_equal(info$format, "h5ad")
  expect_type(info$n_cells, "integer")
  expect_type(info$n_genes, "integer")
  expect_true(info$n_cells > 0)
  expect_true(info$n_genes > 0)
})

test_that("inspect_file returns correct structure for RDS", {
  rds_path <- get_test_rds_path()
  skip_if(is.null(rds_path), "Test RDS file not found")
  
  info <- inspect_file(rds_path)
  
  expect_type(info, "list")
  expect_equal(info$format, "rds")
  expect_type(info$n_cells, "integer")
  expect_type(info$n_genes, "integer")
})

test_that("inspect_file handles file not found", {
  expect_error(
    inspect_file("nonexistent.h5ad"),
    "File not found"
  )
})

test_that("inspect_file handles unsupported format", {
  tmp <- tempfile(fileext = ".csv")
  file.create(tmp)
  on.exit(unlink(tmp))
  
  expect_error(
    inspect_file(tmp),
    "Unsupported file format"
  )
})
