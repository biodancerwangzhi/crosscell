# Test error handling
# Task 13.5: Error handling tests

library(testthat)
library(crosscell)

test_that("read_rds_fast handles file not found", {
  expect_error(
    read_rds_fast("nonexistent.rds"),
    "File not found"
  )
})

test_that("read_rds_fast handles wrong extension", {
  tmp <- tempfile(fileext = ".h5ad")
  file.create(tmp)
  on.exit(unlink(tmp))
  
  # extendr wraps errors in panic messages, so we check for either format
  expect_error(
    read_rds_fast(tmp),
    "Expected .rds|panicked"
  )
})

test_that("write_as_h5ad handles wrong object type", {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  
  expect_error(
    write_as_h5ad(data.frame(x = 1:10), tmp),
    "Seurat or SingleCellExperiment"
  )
})

test_that("write_as_h5ad handles wrong extension", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  library(Matrix)
  
  # Create minimal Seurat object
  counts <- matrix(1:100, nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  obj <- CreateSeuratObject(counts = as(counts, "dgCMatrix"))
  
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  
  # extendr wraps errors in panic messages
  expect_error(
    write_as_h5ad(obj, tmp),
    "Expected .h5ad|panicked"
  )
})

test_that("write_rds_fast handles wrong object type", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  
  expect_error(
    write_rds_fast(data.frame(x = 1:10), tmp),
    "Seurat object"
  )
})

test_that("write_rds_fast handles wrong extension", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  library(Matrix)
  
  counts <- matrix(1:100, nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  obj <- CreateSeuratObject(counts = as(counts, "dgCMatrix"))
  
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  
  # extendr wraps errors in panic messages
  expect_error(
    write_rds_fast(obj, tmp),
    "Expected .rds|panicked"
  )
})
