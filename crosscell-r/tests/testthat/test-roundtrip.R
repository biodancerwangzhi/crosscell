# Test roundtrip conversions
# Task 13.2 & 13.3: Roundtrip tests for H5AD and RDS

library(testthat)
library(crosscell)
library(Seurat)
library(Matrix)

# Helper to create a simple Seurat object for testing
create_test_seurat <- function(n_genes = 50, n_cells = 100) {
  set.seed(42)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Cell", 1:n_cells)
  
  sparse_counts <- as(counts, "dgCMatrix")
  
  obj <- CreateSeuratObject(counts = sparse_counts, project = "test")
  
  # Add some metadata
  obj$cell_type <- sample(c("A", "B", "C"), n_cells, replace = TRUE)
  obj$batch <- sample(1:3, n_cells, replace = TRUE)
  
  obj
}

test_that("Seurat H5AD roundtrip preserves dimensions", {
  skip_on_cran()
  
  obj <- create_test_seurat()
  tmp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_h5ad))
  
  # Write to H5AD
  write_as_h5ad(obj, tmp_h5ad)
  expect_true(file.exists(tmp_h5ad))
  
  # Read back
  obj2 <- read_h5ad_as_seurat(tmp_h5ad)
  
  # Check dimensions
  expect_equal(ncol(obj2), ncol(obj))
  expect_equal(nrow(obj2), nrow(obj))
})

test_that("Seurat H5AD roundtrip preserves expression data", {
  skip_on_cran()
  
  obj <- create_test_seurat(n_genes = 20, n_cells = 30)
  tmp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_h5ad))
  
  write_as_h5ad(obj, tmp_h5ad)
  obj2 <- read_h5ad_as_seurat(tmp_h5ad)
  
  # Get counts matrices
  counts1 <- GetAssayData(obj, layer = "counts")
  counts2 <- GetAssayData(obj2, layer = "counts")
  
  # Check correlation (should be very high)
  cor_val <- cor(as.vector(as.matrix(counts1)), as.vector(as.matrix(counts2)))
  expect_gt(cor_val, 0.999)
})

test_that("Seurat RDS roundtrip preserves dimensions", {
  skip_on_cran()
  
  obj <- create_test_seurat()
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds))
  
  # Write to RDS
  write_rds_fast(obj, tmp_rds)
  expect_true(file.exists(tmp_rds))
  
  # Read back
  obj2 <- read_rds_fast(tmp_rds)
  
  # Check dimensions
  expect_equal(ncol(obj2), ncol(obj))
  expect_equal(nrow(obj2), nrow(obj))
})

test_that("SCE H5AD roundtrip preserves dimensions", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")
  
  library(SingleCellExperiment)
  
  # Create SCE
  set.seed(42)
  n_genes <- 50
  n_cells <- 100
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Cell", 1:n_cells)
  
  sce <- SingleCellExperiment(
    assays = list(counts = as(counts, "dgCMatrix")),
    colData = data.frame(cell_type = sample(c("A", "B"), n_cells, replace = TRUE))
  )
  
  tmp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_h5ad))
  
  # Write to H5AD
  write_as_h5ad(sce, tmp_h5ad)
  expect_true(file.exists(tmp_h5ad))
  
  # Read back
  sce2 <- read_h5ad_as_sce(tmp_h5ad)
  
  # Check dimensions
  expect_equal(ncol(sce2), ncol(sce))
  expect_equal(nrow(sce2), nrow(sce))
})
