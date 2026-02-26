# Test core R logic without crosscell package
# This script validates R object structures and conversions

library(Seurat)
library(SingleCellExperiment)
library(Matrix)

cat("=== Testing Core R Logic for CrossCell ===\n\n")

# ============================================================================
# Test 1: dgCMatrix structure (used in dgcmatrix_to_expression)
# ============================================================================
cat("1. Testing dgCMatrix structure extraction:\n")

set.seed(42)
n_genes <- 10
n_cells <- 5
counts <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Cell", 1:n_cells)
sparse <- as(counts, "dgCMatrix")

# Verify slot access (same as Rust code)
cat("   @i (row indices) length:", length(sparse@i), "\n")
cat("   @p (col pointers) length:", length(sparse@p), "\n")
cat("   @x (data) length:", length(sparse@x), "\n")
cat("   @Dim:", sparse@Dim, "\n")

# Verify CSC format
stopifnot(length(sparse@p) == n_cells + 1)
stopifnot(sparse@Dim[1] == n_genes)
stopifnot(sparse@Dim[2] == n_cells)
cat("   ✓ dgCMatrix structure verified\n\n")

# ============================================================================
# Test 2: Seurat object structure (used in seurat_to_ir)
# ============================================================================
cat("2. Testing Seurat object structure:\n")

seurat_obj <- CreateSeuratObject(counts = sparse, project = "test")
seurat_obj$cell_type <- factor(sample(c("A", "B"), n_cells, replace = TRUE))
seurat_obj$score <- runif(n_cells)

# Verify access patterns (same as Rust code)
cat("   DefaultAssay:", DefaultAssay(seurat_obj), "\n")
cat("   ncol:", ncol(seurat_obj), "\n")
cat("   nrow:", nrow(seurat_obj), "\n")
cat("   meta.data columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
cat("   reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")

# Test GetAssayData
counts_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
cat("   counts class:", class(counts_data), "\n")
stopifnot(inherits(counts_data, "dgCMatrix"))
cat("   ✓ Seurat structure verified\n\n")

# ============================================================================
# Test 3: Seurat with embeddings (used in ir_to_seurat)
# ============================================================================
cat("3. Testing Seurat with embeddings:\n")

# Add PCA
pca_emb <- matrix(rnorm(n_cells * 5), nrow = n_cells, ncol = 5)
rownames(pca_emb) <- colnames(seurat_obj)
colnames(pca_emb) <- paste0("PC_", 1:5)
seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = pca_emb, key = "PC_", assay = "RNA")

cat("   reductions after PCA:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")
cat("   PCA embeddings dim:", dim(seurat_obj@reductions$pca@cell.embeddings), "\n")

# Verify DimReduc access
dr <- seurat_obj@reductions$pca
stopifnot(inherits(dr, "DimReduc"))
stopifnot(all(dim(dr@cell.embeddings) == c(n_cells, 5)))
cat("   ✓ Seurat embeddings verified\n\n")

# ============================================================================
# Test 4: SingleCellExperiment structure (used in sce_to_ir)
# ============================================================================
cat("4. Testing SingleCellExperiment structure:\n")

sce <- SingleCellExperiment(
  assays = list(counts = sparse),
  colData = data.frame(
    cell_type = factor(sample(c("A", "B"), n_cells, replace = TRUE)),
    score = runif(n_cells)
  ),
  rowData = data.frame(
    gene_biotype = sample(c("protein_coding", "lncRNA"), n_genes, replace = TRUE)
  )
)

# Add reducedDims
reducedDim(sce, "PCA") <- pca_emb

cat("   ncol:", ncol(sce), "\n")
cat("   nrow:", nrow(sce), "\n")
cat("   assayNames:", paste(assayNames(sce), collapse = ", "), "\n")
cat("   colData columns:", paste(colnames(colData(sce)), collapse = ", "), "\n")
cat("   rowData columns:", paste(colnames(rowData(sce)), collapse = ", "), "\n")
cat("   reducedDimNames:", paste(reducedDimNames(sce), collapse = ", "), "\n")

# Verify access patterns
stopifnot(inherits(assay(sce, "counts"), "dgCMatrix"))
stopifnot(inherits(as.data.frame(colData(sce)), "data.frame"))
cat("   ✓ SCE structure verified\n\n")

# ============================================================================
# Test 5: data.frame with factor (used in r_to_dataframe)
# ============================================================================
cat("5. Testing data.frame with factor:\n")

test_df <- data.frame(
  int_col = 1:5,
  dbl_col = c(1.1, 2.2, 3.3, 4.4, 5.5),
  str_col = c("a", "b", "c", "d", "e"),
  factor_col = factor(c("x", "y", "x", "y", "z")),
  bool_col = c(TRUE, FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

cat("   Column types:\n")
for (col in colnames(test_df)) {
  cat("     ", col, ":", class(test_df[[col]]), "\n")
}

# Verify factor handling
f <- test_df$factor_col
cat("   Factor levels:", paste(levels(f), collapse = ", "), "\n")
cat("   Factor codes (1-indexed):", paste(as.integer(f), collapse = ", "), "\n")

# Verify NA handling
na_int <- c(1L, NA_integer_, 3L)
cat("   NA integer representation:", na_int[2], "== NA:", is.na(na_int[2]), "\n")
cat("   ✓ data.frame structure verified\n\n")

# ============================================================================
# Test 6: Column-major to row-major conversion
# ============================================================================
cat("6. Testing column-major to row-major conversion:\n")

# R matrices are column-major
r_matrix <- matrix(1:6, nrow = 2, ncol = 3)
cat("   R matrix (column-major):\n")
print(r_matrix)

# Internal storage is column-major
cat("   Internal storage:", paste(as.vector(r_matrix), collapse = ", "), "\n")

# Convert to row-major (as done in Rust)
row_major <- numeric(6)
for (row in 1:2) {
  for (col in 1:3) {
    row_major[(row - 1) * 3 + col] <- r_matrix[row, col]
  }
}
cat("   Row-major conversion:", paste(row_major, collapse = ", "), "\n")
stopifnot(all(row_major == c(1, 3, 5, 2, 4, 6)))
cat("   ✓ Column-major to row-major conversion verified\n\n")

# ============================================================================
# Test 7: Matrix::sparseMatrix construction (used in create_dgcmatrix)
# ============================================================================
cat("7. Testing sparseMatrix construction:\n")

# Method 1: Using new("dgCMatrix", ...) - what Rust actually does
r_i <- c(0L, 1L, 0L)  # 0-indexed row indices
r_p <- c(0L, 1L, 2L, 3L)  # column pointers
r_x <- c(1.0, 3.0, 2.0)  # data
r_dim <- c(2L, 3L)

reconstructed <- new("dgCMatrix",
  i = r_i,
  p = r_p,
  x = r_x,
  Dim = r_dim
)

cat("   Reconstructed matrix (new dgCMatrix):\n")
print(as.matrix(reconstructed))
stopifnot(inherits(reconstructed, "dgCMatrix"))
stopifnot(all(dim(reconstructed) == c(2, 3)))

# Method 2: Using sparseMatrix with i, j, x (triplet format)
# This is an alternative approach
r_i2 <- c(1L, 2L, 1L)  # 1-indexed row indices
r_j2 <- c(1L, 2L, 3L)  # 1-indexed column indices
r_x2 <- c(1.0, 3.0, 2.0)

reconstructed2 <- Matrix::sparseMatrix(
  i = r_i2,
  j = r_j2,
  x = r_x2,
  dims = r_dim
)

cat("   Reconstructed matrix (sparseMatrix triplet):\n")
print(as.matrix(reconstructed2))
stopifnot(inherits(reconstructed2, "dgCMatrix"))
cat("   ✓ sparseMatrix construction verified\n\n")

cat("=== All core logic tests passed ===\n")
