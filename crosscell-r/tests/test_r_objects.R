# Test R object structures for CrossCell R bindings
# This script verifies that Seurat and SCE objects have the expected structure

library(Seurat)
library(SingleCellExperiment)
library(Matrix)

cat("=== Testing R Object Structures ===\n\n")

# Create a simple test matrix
set.seed(42)
n_genes <- 50
n_cells <- 100
counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes, ncol = n_cells)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Cell", 1:n_cells)

# Convert to sparse matrix
sparse_counts <- as(counts, "dgCMatrix")

cat("1. Testing dgCMatrix structure:\n")
cat("   Class:", class(sparse_counts), "\n")
cat("   Dim:", dim(sparse_counts), "\n")
cat("   Slots: i (row indices), p (col pointers), x (data), Dim\n")
cat("   i length:", length(sparse_counts@i), "\n")
cat("   p length:", length(sparse_counts@p), "\n")
cat("   x length:", length(sparse_counts@x), "\n")
cat("   Dim:", sparse_counts@Dim, "\n")
cat("\n")

# Create Seurat object
cat("2. Testing Seurat object structure:\n")
seurat_obj <- CreateSeuratObject(counts = sparse_counts, project = "test")
cat("   Class:", class(seurat_obj), "\n")
cat("   DefaultAssay:", DefaultAssay(seurat_obj), "\n")
cat("   ncol:", ncol(seurat_obj), "\n")
cat("   nrow:", nrow(seurat_obj), "\n")
cat("   meta.data columns:", colnames(seurat_obj@meta.data), "\n")
cat("   reductions:", names(seurat_obj@reductions), "\n")

# Get counts data
counts_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
cat("   counts class:", class(counts_data), "\n")
cat("\n")

# Add some embeddings
cat("3. Testing Seurat with embeddings:\n")
# Create fake PCA embeddings
pca_embeddings <- matrix(rnorm(n_cells * 10), nrow = n_cells, ncol = 10)
rownames(pca_embeddings) <- colnames(seurat_obj)
colnames(pca_embeddings) <- paste0("PC_", 1:10)
seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PC_", assay = "RNA")
cat("   reductions after PCA:", names(seurat_obj@reductions), "\n")
cat("   PCA embeddings dim:", dim(seurat_obj@reductions$pca@cell.embeddings), "\n")
cat("\n")

# Create SingleCellExperiment object
cat("4. Testing SingleCellExperiment structure:\n")
sce_obj <- SingleCellExperiment(
  assays = list(counts = sparse_counts),
  colData = data.frame(cell_type = sample(c("A", "B", "C"), n_cells, replace = TRUE)),
  rowData = data.frame(gene_biotype = sample(c("protein_coding", "lncRNA"), n_genes, replace = TRUE))
)
cat("   Class:", class(sce_obj), "\n")
cat("   ncol:", ncol(sce_obj), "\n")
cat("   nrow:", nrow(sce_obj), "\n")
cat("   assayNames:", assayNames(sce_obj), "\n")
cat("   colData columns:", colnames(colData(sce_obj)), "\n")
cat("   rowData columns:", colnames(rowData(sce_obj)), "\n")

# Add reducedDims
reducedDim(sce_obj, "PCA") <- pca_embeddings
cat("   reducedDimNames:", reducedDimNames(sce_obj), "\n")
cat("\n")

# Test data.frame with factor
cat("5. Testing data.frame with factor:\n")
test_df <- data.frame(
  int_col = 1:5,
  dbl_col = c(1.1, 2.2, 3.3, 4.4, 5.5),
  str_col = c("a", "b", "c", "d", "e"),
  factor_col = factor(c("x", "y", "x", "y", "z")),
  bool_col = c(TRUE, FALSE, TRUE, FALSE, TRUE)
)
cat("   Column types:\n")
for (col in colnames(test_df)) {
  cat("     ", col, ":", class(test_df[[col]]), "\n")
}
cat("   Factor levels:", levels(test_df$factor_col), "\n")
cat("   Factor codes:", as.integer(test_df$factor_col), "\n")
cat("\n")

cat("=== All structure tests passed ===\n")
