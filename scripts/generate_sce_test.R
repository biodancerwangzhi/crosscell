#!/usr/bin/env Rscript
# Generate test SingleCellExperiment objects for CrossCell testing

# Check if SingleCellExperiment is installed
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  message("Installing SingleCellExperiment...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("SingleCellExperiment")
}

library(SingleCellExperiment)
library(Matrix)

# Create output directory
dir.create("tests/data", showWarnings = FALSE, recursive = TRUE)

# ============================================
# 1. Minimal SCE object (counts only)
# ============================================
message("Creating minimal SCE object...")

set.seed(42)
n_cells <- 100
n_genes <- 50

# Create sparse count matrix
counts <- Matrix(
  rpois(n_cells * n_genes, lambda = 5),
  nrow = n_genes,
  ncol = n_cells,
  sparse = TRUE
)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Cell", 1:n_cells)

# Create minimal SCE
sce_minimal <- SingleCellExperiment(
  assays = list(counts = counts)
)

saveRDS(sce_minimal, "tests/data/sce_minimal.rds")
message("  Saved: tests/data/sce_minimal.rds")

# ============================================
# 2. SCE with metadata
# ============================================
message("Creating SCE with metadata...")

# Add cell metadata (colData)
col_data <- DataFrame(
  cell_type = factor(sample(c("T cell", "B cell", "Monocyte"), n_cells, replace = TRUE)),
  batch = factor(sample(c("batch1", "batch2"), n_cells, replace = TRUE)),
  n_genes = colSums(counts > 0),
  total_counts = colSums(counts)
)

# Add gene metadata (rowData)
row_data <- DataFrame(
  gene_name = rownames(counts),
  gene_type = factor(sample(c("protein_coding", "lncRNA"), n_genes, replace = TRUE)),
  mean_expression = rowMeans(as.matrix(counts))
)

sce_with_metadata <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = col_data,
  rowData = row_data
)

saveRDS(sce_with_metadata, "tests/data/sce_with_metadata.rds")
message("  Saved: tests/data/sce_with_metadata.rds")

# ============================================
# 3. SCE with reductions
# ============================================
message("Creating SCE with reductions...")

# Create logcounts
logcounts <- log1p(counts)

# Create PCA (simulated)
pca_coords <- matrix(rnorm(n_cells * 50), nrow = n_cells, ncol = 50)
rownames(pca_coords) <- colnames(counts)
colnames(pca_coords) <- paste0("PC", 1:50)

# Create UMAP (simulated)
umap_coords <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)
rownames(umap_coords) <- colnames(counts)
colnames(umap_coords) <- c("UMAP1", "UMAP2")

# Create tSNE (simulated)
tsne_coords <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)
rownames(tsne_coords) <- colnames(counts)
colnames(tsne_coords) <- c("tSNE1", "tSNE2")

sce_with_reductions <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = logcounts),
  colData = col_data,
  rowData = row_data,
  reducedDims = list(PCA = pca_coords, UMAP = umap_coords, TSNE = tsne_coords)
)

saveRDS(sce_with_reductions, "tests/data/sce_with_reductions.rds")
message("  Saved: tests/data/sce_with_reductions.rds")

# ============================================
# 4. Full SCE object
# ============================================
message("Creating full SCE object...")

# Add metadata (uns equivalent)
metadata(sce_with_reductions) <- list(
  project_name = "CrossCell Test",
  creation_date = Sys.Date(),
  software_version = "1.0.0"
)

saveRDS(sce_with_reductions, "tests/data/sce_full.rds")
message("  Saved: tests/data/sce_full.rds")

# ============================================
# Print structure for analysis
# ============================================
message("\n=== SCE Structure Analysis ===")
message("\nClass hierarchy:")
print(class(sce_with_reductions))

message("\nSlots:")
print(slotNames(sce_with_reductions))

message("\nAssays:")
print(assayNames(sce_with_reductions))

message("\ncolData columns:")
print(colnames(colData(sce_with_reductions)))

message("\nrowData columns:")
print(colnames(rowData(sce_with_reductions)))

message("\nreducedDims:")
print(reducedDimNames(sce_with_reductions))

message("\nDone! Test SCE files created.")
