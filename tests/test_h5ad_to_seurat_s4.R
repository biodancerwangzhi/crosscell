library(SeuratObject)

# Find an H5AD file to test with
h5ad_files <- list.files("tests/data", pattern="\\.h5ad$", full.names=TRUE)
if (length(h5ad_files) == 0) {
  h5ad_files <- list.files("data/generated", pattern="\\.h5ad$", full.names=TRUE)
}
cat("Available H5AD files:", h5ad_files, "\n")

# Test with the CrossCell-generated RDS
rds_file <- "tests/data/rust_generated_seurat.rds"
if (!file.exists(rds_file)) {
  cat("ERROR: RDS file not found\n")
  quit(status=1)
}

obj <- readRDS(rds_file)
cat("\n=== Seurat Object Summary ===\n")
cat("Class:", class(obj), "\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Active assay:", obj@active.assay, "\n")
cat("Assays:", names(obj@assays), "\n")

# Test basic Seurat operations
cat("\n=== Testing Seurat Operations ===\n")

# 1. Access counts
tryCatch({
  counts <- GetAssayData(obj, layer="counts")
  cat("GetAssayData: OK (", class(counts), dim(counts), ")\n")
}, error = function(e) cat("GetAssayData ERROR:", e$message, "\n"))

# 2. Access meta.data
tryCatch({
  md <- obj@meta.data
  cat("meta.data: OK (", nrow(md), "rows,", ncol(md), "cols)\n")
}, error = function(e) cat("meta.data ERROR:", e$message, "\n"))

# 3. Cell names
tryCatch({
  cn <- Cells(obj)
  cat("Cells(): OK (", length(cn), "cells, first 3:", head(cn, 3), ")\n")
}, error = function(e) cat("Cells() ERROR:", e$message, "\n"))

# 4. Feature names
tryCatch({
  fn <- Features(obj)
  cat("Features(): OK (", length(fn), "features, first 3:", head(fn, 3), ")\n")
}, error = function(e) cat("Features() ERROR:", e$message, "\n"))

# 5. Idents
tryCatch({
  id <- Idents(obj)
  cat("Idents(): OK (", length(id), "idents, levels:", levels(id), ")\n")
}, error = function(e) cat("Idents() ERROR:", e$message, "\n"))

# 6. Subset
tryCatch({
  sub <- obj[, 1:3]
  cat("Subset [,1:3]: OK (", ncol(sub), "cells)\n")
}, error = function(e) cat("Subset ERROR:", e$message, "\n"))

cat("\n=== All tests completed ===\n")
