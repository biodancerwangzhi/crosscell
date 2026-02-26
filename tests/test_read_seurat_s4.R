library(SeuratObject)
cat("Reading CrossCell-generated RDS...\n")
obj <- tryCatch(
  readRDS("tests/data/rust_generated_seurat.rds"),
  error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
)
if (!is.null(obj)) {
  cat("class:", class(obj), "\n")
  cat("is S4:", isS4(obj), "\n")
  if (isS4(obj)) {
    cat("slots:", slotNames(obj), "\n")
    cat("active.assay:", obj@active.assay, "\n")
    cat("project.name:", obj@project.name, "\n")
    cat("version:", as.character(obj@version), "\n")
    
    # Try ncol/nrow
    nc <- tryCatch(ncol(obj), error = function(e) paste("ERROR:", e$message))
    nr <- tryCatch(nrow(obj), error = function(e) paste("ERROR:", e$message))
    cat("ncol:", nc, "\n")
    cat("nrow:", nr, "\n")
    
    cat("assay names:", names(obj@assays), "\n")
    a <- obj@assays[["RNA"]]
    cat("assay class:", class(a), "\n")
    cat("assay is S4:", isS4(a), "\n")
    if (isS4(a)) {
      cat("assay key:", a@key, "\n")
      cat("layer names:", names(a@layers), "\n")
      m <- a@layers[["counts"]]
      cat("counts class:", class(m), "\n")
      cat("counts dim:", dim(m), "\n")
      cat("counts rownames (first 5):", head(rownames(m), 5), "\n")
      cat("counts colnames (first 5):", head(colnames(m), 5), "\n")
      cat("cells class:", class(a@cells), "\n")
      cat("features class:", class(a@features), "\n")
    }
    
    # Test active.ident
    cat("active.ident class:", class(obj@active.ident), "\n")
    cat("active.ident levels:", levels(obj@active.ident), "\n")
    cat("active.ident names (first 5):", head(names(obj@active.ident), 5), "\n")
    
    # Test meta.data
    cat("meta.data class:", class(obj@meta.data), "\n")
    cat("meta.data dim:", dim(obj@meta.data), "\n")
    cat("meta.data rownames (first 5):", head(rownames(obj@meta.data), 5), "\n")
    
    cat("\n=== SUCCESS: Seurat S4 object loaded correctly! ===\n")
  }
} else {
  cat("\n=== FAILED: Could not read RDS file ===\n")
}
