# 提取 Seurat 数据集信息
library(Seurat)

seurat_ver <- as.character(packageVersion("Seurat"))
major_ver <- as.integer(strsplit(seurat_ver, "\\.")[[1]][1])
version_prefix <- paste0("v", major_ver)

output_dir <- "/workspace/output"
pattern <- paste0("seurat_", version_prefix, "_.*_raw\\.rds$")
files <- list.files(output_dir, pattern = pattern, full.names = TRUE)

cat("Seurat", seurat_ver, "- Found", length(files), "raw RDS files\n\n")

results <- lapply(files, function(f) {
  cat("Inspecting:", basename(f), "\n")
  obj <- readRDS(f)
  
  # 安全获取 Images
  has_spatial <- tryCatch({
    length(Images(obj)) > 0
  }, error = function(e) FALSE)
  
  # 安全获取 Reductions
  reds <- tryCatch({
    r <- Reductions(obj)
    if (length(r) > 0) paste(r, collapse = ",") else ""
  }, error = function(e) "")
  
  list(
    file = basename(f),
    version = paste0("V", major_ver),
    dataset = gsub(paste0("seurat_", version_prefix, "_(.*)_raw\\.rds"), "\\1", basename(f)),
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    assays = paste(Assays(obj), collapse = ","),
    default_assay = DefaultAssay(obj),
    assay_class = class(obj[[DefaultAssay(obj)]])[1],
    has_spatial = has_spatial,
    reductions = reds,
    file_size_mb = round(file.size(f) / 1024^2, 2)
  )
})

df <- do.call(rbind, lapply(results, as.data.frame))
csv_file <- file.path(output_dir, paste0("seurat_", version_prefix, "_summary.csv"))
write.csv(df, csv_file, row.names = FALSE)

cat("\n=== Summary ===\n")
print(df[, c("dataset", "n_cells", "n_genes", "assay_class", "file_size_mb")])
cat("\nSaved to:", csv_file, "\n")
