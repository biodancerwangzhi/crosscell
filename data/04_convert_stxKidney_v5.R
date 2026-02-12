# 专门转换 stxKidney 空间数据到 V5 格式（保留空间坐标）
library(Seurat)
options(Seurat.object.assay.version = "v5")

output_dir <- "/workspace/output"

# 读取 V4 版本
v4_file <- file.path(output_dir, "seurat_v4_stxKidney_raw.rds")
cat("Reading V4 file:", v4_file, "\n")
obj <- readRDS(v4_file)

cat("Original object:\n")
cat("  Cells:", ncol(obj), "\n")
cat("  Genes:", nrow(obj), "\n")
cat("  Assay class:", class(obj[["Spatial"]])[1], "\n")
cat("  Images:", paste(Images(obj), collapse = ", "), "\n")

# 提取数据
counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
meta <- obj@meta.data

# 尝试提取空间坐标
cat("\nExtracting spatial coordinates...\n")
img <- obj@images$image
coords <- NULL

tryCatch({
  # SliceImage 的坐标存储方式
  if ("coordinates" %in% slotNames(img)) {
    coords <- img@coordinates
    cat("  Found coordinates:", nrow(coords), "spots\n")
  }
}, error = function(e) {
  cat("  Could not extract coordinates:", e$message, "\n")
})

# 创建新的 V5 对象
cat("\nCreating V5 object...\n")
new_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  assay = "Spatial"
)

# 如果有坐标，添加到 meta.data
if (!is.null(coords) && nrow(coords) > 0) {
  # 确保行名匹配
  common_cells <- intersect(rownames(coords), colnames(new_obj))
  if (length(common_cells) > 0) {
    new_obj$spatial_x <- coords[colnames(new_obj), "imagerow"]
    new_obj$spatial_y <- coords[colnames(new_obj), "imagecol"]
    cat("  Added spatial coordinates to meta.data\n")
  }
}

cat("New object:\n")
cat("  Assay class:", class(new_obj[["Spatial"]])[1], "\n")

# 保存
v5_file <- file.path(output_dir, "seurat_v5_stxKidney_raw.rds")
saveRDS(new_obj, v5_file)
cat("\nSaved:", v5_file, "\n")
cat("File size:", round(file.size(v5_file) / 1024^2, 2), "MB\n")

if (!is.null(coords)) {
  cat("\nNote: Spatial coordinates preserved in meta.data (spatial_x, spatial_y)\n")
  cat("      Original SliceImage not compatible with V5\n")
} else {
  cat("\nNote: Spatial images not preserved (SliceImage incompatible with V5)\n")
}
