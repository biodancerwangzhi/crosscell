# Seurat 数据生成脚本 - 自动检测版本生成对应格式
library(Seurat)

seurat_ver <- as.character(packageVersion("Seurat"))
major_ver <- as.integer(strsplit(seurat_ver, "\\.")[[1]][1])
prefix <- paste0("seurat_v", major_ver)
is_v5 <- major_ver >= 5

cat("=== Seurat", seurat_ver, "(V", major_ver, ") ===\n")
if (is_v5) options(Seurat.object.assay.version = "v5")

output_dir <- "/workspace/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_obj <- function(obj, name) {
  f <- file.path(output_dir, paste0(prefix, "_", name, ".rds"))
  saveRDS(obj, f)
  cat("  Saved:", basename(f), "\n")
}

process <- function(obj) {
  tryCatch({
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
    obj <- FindClusters(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
    obj
  }, error = function(e) NULL)
}

convert_to_v5 <- function(obj) {
  tryCatch({
    obj <- UpdateSeuratObject(obj)
    
    # 获取默认 assay 名称
    default_assay <- DefaultAssay(obj)
    counts <- GetAssayData(obj, assay = default_assay, layer = "counts")
    if (is.null(counts) || ncol(counts) == 0) {
      counts <- GetAssayData(obj, assay = default_assay, layer = "data")
    }
    
    # 创建新对象
    new_obj <- CreateSeuratObject(
      counts = counts,
      meta.data = obj@meta.data,
      project = obj@project.name,
      assay = default_assay
    )
    
    # 保留空间图像数据
    if (length(Images(obj)) > 0) {
      for (img_name in Images(obj)) {
        tryCatch({
          new_obj[[img_name]] <- obj[[img_name]]
          cat("  Preserved spatial image:", img_name, "\n")
        }, error = function(e) {
          cat("  Warning: Could not preserve image", img_name, "\n")
        })
      }
    }
    
    return(new_obj)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
}

load_dataset <- function(pkg_name) {
  data_name <- sub("\\.SeuratData$", "", pkg_name)
  tryCatch({
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      e <- new.env()
      data(list = data_name, package = pkg_name, envir = e)
      if (exists(data_name, envir = e)) return(list(name = data_name, obj = get(data_name, envir = e)))
    }
    NULL
  }, error = function(e) NULL)
}

# 处理所有已安装的数据包
installed_pkgs <- installed.packages()[, "Package"]
seurat_data_pkgs <- grep("\\.SeuratData$", installed_pkgs, value = TRUE)
cat("Found", length(seurat_data_pkgs), "data packages\n\n")

for (pkg in seurat_data_pkgs) {
  data_name <- sub("\\.SeuratData$", "", pkg)
  cat("Processing", data_name, "...\n")
  
  result <- load_dataset(pkg)
  if (is.null(result) || !inherits(result$obj, "Seurat")) {
    cat("  Skipped (not Seurat object)\n")
    next
  }
  
  obj <- result$obj
  if (is_v5) {
    obj <- convert_to_v5(obj)
    if (is.null(obj)) { cat("  Skipped (conversion failed)\n"); next }
  }
  
  save_obj(obj, paste0(data_name, "_raw"))
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  processed <- process(obj)
  if (!is.null(processed)) save_obj(processed, paste0(data_name, "_processed"))
}

cat("\n=== Done ===\n")
