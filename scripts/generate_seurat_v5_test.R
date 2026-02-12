#!/usr/bin/env Rscript
# 生成 Seurat V5 Assay5 测试数据（SimplifiedSeurat 格式）
# 用于测试 CrossCell 的 Seurat V5 支持

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cat("📊 生成 Seurat V5 Assay5 测试数据（SimplifiedSeurat 格式）\n")
cat("=" , rep("=", 50), "\n", sep="")

# 检查 Seurat 版本
seurat_version <- packageVersion("Seurat")
cat(sprintf("Seurat 版本: %s\n", seurat_version))

if (seurat_version < "5.0.0") {
  stop("需要 Seurat >= 5.0.0 来创建 Assay5 对象")
}

# 设置输出目录
output_dir <- "tests/data/real_datasets"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# 辅助函数：将 Seurat V5 对象转换为 SimplifiedSeurat 格式
# ============================================================================
simplify_seurat_v5 <- function(seurat_obj) {
  # 获取 active assay
  active_assay <- DefaultAssay(seurat_obj)
  
  # 提取 assays
  assays_list <- list()
  for (assay_name in Assays(seurat_obj)) {
    assay <- seurat_obj[[assay_name]]
    assay_class <- class(assay)[1]
    
    if (assay_class == "Assay5") {
      # Assay5 格式
      layers_list <- list()
      for (layer_name in Layers(assay)) {
        layer_data <- LayerData(assay, layer = layer_name)
        if (!is.null(layer_data) && length(layer_data) > 0) {
          # 确保是 dgCMatrix
          if (!inherits(layer_data, "dgCMatrix")) {
            layer_data <- as(layer_data, "dgCMatrix")
          }
          layers_list[[layer_name]] <- layer_data
        }
      }
      
      assays_list[[assay_name]] <- list(
        class = "Assay5",
        layers = layers_list,
        features = rownames(assay),
        cells = colnames(assay)
      )
    } else {
      # 传统 Assay 格式
      counts_data <- NULL
      data_data <- NULL
      
      tryCatch({
        counts_data <- GetAssayData(assay, slot = "counts")
        if (!inherits(counts_data, "dgCMatrix")) {
          counts_data <- as(counts_data, "dgCMatrix")
        }
      }, error = function(e) {})
      
      tryCatch({
        data_data <- GetAssayData(assay, slot = "data")
        if (!inherits(data_data, "dgCMatrix")) {
          data_data <- as(data_data, "dgCMatrix")
        }
      }, error = function(e) {})
      
      assays_list[[assay_name]] <- list(
        class = "Assay",
        counts = counts_data,
        data = data_data,
        features = rownames(assay),
        cells = colnames(assay)
      )
    }
  }
  
  # 提取元数据
  meta_data <- seurat_obj@meta.data
  
  # 提取降维结果
  reductions_list <- list()
  for (red_name in names(seurat_obj@reductions)) {
    red <- seurat_obj@reductions[[red_name]]
    reductions_list[[red_name]] <- list(
      cell_embeddings = Embeddings(red),
      feature_loadings = if (length(Loadings(red)) > 0) Loadings(red) else NULL,
      key = Key(red)
    )
  }
  
  # 构建 SimplifiedSeurat 结构
  simplified <- list(
    class = "SimplifiedSeurat",
    project_name = Project(seurat_obj),
    active_assay = active_assay,
    assays = assays_list,
    meta_data = meta_data,
    reductions = reductions_list
  )
  
  return(simplified)
}

# ============================================================================
# 1. 创建小型 Seurat V5 对象（用于快速测试）
# ============================================================================
cat("\n[1/3] 创建小型 Seurat V5 对象...\n")

set.seed(42)
n_cells <- 100
n_genes <- 50

# 创建稀疏计数矩阵
counts <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1)
counts <- abs(counts) * 10  # 确保非负
counts <- as(counts, "dgCMatrix")
rownames(counts) <- paste0("Gene-", 1:n_genes)
colnames(counts) <- paste0("Cell-", 1:n_cells)

# 创建 Seurat V5 对象
seurat_v5_small <- CreateSeuratObject(
  counts = counts,
  project = "SeuratV5_Small",
  assay = "RNA"
)

# 添加元数据
seurat_v5_small$cell_type <- sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)
seurat_v5_small$batch <- sample(c("Batch1", "Batch2"), n_cells, replace = TRUE)
seurat_v5_small$n_counts <- colSums(counts)

# 标准化数据（创建 data layer）
seurat_v5_small <- NormalizeData(seurat_v5_small, verbose = FALSE)

# 检查 Assay 类型
assay_class <- class(seurat_v5_small[["RNA"]])[1]
cat(sprintf("  Assay 类型: %s\n", assay_class))

if (assay_class == "Assay5") {
  layers <- Layers(seurat_v5_small[["RNA"]])
  cat(sprintf("  Layers: %s\n", paste(layers, collapse = ", ")))
}

# 转换为 SimplifiedSeurat 格式
simplified_small <- simplify_seurat_v5(seurat_v5_small)

# 保存
output_path_small <- file.path(output_dir, "seurat_v5_small.rds")
saveRDS(simplified_small, output_path_small)
cat(sprintf("  ✅ 保存到: %s\n", output_path_small))
cat(sprintf("  文件大小: %.2f KB\n", file.size(output_path_small) / 1024))

# ============================================================================
# 2. 创建中型 Seurat V5 对象（包含降维）
# ============================================================================
cat("\n[2/3] 创建中型 Seurat V5 对象（包含降维）...\n")

n_cells_med <- 500
n_genes_med <- 200

# 创建稀疏计数矩阵
counts_med <- Matrix::rsparsematrix(n_genes_med, n_cells_med, density = 0.15)
counts_med <- abs(counts_med) * 10
counts_med <- as(counts_med, "dgCMatrix")
rownames(counts_med) <- paste0("Gene-", 1:n_genes_med)
colnames(counts_med) <- paste0("Cell-", 1:n_cells_med)

# 创建 Seurat V5 对象
seurat_v5_medium <- CreateSeuratObject(
  counts = counts_med,
  project = "SeuratV5_Medium",
  assay = "RNA"
)

# 添加元数据
seurat_v5_medium$cell_type <- sample(c("TypeA", "TypeB", "TypeC", "TypeD"), n_cells_med, replace = TRUE)
seurat_v5_medium$batch <- sample(c("Batch1", "Batch2", "Batch3"), n_cells_med, replace = TRUE)
seurat_v5_medium$sample <- sample(c("Sample1", "Sample2"), n_cells_med, replace = TRUE)

# 标准化和缩放
seurat_v5_medium <- NormalizeData(seurat_v5_medium, verbose = FALSE)
seurat_v5_medium <- FindVariableFeatures(seurat_v5_medium, verbose = FALSE)
seurat_v5_medium <- ScaleData(seurat_v5_medium, verbose = FALSE)

# PCA
seurat_v5_medium <- RunPCA(seurat_v5_medium, npcs = 20, verbose = FALSE)

# UMAP
seurat_v5_medium <- RunUMAP(seurat_v5_medium, dims = 1:10, verbose = FALSE)

# 检查 Assay 类型和 layers
assay_class_med <- class(seurat_v5_medium[["RNA"]])[1]
cat(sprintf("  Assay 类型: %s\n", assay_class_med))

if (assay_class_med == "Assay5") {
  layers_med <- Layers(seurat_v5_medium[["RNA"]])
  cat(sprintf("  Layers: %s\n", paste(layers_med, collapse = ", ")))
}

# 检查降维
reductions <- names(seurat_v5_medium@reductions)
cat(sprintf("  降维: %s\n", paste(reductions, collapse = ", ")))

# 转换为 SimplifiedSeurat 格式
simplified_medium <- simplify_seurat_v5(seurat_v5_medium)

# 保存
output_path_medium <- file.path(output_dir, "seurat_v5_medium.rds")
saveRDS(simplified_medium, output_path_medium)
cat(sprintf("  ✅ 保存到: %s\n", output_path_medium))
cat(sprintf("  文件大小: %.2f KB\n", file.size(output_path_medium) / 1024))

# ============================================================================
# 3. 创建多 Assay Seurat V5 对象
# ============================================================================
cat("\n[3/3] 创建多 Assay Seurat V5 对象...\n")

n_cells_multi <- 200
n_genes_rna <- 100
n_genes_adt <- 20

# RNA counts
counts_rna <- Matrix::rsparsematrix(n_genes_rna, n_cells_multi, density = 0.12)
counts_rna <- abs(counts_rna) * 10
counts_rna <- as(counts_rna, "dgCMatrix")
rownames(counts_rna) <- paste0("Gene-", 1:n_genes_rna)
colnames(counts_rna) <- paste0("Cell-", 1:n_cells_multi)

# ADT counts (蛋白质)
counts_adt <- Matrix::rsparsematrix(n_genes_adt, n_cells_multi, density = 0.5)
counts_adt <- abs(counts_adt) * 100
counts_adt <- as(counts_adt, "dgCMatrix")
rownames(counts_adt) <- paste0("ADT-", 1:n_genes_adt)
colnames(counts_adt) <- paste0("Cell-", 1:n_cells_multi)

# 创建 Seurat V5 对象
seurat_v5_multi <- CreateSeuratObject(
  counts = counts_rna,
  project = "SeuratV5_MultiAssay",
  assay = "RNA"
)

# 添加 ADT assay
seurat_v5_multi[["ADT"]] <- CreateAssay5Object(counts = counts_adt)

# 添加元数据
seurat_v5_multi$cell_type <- sample(c("T_cell", "B_cell", "Monocyte"), n_cells_multi, replace = TRUE)
seurat_v5_multi$donor <- sample(c("Donor1", "Donor2"), n_cells_multi, replace = TRUE)

# 标准化 RNA
seurat_v5_multi <- NormalizeData(seurat_v5_multi, assay = "RNA", verbose = FALSE)

# 标准化 ADT (CLR)
seurat_v5_multi <- NormalizeData(seurat_v5_multi, assay = "ADT", normalization.method = "CLR", verbose = FALSE)

# 检查 Assay 类型
for (assay_name in Assays(seurat_v5_multi)) {
  assay_class_multi <- class(seurat_v5_multi[[assay_name]])[1]
  cat(sprintf("  %s Assay 类型: %s\n", assay_name, assay_class_multi))
  if (assay_class_multi == "Assay5") {
    layers_multi <- Layers(seurat_v5_multi[[assay_name]])
    cat(sprintf("    Layers: %s\n", paste(layers_multi, collapse = ", ")))
  }
}

# 转换为 SimplifiedSeurat 格式
simplified_multi <- simplify_seurat_v5(seurat_v5_multi)

# 保存
output_path_multi <- file.path(output_dir, "seurat_v5_multi_assay.rds")
saveRDS(simplified_multi, output_path_multi)
cat(sprintf("  ✅ 保存到: %s\n", output_path_multi))
cat(sprintf("  文件大小: %.2f KB\n", file.size(output_path_multi) / 1024))

# ============================================================================
# 总结
# ============================================================================
cat("\n")
cat("=" , rep("=", 50), "\n", sep="")
cat("✅ Seurat V5 测试数据生成完成！\n\n")

cat("生成的文件（SimplifiedSeurat 格式）:\n")
cat(sprintf("  1. %s (小型, %d cells × %d genes)\n", 
            output_path_small, n_cells, n_genes))
cat(sprintf("  2. %s (中型, %d cells × %d genes, 含 PCA/UMAP)\n", 
            output_path_medium, n_cells_med, n_genes_med))
cat(sprintf("  3. %s (多 Assay, %d cells, RNA + ADT)\n", 
            output_path_multi, n_cells_multi))

cat("\n测试命令:\n")
cat("  docker-compose run --rm dev cargo test seurat_v5 -- --nocapture\n")
