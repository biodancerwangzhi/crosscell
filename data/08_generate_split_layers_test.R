#!/usr/bin/env Rscript
# 生成包含 Split Layers 的 Seurat V5 测试数据
# 用于测试 CrossCell 的 Split Layer 自动合并功能
#
# Split Layers 是 Seurat V5 中按批次拆分的层结构：
#   counts.1, counts.2, ... (每个子层包含不同细胞子集，相同基因集)
#
# 使用方法：
#   docker-compose run --rm dev Rscript data/08_generate_split_layers_test.R

library(Seurat)
library(SeuratObject)

cat("=== 生成 Split Layers 测试数据 ===\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version:", as.character(packageVersion("SeuratObject")), "\n\n")

# ============================================
# 方法 1：使用 ifnb 数据集（如果 SeuratData 可用）
# ============================================
generate_from_seuratdata <- function() {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    cat("SeuratData not available, using synthetic approach\n")
    return(NULL)
  }
  
  tryCatch({
    library(SeuratData)
    data("ifnb")
    
    # 确保是 V5 Assay
    ifnb[["RNA"]] <- as(ifnb[["RNA"]], "Assay5")
    
    # 按 stim 列拆分 → 产生 counts.CTRL, counts.STIM, data.CTRL, data.STIM
    ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
    
    cat("Split layers from ifnb:\n")
    print(Layers(ifnb))
    
    return(ifnb)
  }, error = function(e) {
    cat("Failed to load ifnb:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ============================================
# 方法 2：手动构造 Split Layers 对象
# ============================================
generate_synthetic_split <- function() {
  cat("Generating synthetic split layers data...\n")
  
  set.seed(42)
  n_genes <- 200
  n_cells_batch1 <- 100
  n_cells_batch2 <- 150
  n_cells_total <- n_cells_batch1 + n_cells_batch2
  gene_names <- paste0("Gene_", seq_len(n_genes))
  
  # 创建两个批次的计数矩阵
  counts1 <- Matrix::rsparsematrix(n_genes, n_cells_batch1, density = 0.1)
  counts1@x <- abs(counts1@x) * 10  # 确保非负
  counts1@x <- round(counts1@x)
  rownames(counts1) <- gene_names
  colnames(counts1) <- paste0("Batch1_Cell_", seq_len(n_cells_batch1))
  
  counts2 <- Matrix::rsparsematrix(n_genes, n_cells_batch2, density = 0.1)
  counts2@x <- abs(counts2@x) * 10
  counts2@x <- round(counts2@x)
  rownames(counts2) <- gene_names
  colnames(counts2) <- paste0("Batch2_Cell_", seq_len(n_cells_batch2))
  
  # 合并为完整矩阵，创建 Seurat 对象
  counts_full <- cbind(counts1, counts2)
  obj <- CreateSeuratObject(counts = counts_full)
  
  # 添加批次信息
  obj$batch <- c(
    rep("batch1", n_cells_batch1),
    rep("batch2", n_cells_batch2)
  )
  
  # 添加一些元数据（用于测试类型保真）
  obj$n_counts <- as.integer(colSums(counts_full))
  obj$cell_type <- factor(sample(c("T_cell", "B_cell", "Monocyte"), n_cells_total, replace = TRUE))
  
  # 确保是 V5 Assay
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay5")
  
  # 按 batch 拆分 → 产生 counts.batch1, counts.batch2
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  
  cat("Generated object:\n")
  cat("  Cells:", ncol(obj), "\n")
  cat("  Genes:", nrow(obj), "\n")
  cat("  Layers:", paste(Layers(obj), collapse = ", "), "\n")
  cat("  Batch 1 cells:", n_cells_batch1, "\n")
  cat("  Batch 2 cells:", n_cells_batch2, "\n")
  
  return(obj)
}

# ============================================
# 主逻辑
# ============================================

# 先尝试 SeuratData，失败则用合成数据
obj <- generate_from_seuratdata()
if (is.null(obj)) {
  obj <- generate_synthetic_split()
}

# 验证 split layers 存在
layers <- Layers(obj)
cat("\nFinal layers:\n")
print(layers)

# 检查是否有 split 模式的层名
split_pattern <- grep("\\.", layers, value = TRUE)
if (length(split_pattern) == 0) {
  stop("ERROR: No split layers detected! Something went wrong.")
}
cat("\nSplit layers detected:", paste(split_pattern, collapse = ", "), "\n")

# 保存
output_path <- "tests/data/seurat_v5_split_layers.rds"
saveRDS(obj, output_path)
cat("\n✅ Saved to:", output_path, "\n")
cat("File size:", file.size(output_path), "bytes\n")

# 验证可以重新读取
obj2 <- readRDS(output_path)
cat("\nVerification - re-read layers:", paste(Layers(obj2), collapse = ", "), "\n")
cat("Verification - cells:", ncol(obj2), "\n")
cat("Verification - genes:", nrow(obj2), "\n")

cat("\n=== 完成 ===\n")
