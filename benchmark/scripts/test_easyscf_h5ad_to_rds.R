#!/usr/bin/env Rscript
# easySCF H5AD→RDS 测试 — 步骤2: .h5 (自定义格式) → Seurat → RDS
#
# 前置: 先运行 test_easyscf_h5ad_to_rds.py 生成 .h5 文件

suppressPackageStartupMessages({
  library(easySCFr)
  library(Seurat)
  library(jsonlite)
})

TMP_DIR <- "/benchmark/results/easyscf_h2r_tmp"
RESULTS_FILE <- "/benchmark/results/easyscf_h5ad_to_rds_test.json"

# 读取步骤1结果
step1_file <- file.path(TMP_DIR, "step1_results.json")
if (!file.exists(step1_file)) {
  stop("Step 1 results not found. Run test_easyscf_h5ad_to_rds.py first.")
}
step1 <- fromJSON(step1_file)

results <- list()

for (i in seq_len(nrow(step1))) {
  test_id <- step1$test_id[i]
  step1_ok <- step1$step1_ok[i]
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("Testing: %s\n", test_id))
  
  if (!step1_ok) {
    cat("  ⏭️  Step 1 failed, skipping\n")
    results[[length(results) + 1]] <- list(
      test_id = test_id,
      step1_ok = FALSE,
      step2_ok = FALSE,
      step2_error = "step1 failed"
    )
    next
  }
  
  h5_path <- file.path(TMP_DIR, paste0(test_id, "_easyscf.h5"))
  rds_path <- file.path(TMP_DIR, paste0(test_id, "_easyscf.rds"))
  
  if (!file.exists(h5_path)) {
    cat(sprintf("  ❌ .h5 file not found: %s\n", h5_path))
    results[[length(results) + 1]] <- list(
      test_id = test_id,
      step1_ok = TRUE,
      step2_ok = FALSE,
      step2_error = ".h5 file not found"
    )
    next
  }
  
  # 尝试 loadH5 → saveRDS
  tryCatch({
    t0 <- proc.time()["elapsed"]
    obj <- loadH5(h5_path)
    load_time <- proc.time()["elapsed"] - t0
    
    cat(sprintf("  loadH5: ✅ (%s, %.1fs)\n", class(obj)[1], load_time))
    cat(sprintf("  Cells: %d, Genes: %d\n", ncol(obj), nrow(obj)))
    cat(sprintf("  Assays: %s\n", paste(names(obj@assays), collapse = ", ")))
    
    if (length(obj@reductions) > 0) {
      cat(sprintf("  Reductions: %s\n", paste(names(obj@reductions), collapse = ", ")))
    }
    
    t1 <- proc.time()["elapsed"]
    saveRDS(obj, rds_path)
    save_time <- proc.time()["elapsed"] - t1
    
    rds_size <- file.info(rds_path)$size / 1024 / 1024
    cat(sprintf("  saveRDS: ✅ (%.1fs, %.1f MB)\n", save_time, rds_size))
    
    results[[length(results) + 1]] <- list(
      test_id = test_id,
      step1_ok = TRUE,
      step2_ok = TRUE,
      load_time = round(load_time, 2),
      save_time = round(save_time, 2),
      n_cells = ncol(obj),
      n_genes = nrow(obj),
      n_reductions = length(obj@reductions),
      assays = paste(names(obj@assays), collapse = ","),
      rds_size_mb = round(rds_size, 1)
    )
  }, error = function(e) {
    cat(sprintf("  ❌ Error: %s\n", e$message))
    results[[length(results) + 1]] <<- list(
      test_id = test_id,
      step1_ok = TRUE,
      step2_ok = FALSE,
      step2_error = e$message
    )
  })
}

# 汇总
cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
cat("Step 2 Summary (.h5 → Seurat → RDS):\n")
ok_count <- sum(sapply(results, function(r) isTRUE(r$step2_ok)))
total <- length(results)
cat(sprintf("  %d/%d succeeded\n", ok_count, total))
for (r in results) {
  status <- if (isTRUE(r$step2_ok)) "✅" else "❌"
  err <- if (!isTRUE(r$step2_ok)) paste0(" — ", r$step2_error) else ""
  cat(sprintf("  %s %s%s\n", status, r$test_id, err))
}

# 保存最终结果
write(toJSON(results, auto_unbox = TRUE, pretty = TRUE), RESULTS_FILE)
cat(sprintf("\nResults saved to %s\n", RESULTS_FILE))
