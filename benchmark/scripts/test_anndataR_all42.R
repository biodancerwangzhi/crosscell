#!/usr/bin/env Rscript
# Test anndataR as_AnnData + write_h5ad on ALL 42 RDS datasets

suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
  library(jsonlite)
})

cat("anndataR version:", as.character(packageVersion("anndataR")), "\n")

config <- fromJSON("/benchmark/config/benchmark_testcases.json")
cases <- config$test_cases$rds_to_h5ad

results <- list()
ok_count <- 0
fail_count <- 0

for (i in seq_len(nrow(cases))) {
  test_id <- cases$test_id[i]
  fname <- cases$file[i]
  path <- file.path("/benchmark/data/generated", fname)
  out <- file.path("/tmp", sub(".rds", "_anndataR.h5ad", fname))
  
  cat(sprintf("\n[%d/%d] %s\n", i, nrow(cases), test_id))
  
  tryCatch({
    obj <- readRDS(path)
    tryCatch({obj <- UpdateSeuratObject(obj)}, error=function(e){})
    
    t0 <- proc.time()["elapsed"]
    ad <- as_AnnData(obj)
    write_h5ad(ad, out)
    elapsed <- proc.time()["elapsed"] - t0
    
    cat(sprintf("  ✅ %.1fs\n", elapsed))
    ok_count <- ok_count + 1
    results[[length(results) + 1]] <- list(test_id=test_id, ok=TRUE, time=round(elapsed,2))
    
    # Clean up
    rm(obj, ad); gc(verbose=FALSE)
    file.remove(out)
    
  }, error=function(e) {
    cat(sprintf("  ❌ %s\n", e$message))
    fail_count <<- fail_count + 1
    results[[length(results) + 1]] <<- list(test_id=test_id, ok=FALSE, error=e$message)
    gc(verbose=FALSE)
  })
}

cat(sprintf("\n\n=== Summary: %d/%d passed, %d failed ===\n", ok_count, nrow(cases), fail_count))
for (r in results) {
  if (!isTRUE(r$ok)) {
    cat(sprintf("  ❌ %s: %s\n", r$test_id, r$error))
  }
}
