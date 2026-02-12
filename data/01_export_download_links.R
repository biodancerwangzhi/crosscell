# 导出 SeuratData 下载链接（用于 IDM 手动下载）
library(SeuratData)

available <- AvailableData()
base_url <- "http://seurat.nygenome.org/src/contrib/"

links <- data.frame(
  dataset = rownames(available),
  version = available$Version,
  url = paste0(base_url, rownames(available), "_", available$Version, ".tar.gz"),
  stringsAsFactors = FALSE
)

output_dir <- "/workspace/packages"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

writeLines(links$url, file.path(output_dir, "seuratdata_download_links.txt"))
write.csv(links, file.path(output_dir, "seuratdata_datasets.csv"), row.names = FALSE)

cat("Exported", nrow(links), "download links to", output_dir, "\n")
