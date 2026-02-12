cat("=== .libPaths() ===\n")
cat(paste(.libPaths(), collapse="\n"), "\n")

cat("\n=== Seurat location ===\n")
cat(system.file(package="Seurat"), "\n")

cat("\n=== anndata location ===\n")
cat(system.file(package="anndata"), "\n")

cat("\n=== convert2anndata location ===\n")
cat(system.file(package="convert2anndata"), "\n")

cat("\n=== Seurat loadable? ===\n")
cat(requireNamespace("Seurat", quietly=TRUE), "\n")

cat("\n=== anndata loadable? ===\n")
cat(requireNamespace("anndata", quietly=TRUE), "\n")
