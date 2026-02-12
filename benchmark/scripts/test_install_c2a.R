#!/usr/bin/env Rscript
# Test install convert2anndata from local source
cat("R library paths:\n")
print(.libPaths())

cat("\nInstalling optparse...\n")
install.packages("optparse", repos="https://cloud.r-project.org")

cat("\nInstalling convert2anndata from /tmp/c2a/convert2anndata-main...\n")
install.packages("/tmp/c2a/convert2anndata-main", repos=NULL, type="source")

cat("\nVerifying...\n")
library(convert2anndata)
cat("convert2anndata version:", as.character(packageVersion("convert2anndata")), "\n")
