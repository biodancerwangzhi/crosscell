#' Read H5AD file as Seurat object
#'
#' @param path Path to the .h5ad file.
#' @param normalize Logical. Apply library-size normalization + log1p. Default FALSE.
#' @param top_genes Integer or NULL. Select top N variable genes. Default NULL.
#' @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
#' @return A Seurat object.
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- read_h5ad_as_seurat("data.h5ad")
#' seurat_obj <- read_h5ad_as_seurat("data.h5ad", normalize = TRUE, top_genes = 2000)
#' }
read_h5ad_as_seurat <- function(path, normalize = FALSE, top_genes = NULL, gene_id_column = NULL) {
    if (!file.exists(path)) {
        stop("File not found: ", path)
    }
    .Call(wrap__read_h5ad_as_seurat, path, normalize, top_genes, gene_id_column)
}

#' Read H5AD file as SingleCellExperiment object
#'
#' @param path Path to the .h5ad file.
#' @return A SingleCellExperiment object.
#' @export
#' @examples
#' \dontrun{
#' sce_obj <- read_h5ad_as_sce("data.h5ad")
#' }
read_h5ad_as_sce <- function(path) {
    if (!file.exists(path)) {
        stop("File not found: ", path)
    }
    .Call(wrap__read_h5ad_as_sce, path)
}

#' Fast read of Seurat RDS file
#'
#' Uses Rust-based RDS parser for faster reading of Seurat objects.
#'
#' @param path Path to the .rds file.
#' @param normalize Logical. Apply library-size normalization + log1p. Default FALSE.
#' @param top_genes Integer or NULL. Select top N variable genes. Default NULL.
#' @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
#' @param keep_layers Logical. Preserve Seurat V5 split layers. Default FALSE.
#' @return A Seurat object.
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- read_rds_fast("seurat.rds")
#' seurat_obj <- read_rds_fast("seurat.rds", normalize = TRUE, top_genes = 2000)
#' seurat_obj <- read_rds_fast("seurat.rds", keep_layers = TRUE)
#' }
read_rds_fast <- function(path, normalize = FALSE, top_genes = NULL,
                          gene_id_column = NULL, keep_layers = FALSE) {
    if (!file.exists(path)) {
        stop("File not found: ", path)
    }
    .Call(wrap__read_rds_fast, path, normalize, top_genes, gene_id_column, keep_layers)
}
