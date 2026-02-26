#' Write Seurat/SCE object to H5AD format
#'
#' @param obj A Seurat or SingleCellExperiment object.
#' @param path Output file path (.h5ad).
#' @param normalize Logical. Apply library-size normalization + log1p before writing. Default FALSE.
#' @param top_genes Integer or NULL. Select top N variable genes before writing. Default NULL.
#' @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
#' @export
#' @examples
#' \dontrun{
#' write_as_h5ad(seurat_obj, "output.h5ad")
#' write_as_h5ad(seurat_obj, "output.h5ad", normalize = TRUE, top_genes = 2000)
#' }
write_as_h5ad <- function(obj, path, normalize = FALSE, top_genes = NULL,
                          gene_id_column = NULL) {
    if (!inherits(obj, c("Seurat", "SingleCellExperiment"))) {
        stop("obj must be a Seurat or SingleCellExperiment object")
    }
    .Call(wrap__write_as_h5ad, obj, path, normalize, top_genes, gene_id_column)
    invisible(NULL)
}

#' Fast write of Seurat object to RDS format
#'
#' Uses Rust-based RDS writer for faster serialization.
#'
#' @param obj A Seurat object.
#' @param path Output file path (.rds).
#' @export
#' @examples
#' \dontrun{
#' write_rds_fast(seurat_obj, "output.rds")
#' }
write_rds_fast <- function(obj, path) {
    if (!inherits(obj, "Seurat")) {
        stop("obj must be a Seurat object")
    }
    .Call(wrap__write_rds_fast, obj, path)
    invisible(NULL)
}

#' Inspect file metadata
#'
#' @param path Path to an .h5ad or .rds file.
#' @return A list with file metadata including cell/gene counts, format info, etc.
#' @export
#' @examples
#' \dontrun{
#' info <- inspect_file("data.h5ad")
#' print(info$n_cells)
#' }
inspect_file <- function(path) {
    if (!file.exists(path)) {
        stop("File not found: ", path)
    }
    ext <- tolower(tools::file_ext(path))
    if (!ext %in% c("h5ad", "rds")) {
        stop("Unsupported file format: .", ext, ". Supported: .h5ad, .rds")
    }
    .Call(wrap__inspect_file, path)
}

#' Validate conversion fidelity between two files
#'
#' Compares an original file with a converted file and returns a validation
#' report including Pearson R, MSE, RMSE, exact match percentage, and
#' optionally ARI/NMI for cluster labels.
#'
#' @param original Path to the original file (.h5ad or .rds).
#' @param converted Path to the converted file (.h5ad or .rds).
#' @param tolerance Numeric. Numerical tolerance. Default 1e-7.
#' @param cluster_column Character or NULL. Column in obs for cluster labels. Default NULL.
#' @return A list with validation results (passed, summary, detailed_report, and optionally ari, nmi).
#' @export
#' @examples
#' \dontrun{
#' result <- validate_conversion("original.h5ad", "converted.h5ad")
#' if (result$passed) message("Validation passed!")
#' }
validate_conversion <- function(original, converted, tolerance = 1e-7,
                                cluster_column = NULL) {
    if (!file.exists(original)) {
        stop("Original file not found: ", original)
    }
    if (!file.exists(converted)) {
        stop("Converted file not found: ", converted)
    }
    .Call(wrap__validate_conversion, original, converted, tolerance, cluster_column)
}
