//! CrossCell R bindings
//!
//! Provides R functions for reading/writing single-cell data formats.

use extendr_api::prelude::*;
use std::path::Path;

pub mod convert_to_r;
pub mod convert_from_r;

#[cfg(test)]
mod tests;

use convert_to_r::{ir_to_seurat, ir_to_sce};
use convert_from_r::{seurat_to_ir, sce_to_ir};

// ============================================================================
// read_h5ad_as_seurat and read_h5ad_as_sce
// ============================================================================

/// Read an H5AD file and return a Seurat object.
///
/// @param path Path to the .h5ad file.
/// @param normalize Logical. Apply library-size normalization + log1p. Default FALSE.
/// @param top_genes Integer or NULL. Select top N variable genes. Default NULL.
/// @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
/// @return A Seurat object.
/// @export
#[extendr]
fn read_h5ad_as_seurat(
    path: &str,
    normalize: Option<bool>,
    top_genes: Option<i32>,
    gene_id_column: Option<String>,
) -> Result<Robj> {
    if !Path::new(path).exists() {
        return Err(format!("File not found: {}", path).into());
    }
    if !path.to_lowercase().ends_with(".h5ad") {
        return Err(format!("Expected .h5ad file, got: {}", path).into());
    }

    let mut ir = crosscell::anndata::read_h5ad(path)
        .map_err(|e| Error::Other(format!("Failed to read H5AD: {}", e)))?;

    apply_transforms(&mut ir, normalize.unwrap_or(false), top_genes, gene_id_column)?;

    ir_to_seurat(&ir)
}

/// Read an H5AD file and return a SingleCellExperiment object.
///
/// @param path Path to the .h5ad file.
/// @return A SingleCellExperiment object.
/// @export
#[extendr]
fn read_h5ad_as_sce(path: &str) -> Result<Robj> {
    if !Path::new(path).exists() {
        return Err(format!("File not found: {}", path).into());
    }
    if !path.to_lowercase().ends_with(".h5ad") {
        return Err(format!("Expected .h5ad file, got: {}", path).into());
    }

    let ir = crosscell::anndata::read_h5ad(path)
        .map_err(|e| Error::Other(format!("Failed to read H5AD: {}", e)))?;

    ir_to_sce(&ir)
}

// ============================================================================
// read_rds_fast
// ============================================================================

/// Fast read of a Seurat RDS file.
///
/// @param path Path to the .rds file.
/// @param normalize Logical. Apply library-size normalization + log1p. Default FALSE.
/// @param top_genes Integer or NULL. Select top N variable genes. Default NULL.
/// @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
/// @param keep_layers Logical. Preserve Seurat V5 split layers. Default FALSE.
/// @return A Seurat object.
/// @export
#[extendr]
fn read_rds_fast(
    path: &str,
    normalize: Option<bool>,
    top_genes: Option<i32>,
    gene_id_column: Option<String>,
    keep_layers: Option<bool>,
) -> Result<Robj> {
    if !Path::new(path).exists() {
        return Err(format!("File not found: {}", path).into());
    }
    if !path.to_lowercase().ends_with(".rds") {
        return Err(format!("Expected .rds file, got: {}", path).into());
    }

    let kl = keep_layers.unwrap_or(false);
    let result = crosscell::seurat::direct_read::read_seurat_direct(path, kl)
        .map_err(|e| Error::Other(format!("Failed to read RDS: {}", e)))?;

    let mut data = result.data;
    apply_transforms(&mut data, normalize.unwrap_or(false), top_genes, gene_id_column)?;

    ir_to_seurat(&data)
}

// ============================================================================
// write_as_h5ad
// ============================================================================

/// Write a Seurat or SCE object to H5AD format.
///
/// @param obj A Seurat or SingleCellExperiment object.
/// @param path Output file path (.h5ad).
/// @param normalize Logical. Apply library-size normalization + log1p before writing. Default FALSE.
/// @param top_genes Integer or NULL. Select top N variable genes before writing. Default NULL.
/// @param gene_id_column Character or NULL. Use specified var column as gene IDs. Default NULL.
/// @export
#[extendr]
fn write_as_h5ad(
    obj: Robj,
    path: &str,
    normalize: Option<bool>,
    top_genes: Option<i32>,
    gene_id_column: Option<String>,
) -> Result<()> {
    if !path.to_lowercase().ends_with(".h5ad") {
        return Err(format!("Expected .h5ad output path, got: {}", path).into());
    }

    let mut ir = if obj.inherits("Seurat") {
        seurat_to_ir(obj)?
    } else if obj.inherits("SingleCellExperiment") {
        sce_to_ir(obj)?
    } else {
        return Err("Expected Seurat or SingleCellExperiment object".into());
    };

    apply_transforms(&mut ir, normalize.unwrap_or(false), top_genes, gene_id_column)?;

    crosscell::anndata::write_h5ad(&ir, path)
        .map_err(|e| Error::Other(format!("Failed to write H5AD: {}", e)))?;

    Ok(())
}

// ============================================================================
// write_rds_fast
// ============================================================================

/// Fast write of a Seurat object to RDS format.
///
/// @param obj A Seurat object.
/// @param path Output file path (.rds).
/// @export
#[extendr]
fn write_rds_fast(obj: Robj, path: &str) -> Result<()> {
    if !path.to_lowercase().ends_with(".rds") {
        return Err(format!("Expected .rds output path, got: {}", path).into());
    }

    if !obj.inherits("Seurat") {
        return Err("Expected Seurat object".into());
    }

    let ir = seurat_to_ir(obj)?;

    crosscell::seurat::ir_to_seurat::write_seurat_rds(&ir, path)
        .map_err(|e| Error::Other(format!("Failed to write RDS: {}", e)))?;

    Ok(())
}

// ============================================================================
// inspect_file
// ============================================================================

/// Inspect a file and return metadata.
///
/// @param path Path to an .h5ad or .rds file.
/// @return A list with file metadata.
/// @export
#[extendr]
fn inspect_file(path: &str) -> Result<Robj> {
    if !Path::new(path).exists() {
        return Err(format!("File not found: {}", path).into());
    }

    let path_lower = path.to_lowercase();

    if path_lower.ends_with(".h5ad") {
        inspect_h5ad(path)
    } else if path_lower.ends_with(".rds") {
        inspect_rds(path)
    } else {
        Err(format!("Unsupported file format. Expected .h5ad or .rds, got: {}", path).into())
    }
}

// ============================================================================
// validate_conversion
// ============================================================================

/// Validate conversion fidelity between two files.
///
/// Compares an original file with a converted file and returns a validation
/// report including Pearson R, MSE, RMSE, exact match percentage, and
/// optionally ARI/NMI for cluster labels.
///
/// @param original Path to the original file (.h5ad or .rds).
/// @param converted Path to the converted file (.h5ad or .rds).
/// @param tolerance Numeric. Numerical tolerance. Default 1e-7.
/// @param cluster_column Character or NULL. Column in obs for cluster labels. Default NULL.
/// @return A list with validation results.
/// @export
#[extendr]
fn validate_conversion(
    original: &str,
    converted: &str,
    tolerance: Option<f64>,
    cluster_column: Option<String>,
) -> Result<Robj> {
    if !Path::new(original).exists() {
        return Err(format!("Original file not found: {}", original).into());
    }
    if !Path::new(converted).exists() {
        return Err(format!("Converted file not found: {}", converted).into());
    }

    let tol = tolerance.unwrap_or(1e-7);

    let orig_data = load_ir_from_path(original)?;
    let conv_data = load_ir_from_path(converted)?;

    let report = crosscell::validation::roundtrip::validate_roundtrip(
        &orig_data, &conv_data, tol,
    );

    let passed = report.passed();
    let summary = report.summary();
    let detailed = report.detailed_report();

    if let Some(ref col) = cluster_column {
        let orig_labels = extract_string_labels_from_ir(&orig_data, col)?;
        let conv_labels = extract_string_labels_from_ir(&conv_data, col)?;
        let ari = crosscell::validation::calculate_ari(&orig_labels, &conv_labels);
        let nmi = crosscell::validation::calculate_nmi(&orig_labels, &conv_labels);

        R!("list(
            passed = {{passed}},
            summary = {{summary}},
            detailed_report = {{detailed}},
            ari = {{ari}},
            nmi = {{nmi}}
        )")
    } else {
        R!("list(
            passed = {{passed}},
            summary = {{summary}},
            detailed_report = {{detailed}}
        )")
    }
}

// ============================================================================
// Internal helpers
// ============================================================================

fn apply_transforms(
    data: &mut crosscell::ir::SingleCellData,
    normalize: bool,
    top_genes: Option<i32>,
    gene_id_column: Option<String>,
) -> Result<()> {
    if normalize {
        data.expression = crosscell::transform::normalize_library_size(&data.expression)
            .map_err(|e| Error::Other(format!("Normalization failed: {}", e)))?;
    }

    if let Some(n) = top_genes {
        if n <= 0 {
            return Err("top_genes must be a positive integer".into());
        }
        crosscell::transform::select_top_variable_genes(data, n as usize)
            .map_err(|e| Error::Other(format!("Gene selection failed: {}", e)))?;
    }

    if let Some(ref col) = gene_id_column {
        crosscell::transform::apply_gene_id_column(data, col)
            .map_err(|e| Error::Other(format!("Gene ID replacement failed: {}", e)))?;
    }

    Ok(())
}

fn load_ir_from_path(path: &str) -> Result<crosscell::ir::SingleCellData> {
    let path_lower = path.to_lowercase();
    if path_lower.ends_with(".h5ad") {
        crosscell::anndata::read_h5ad(path)
            .map_err(|e| Error::Other(format!("Failed to read H5AD: {}", e)))
    } else if path_lower.ends_with(".rds") {
        let result = crosscell::seurat::direct_read::read_seurat_direct(path, false)
            .map_err(|e| Error::Other(format!("Failed to read RDS: {}", e)))?;
        Ok(result.data)
    } else {
        Err(format!("Unsupported format: {}", path).into())
    }
}

fn extract_string_labels_from_ir(
    data: &crosscell::ir::SingleCellData,
    column: &str,
) -> Result<Vec<String>> {
    use arrow::array::{Array, StringArray, LargeStringArray, DictionaryArray};
    use arrow::datatypes::Int32Type;

    let col_idx = data
        .cell_metadata
        .column_index(column)
        .ok_or_else(|| Error::Other(format!(
            "Column '{}' not found in obs. Available: {:?}",
            column, data.cell_metadata.columns
        )))?;

    let col_data = &data.cell_metadata.data[col_idx];

    if let Some(arr) = col_data.as_any().downcast_ref::<StringArray>() {
        Ok((0..arr.len())
            .map(|i| if arr.is_null(i) { String::new() } else { arr.value(i).to_string() })
            .collect())
    } else if let Some(arr) = col_data.as_any().downcast_ref::<LargeStringArray>() {
        Ok((0..arr.len())
            .map(|i| if arr.is_null(i) { String::new() } else { arr.value(i).to_string() })
            .collect())
    } else if let Some(arr) = col_data.as_any().downcast_ref::<DictionaryArray<Int32Type>>() {
        let values = arr.values().as_any().downcast_ref::<StringArray>()
            .ok_or_else(|| Error::Other("Dictionary values are not strings".to_string()))?;
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    let key = arr.keys().value(i) as usize;
                    values.value(key).to_string()
                }
            })
            .collect())
    } else {
        Err(Error::Other(format!(
            "Column '{}' is not a string or categorical type",
            column
        )))
    }
}

fn inspect_h5ad(path: &str) -> Result<Robj> {
    let ir = crosscell::anndata::read_h5ad(path)
        .map_err(|e| Error::Other(format!("Failed to read H5AD: {}", e)))?;

    let format = "h5ad".to_string();
    let n_cells = ir.metadata.n_cells as i32;
    let n_genes = ir.metadata.n_genes as i32;
    let obs_cols: Vec<String> = ir.cell_metadata.columns.clone();
    let var_cols: Vec<String> = ir.gene_metadata.columns.clone();
    let has_embeddings = ir.embeddings.is_some();
    let has_layers = ir.layers.is_some();

    let embedding_names: Vec<String> = ir.embeddings
        .as_ref()
        .map(|e| e.keys().cloned().collect())
        .unwrap_or_default();

    let layer_names: Vec<String> = ir.layers
        .as_ref()
        .map(|l| l.keys().cloned().collect())
        .unwrap_or_default();

    R!("list(
        format = {{format}},
        n_cells = {{n_cells}},
        n_genes = {{n_genes}},
        obs_columns = {{obs_cols}},
        var_columns = {{var_cols}},
        has_embeddings = {{has_embeddings}},
        has_layers = {{has_layers}},
        embedding_names = {{embedding_names}},
        layer_names = {{layer_names}}
    )")
}

fn inspect_rds(path: &str) -> Result<Robj> {
    let result = crosscell::seurat::direct_read::read_seurat_direct(path, false)
        .map_err(|e| Error::Other(format!("Failed to read RDS: {}", e)))?;
    let ir = &result.data;

    let format = "rds".to_string();
    let n_cells = ir.metadata.n_cells as i32;
    let n_genes = ir.metadata.n_genes as i32;
    let obs_cols: Vec<String> = ir.cell_metadata.columns.clone();
    let var_cols: Vec<String> = ir.gene_metadata.columns.clone();
    let has_embeddings = ir.embeddings.is_some();
    let has_layers = ir.layers.is_some();

    let embedding_names: Vec<String> = ir.embeddings
        .as_ref()
        .map(|e| e.keys().cloned().collect())
        .unwrap_or_default();

    let layer_names: Vec<String> = ir.layers
        .as_ref()
        .map(|l| l.keys().cloned().collect())
        .unwrap_or_default();

    R!("list(
        format = {{format}},
        n_cells = {{n_cells}},
        n_genes = {{n_genes}},
        obs_columns = {{obs_cols}},
        var_columns = {{var_cols}},
        has_embeddings = {{has_embeddings}},
        has_layers = {{has_layers}},
        embedding_names = {{embedding_names}},
        layer_names = {{layer_names}}
    )")
}

// Macro to generate exports
extendr_module! {
    mod crosscell;
    fn read_h5ad_as_seurat;
    fn read_h5ad_as_sce;
    fn read_rds_fast;
    fn write_as_h5ad;
    fn write_rds_fast;
    fn inspect_file;
    fn validate_conversion;
}
