use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::Path;

use arrow::array::{Array, DictionaryArray, LargeStringArray, StringArray};
use arrow::datatypes::Int32Type;

pub mod convert_from_anndata;
pub mod convert_to_anndata;

/// CrossCell Python bindings module
#[pymodule]
fn _crosscell_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(read_rds, m)?)?;
    m.add_function(wrap_pyfunction!(read_h5ad, m)?)?;
    m.add_function(wrap_pyfunction!(write_h5ad, m)?)?;
    m.add_function(wrap_pyfunction!(write_rds, m)?)?;
    m.add_function(wrap_pyfunction!(inspect, m)?)?;
    m.add_function(wrap_pyfunction!(validate, m)?)?;
    Ok(())
}

/// Read a Seurat RDS file and return an AnnData object.
///
/// Parameters
/// ----------
/// path : str
///     Path to the .rds file containing a Seurat object.
/// normalize : bool, optional
///     Apply library-size normalization + log1p transform. Default False.
/// top_genes : int, optional
///     Select top N highly variable genes by variance. Default None (keep all).
/// gene_id_column : str, optional
///     Use specified column from var as gene identifiers. Default None.
/// keep_layers : bool, optional
///     Preserve Seurat V5 split layers as separate AnnData layers
///     instead of merging them. Default False.
///
/// Returns
/// -------
/// anndata.AnnData
///     The converted AnnData object.
///
/// Raises
/// ------
/// FileNotFoundError
///     If the file does not exist.
/// ValueError
///     If the file is not a valid Seurat object.
#[pyfunction]
#[pyo3(signature = (path, *, normalize=false, top_genes=None, gene_id_column=None, keep_layers=false))]
fn read_rds(
    py: Python<'_>,
    path: &str,
    normalize: bool,
    top_genes: Option<usize>,
    gene_id_column: Option<String>,
    keep_layers: bool,
) -> PyResult<PyObject> {
    let p = Path::new(path);
    if !p.exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "File not found: {}",
            path
        )));
    }

    let result = crosscell::seurat::read_seurat_direct(p, keep_layers).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to read Seurat RDS: {}", e))
    })?;

    let mut data = result.data;

    // Apply AI-Ready transforms in the same order as CLI
    if normalize {
        data.expression =
            crosscell::transform::normalize_library_size(&data.expression).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!("Normalization failed: {}", e))
            })?;
    }

    if let Some(n) = top_genes {
        crosscell::transform::select_top_variable_genes(&mut data, n).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Gene selection failed: {}", e))
        })?;
    }

    if let Some(ref col) = gene_id_column {
        crosscell::transform::apply_gene_id_column(&mut data, col).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Gene ID replacement failed: {}", e))
        })?;
    }

    convert_to_anndata::ir_to_anndata(py, &data)
}

/// Read an H5AD file and return an AnnData object.
///
/// Parameters
/// ----------
/// path : str
///     Path to the .h5ad file.
/// normalize : bool, optional
///     Apply library-size normalization + log1p transform. Default False.
/// top_genes : int, optional
///     Select top N highly variable genes by variance. Default None (keep all).
/// gene_id_column : str, optional
///     Use specified column from var as gene identifiers. Default None.
///
/// Returns
/// -------
/// anndata.AnnData
///     The loaded AnnData object.
///
/// Raises
/// ------
/// FileNotFoundError
///     If the file does not exist.
/// ValueError
///     If the file is not a valid H5AD file.
#[pyfunction]
#[pyo3(signature = (path, *, normalize=false, top_genes=None, gene_id_column=None))]
fn read_h5ad(
    py: Python<'_>,
    path: &str,
    normalize: bool,
    top_genes: Option<usize>,
    gene_id_column: Option<String>,
) -> PyResult<PyObject> {
    let p = Path::new(path);
    if !p.exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "File not found: {}",
            path
        )));
    }

    let mut data = crosscell::anndata::read_h5ad(p).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to read H5AD: {}", e))
    })?;

    // Apply transforms
    if normalize {
        data.expression =
            crosscell::transform::normalize_library_size(&data.expression).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!("Normalization failed: {}", e))
            })?;
    }

    if let Some(n) = top_genes {
        crosscell::transform::select_top_variable_genes(&mut data, n).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Gene selection failed: {}", e))
        })?;
    }

    if let Some(ref col) = gene_id_column {
        crosscell::transform::apply_gene_id_column(&mut data, col).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Gene ID replacement failed: {}", e))
        })?;
    }

    convert_to_anndata::ir_to_anndata(py, &data)
}

/// Write an AnnData object to an H5AD file.
///
/// Parameters
/// ----------
/// adata : anndata.AnnData
///     The AnnData object to write.
/// path : str
///     Output file path (.h5ad).
///
/// Raises
/// ------
/// ValueError
///     If the AnnData object cannot be converted.
#[pyfunction]
fn write_h5ad(py: Python<'_>, adata: &Bound<'_, PyAny>, path: &str) -> PyResult<()> {
    let ir = convert_from_anndata::anndata_to_ir(py, adata)?;

    crosscell::anndata::write_h5ad(&ir, path).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to write H5AD: {}", e))
    })
}

/// Write an AnnData object to a Seurat RDS file.
///
/// Parameters
/// ----------
/// adata : anndata.AnnData
///     The AnnData object to convert and write.
/// path : str
///     Output file path (.rds).
///
/// Raises
/// ------
/// ValueError
///     If the AnnData object is empty or cannot be converted.
#[pyfunction]
fn write_rds(py: Python<'_>, adata: &Bound<'_, PyAny>, path: &str) -> PyResult<()> {
    let ir = convert_from_anndata::anndata_to_ir(py, adata)?;

    // Validate non-empty
    if ir.metadata.n_cells == 0 || ir.metadata.n_genes == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Cannot write empty AnnData (0 cells or 0 genes) to RDS",
        ));
    }

    crosscell::seurat::write_seurat_rds(&ir, path).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to write Seurat RDS: {}", e))
    })
}

/// Inspect a file and return metadata information.
///
/// Parameters
/// ----------
/// path : str
///     Path to an .h5ad or .rds file.
///
/// Returns
/// -------
/// dict
///     File metadata including cell/gene counts, format info, etc.
///
/// Raises
/// ------
/// FileNotFoundError
///     If the file does not exist.
/// ValueError
///     If the file format is not supported.
#[pyfunction]
fn inspect(py: Python<'_>, path: &str) -> PyResult<PyObject> {
    let p = Path::new(path);
    if !p.exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "File not found: {}",
            path
        )));
    }

    let ext = p
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "h5ad" => inspect_h5ad_file(py, p),
        "rds" => inspect_rds_file(py, p),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unsupported file format: .{}. Supported: .h5ad, .rds",
            ext
        ))),
    }
}

/// Validate conversion fidelity between two files.
///
/// Compares an original file with a converted file and returns a detailed
/// validation report including Pearson R, MSE, RMSE, exact match %, and
/// optionally ARI/NMI for cluster labels.
///
/// Parameters
/// ----------
/// original : str
///     Path to the original file (.h5ad or .rds).
/// converted : str
///     Path to the converted file (.h5ad or .rds).
/// tolerance : float, optional
///     Numerical tolerance for floating-point comparisons. Default 1e-7.
/// cluster_column : str, optional
///     Column name in obs for cluster labels. When provided, ARI and NMI
///     are calculated. Default None.
///
/// Returns
/// -------
/// dict
///     Validation report with keys: passed (bool), summary (str),
///     detailed_report (str), and optionally ari/nmi (float).
///
/// Raises
/// ------
/// FileNotFoundError
///     If either file does not exist.
/// ValueError
///     If file formats are unsupported or validation fails to run.
#[pyfunction]
#[pyo3(signature = (original, converted, *, tolerance=1e-7, cluster_column=None))]
fn validate(
    py: Python<'_>,
    original: &str,
    converted: &str,
    tolerance: f64,
    cluster_column: Option<String>,
) -> PyResult<PyObject> {
    // Check files exist
    if !Path::new(original).exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "Original file not found: {}",
            original
        )));
    }
    if !Path::new(converted).exists() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "Converted file not found: {}",
            converted
        )));
    }

    let orig_ext = Path::new(original)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();
    let conv_ext = Path::new(converted)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    // Load data based on format
    let orig_data = load_ir(original, &orig_ext)?;
    let conv_data = load_ir(converted, &conv_ext)?;

    // Run validation
    let report =
        crosscell::validation::roundtrip::validate_roundtrip(&orig_data, &conv_data, tolerance);

    let dict = PyDict::new_bound(py);
    dict.set_item("passed", report.passed())?;
    dict.set_item("summary", report.summary())?;
    dict.set_item("detailed_report", report.detailed_report())?;

    // Calculate ARI/NMI if cluster_column provided
    if let Some(ref col) = cluster_column {
        let orig_labels = extract_labels(&orig_data, col)?;
        let conv_labels = extract_labels(&conv_data, col)?;
        let ari = crosscell::validation::calculate_ari(&orig_labels, &conv_labels);
        let nmi = crosscell::validation::calculate_nmi(&orig_labels, &conv_labels);
        dict.set_item("ari", ari)?;
        dict.set_item("nmi", nmi)?;
    }

    Ok(dict.into())
}

// ============================================================================
// Helper functions
// ============================================================================

fn load_ir(path: &str, ext: &str) -> PyResult<crosscell::ir::SingleCellData> {
    match ext {
        "h5ad" => crosscell::anndata::read_h5ad(path).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Failed to read H5AD: {}", e))
        }),
        "rds" => {
            let result = crosscell::seurat::read_seurat_direct(path, false).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!("Failed to read RDS: {}", e))
            })?;
            Ok(result.data)
        }
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unsupported format: .{}",
            ext
        ))),
    }
}

fn extract_labels(data: &crosscell::ir::SingleCellData, column: &str) -> PyResult<Vec<String>> {
    let col_idx = data.cell_metadata.column_index(column).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "Column '{}' not found in obs. Available: {:?}",
            column, data.cell_metadata.columns
        ))
    })?;

    let col_data = &data.cell_metadata.data[col_idx];

    if let Some(arr) = col_data.as_any().downcast_ref::<StringArray>() {
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            })
            .collect())
    } else if let Some(arr) = col_data.as_any().downcast_ref::<LargeStringArray>() {
        Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            })
            .collect())
    } else if let Some(arr) = col_data
        .as_any()
        .downcast_ref::<DictionaryArray<Int32Type>>()
    {
        let values = arr
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err("Dictionary values are not strings")
            })?;
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
        Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Column '{}' is not a string or categorical type",
            column
        )))
    }
}

fn inspect_h5ad_file(py: Python<'_>, path: &Path) -> PyResult<PyObject> {
    let info = crosscell::anndata::inspect_h5ad(path).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to inspect H5AD: {}", e))
    })?;

    let dict = PyDict::new_bound(py);
    dict.set_item("format", "h5ad")?;
    dict.set_item("n_cells", info.n_cells)?;
    dict.set_item("n_genes", info.n_genes)?;
    dict.set_item("is_sparse", info.is_sparse)?;
    dict.set_item("sparse_format", info.sparse_format)?;
    dict.set_item("nnz", info.nnz)?;
    dict.set_item("estimated_matrix_size", info.estimated_matrix_size)?;
    dict.set_item("n_obs_columns", info.n_obs_columns)?;
    dict.set_item("n_var_columns", info.n_var_columns)?;
    dict.set_item("embedding_names", info.embedding_names)?;
    dict.set_item("layer_names", info.layer_names)?;
    dict.set_item("obsp_names", info.obsp_names)?;
    dict.set_item("varp_names", info.varp_names)?;
    dict.set_item("has_spatial", info.has_spatial)?;
    dict.set_item("file_size", info.file_size)?;
    Ok(dict.into())
}

fn inspect_rds_file(py: Python<'_>, path: &Path) -> PyResult<PyObject> {
    let result = crosscell::seurat::read_seurat_direct(path, false).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Failed to inspect RDS: {}", e))
    })?;

    let dict = PyDict::new_bound(py);
    dict.set_item("format", "rds")?;
    dict.set_item("n_cells", result.data.metadata.n_cells)?;
    dict.set_item("n_genes", result.data.metadata.n_genes)?;
    dict.set_item("seurat_version", result.version.to_string())?;

    let (n_rows, n_cols) = result.data.expression.shape();
    let is_sparse = !matches!(
        result.data.expression,
        crosscell::ir::ExpressionMatrix::Dense(_)
    );
    dict.set_item("is_sparse", is_sparse)?;
    dict.set_item("matrix_shape", (n_rows, n_cols))?;

    // Embedding names
    let embedding_names: Vec<String> = result
        .data
        .embeddings
        .as_ref()
        .map(|e| e.keys().cloned().collect())
        .unwrap_or_default();
    dict.set_item("embedding_names", embedding_names)?;

    // Layer names
    let layer_names: Vec<String> = result
        .data
        .layers
        .as_ref()
        .map(|l| l.keys().cloned().collect())
        .unwrap_or_default();
    dict.set_item("layer_names", layer_names)?;

    // Skipped components
    if result.skipped.has_skipped() {
        dict.set_item("skipped_components", result.skipped.total_skipped())?;
    }

    Ok(dict.into())
}
