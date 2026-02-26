//! SingleCellExperiment to IR conversion
//!
//! This module provides functions to convert SingleCellExperiment objects to CrossCell IR format.

use crate::ir::{
    DataFrame, DatasetMetadata, Embedding, ExpressionMatrix, SingleCellData, SparseMatrixCSR,
};
use crate::rds::{parse_rds, RObject, RdsFile};
use crate::sce::error::{Result, SceError};
use crate::sce::extract::{
    ensure_sce_class, extract_integer_vector, extract_real_vector, extract_slot,
    extract_slot_optional, get_named_list_items,
};
use crate::sparse::convert::csr_to_csc;
use arrow::datatypes::Int32Type;
use std::collections::HashMap;
use std::sync::Arc;

/// Read SingleCellExperiment RDS file and convert to IR
///
/// This function reads a SingleCellExperiment object and converts it to CrossCell IR format.
pub fn sce_rds_to_ir(path: &str) -> Result<SingleCellData> {
    use std::path::Path;

    let file = parse_rds(Path::new(path))
        .map_err(|e| SceError::ParseError(format!("Failed to read RDS: {}", e)))?;

    // Parse SCE structure
    parse_sce(&file.object, &file)
}

/// Parse SingleCellExperiment object to IR
fn parse_sce(robj: &RObject, file: &RdsFile) -> Result<SingleCellData> {
    // Verify this is a SCE object
    ensure_sce_class(robj)?;

    log::debug!("Parsing SingleCellExperiment object");

    // Extract assays (required)
    let assays_slot = extract_slot(robj, "assays", file)?;
    let (expression, layers) = parse_assays(assays_slot, file)?;
    let (n_cells, n_genes) = expression.shape();

    log::debug!(
        "Extracted expression matrix: {} cells × {} genes",
        n_cells,
        n_genes
    );

    // Extract colData (cell metadata)
    let col_data_slot = extract_slot_optional(robj, "colData", file)?;
    let cell_metadata = if let Some(col_data) = col_data_slot {
        parse_col_data(col_data, n_cells, file)?
    } else {
        DataFrame::empty(n_cells)
    };

    log::debug!(
        "Extracted cell metadata: {} columns",
        cell_metadata.n_cols()
    );

    // Extract rowData (gene metadata) - stored in elementMetadata slot
    let row_data_slot = extract_slot_optional(robj, "elementMetadata", file)?;
    let gene_metadata = if let Some(row_data) = row_data_slot {
        parse_row_data(row_data, n_genes, file)?
    } else {
        DataFrame::empty(n_genes)
    };

    log::debug!(
        "Extracted gene metadata: {} columns",
        gene_metadata.n_cols()
    );

    // Extract reducedDims (from int_colData slot)
    let int_col_data_slot = extract_slot_optional(robj, "int_colData", file)?;
    let embeddings = if let Some(int_col_data) = int_col_data_slot {
        parse_reduced_dims(int_col_data, n_cells, file)?
    } else {
        None
    };

    if let Some(ref embs) = embeddings {
        log::debug!("Extracted {} embeddings", embs.len());
    }

    // Create dataset metadata
    let mut metadata = DatasetMetadata::new(n_cells, n_genes, "sce".to_string());
    metadata.source_version = Some("1.0".to_string());

    // Create SingleCellData
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .map_err(|e| SceError::ValidationError(e))?;

    // Add optional fields
    data.layers = layers;
    data.embeddings = embeddings;

    Ok(data)
}

/// Parse assays slot from SCE
///
/// SCE assays are stored in a SimpleList within an Assays S4 object
/// Structure: assays (SimpleAssays) -> data (SimpleList) -> listData (List)
fn parse_assays(
    assays_slot: &RObject,
    _file: &RdsFile,
) -> Result<(ExpressionMatrix, Option<HashMap<String, ExpressionMatrix>>)> {
    // assays slot is an S4 object (SimpleAssays class) containing a SimpleList
    let simple_list = match assays_slot {
        RObject::S4Object(s4) => {
            // Get the data slot which contains the SimpleList
            s4.attributes
                .get("data")
                .ok_or_else(|| SceError::InvalidAssaysStructure)?
        }
        RObject::GenericVector(_) => assays_slot,
        _ => return Err(SceError::InvalidAssaysStructure),
    };

    // Now get the listData from the SimpleList
    let assay_data = match simple_list {
        RObject::S4Object(s4) => {
            // SimpleList has listData attribute
            s4.attributes
                .get("listData")
                .ok_or_else(|| SceError::InvalidAssaysStructure)?
        }
        RObject::GenericVector(_) => simple_list,
        _ => return Err(SceError::InvalidAssaysStructure),
    };

    // Get named list of assays
    let assay_list =
        get_named_list_items(assay_data).ok_or_else(|| SceError::InvalidAssaysStructure)?;

    if assay_list.is_empty() {
        return Err(SceError::MissingAssay("No assays found".to_string()));
    }

    log::debug!(
        "Found {} assays: {:?}",
        assay_list.len(),
        assay_list.iter().map(|(n, _)| *n).collect::<Vec<_>>()
    );

    // Find primary assay (counts or first available)
    let primary_assay = assay_list
        .iter()
        .find(|(name, _)| *name == "counts")
        .or_else(|| assay_list.iter().find(|(name, _)| *name == "logcounts"))
        .or_else(|| assay_list.first())
        .ok_or_else(|| SceError::MissingAssay("No valid assay found".to_string()))?;

    let (primary_name, primary_matrix) = primary_assay;
    log::debug!("Using '{}' as primary assay", primary_name);

    let expression = parse_sparse_matrix(primary_matrix)?;

    // Parse additional assays as layers
    let layers = if assay_list.len() > 1 {
        let mut layer_map = HashMap::new();
        for (name, matrix) in &assay_list {
            if *name != *primary_name {
                match parse_sparse_matrix(matrix) {
                    Ok(layer_expr) => {
                        layer_map.insert(name.to_string(), layer_expr);
                    }
                    Err(e) => {
                        log::warn!("Failed to parse layer '{}': {}", name, e);
                    }
                }
            }
        }
        if layer_map.is_empty() {
            None
        } else {
            Some(layer_map)
        }
    } else {
        None
    };

    Ok((expression, layers))
}

/// Parse sparse matrix (dgCMatrix) from SCE
fn parse_sparse_matrix(robj: &RObject) -> Result<ExpressionMatrix> {
    // dgCMatrix is an S4 object with slots: i, p, x, Dim
    let (i, p, x, dim) =
        match robj {
            RObject::S4Object(s4) => {
                let i = s4.attributes.get("i").ok_or_else(|| {
                    SceError::ParseError("Missing 'i' slot in dgCMatrix".to_string())
                })?;
                let p = s4.attributes.get("p").ok_or_else(|| {
                    SceError::ParseError("Missing 'p' slot in dgCMatrix".to_string())
                })?;
                let x = s4.attributes.get("x").ok_or_else(|| {
                    SceError::ParseError("Missing 'x' slot in dgCMatrix".to_string())
                })?;
                let dim = s4.attributes.get("Dim").ok_or_else(|| {
                    SceError::ParseError("Missing 'Dim' slot in dgCMatrix".to_string())
                })?;
                (i, p, x, dim)
            }
            _ => {
                return Err(SceError::ParseError(format!(
                    "Expected S4 object for dgCMatrix, got {}",
                    robj.type_name()
                )))
            }
        };

    // Extract vectors
    let r_row_indices = extract_integer_vector(i)?;
    let r_col_ptrs = extract_integer_vector(p)?;
    let data = extract_real_vector(x)?;
    let dimensions = extract_integer_vector(dim)?;

    if dimensions.len() != 2 {
        return Err(SceError::ParseError(format!(
            "Expected 2D dimensions, got {}",
            dimensions.len()
        )));
    }

    // R dgCMatrix dimensions: [n_genes, n_cells] (genes as rows, cells as columns)
    let n_genes = dimensions[0] as usize;
    let n_cells = dimensions[1] as usize;

    log::debug!(
        "dgCMatrix dimensions from R: {} genes × {} cells",
        n_genes,
        n_cells
    );

    // Convert to usize
    let r_row_indices: Vec<usize> = r_row_indices.iter().map(|&x| x as usize).collect();
    let r_col_ptrs: Vec<usize> = r_col_ptrs.iter().map(|&x| x as usize).collect();

    // Interpret R dgCMatrix as CSR (cells × genes)
    // R stores genes×cells in CSC format, which is equivalent to cells×genes in CSR
    let csr = SparseMatrixCSR {
        n_rows: n_cells,
        n_cols: n_genes,
        indptr: r_col_ptrs,
        indices: r_row_indices,
        data,
    };

    log::debug!(
        "Interpreted as CSR: {} cells × {} genes",
        csr.n_rows,
        csr.n_cols
    );

    // Convert CSR to CSC
    let csc = csr_to_csc(&csr);

    log::debug!(
        "Converted to CSC: {} cells × {} genes",
        csc.n_rows,
        csc.n_cols
    );

    Ok(ExpressionMatrix::SparseCSC(csc))
}

/// Parse colData (cell metadata) from SCE
fn parse_col_data(col_data: &RObject, expected_rows: usize, _file: &RdsFile) -> Result<DataFrame> {
    // colData is a DFrame (S4 object)
    let list_data = match col_data {
        RObject::S4Object(s4) => s4
            .attributes
            .get("listData")
            .or_else(|| s4.attributes.get("data")),
        RObject::GenericVector(_) => Some(col_data),
        _ => None,
    };

    let list_data = match list_data {
        Some(data) => data,
        None => return Ok(DataFrame::empty(expected_rows)),
    };

    // Get named columns
    let columns = get_named_list_items(list_data).unwrap_or_default();

    if columns.is_empty() {
        return Ok(DataFrame::empty(expected_rows));
    }

    // Convert R columns to Arrow arrays
    use arrow::array::ArrayRef;

    let mut col_names: Vec<String> = Vec::new();
    let mut col_arrays: Vec<ArrayRef> = Vec::new();

    for (col_name, col_value) in columns {
        // Skip internal columns
        if col_name.starts_with('_') {
            continue;
        }

        let array = robject_to_arrow_array(col_value, col_name, expected_rows)?;
        col_names.push(col_name.to_string());
        col_arrays.push(array);
    }

    if col_names.is_empty() {
        return Ok(DataFrame::empty(expected_rows));
    }

    let df = DataFrame::new(col_names, col_arrays, expected_rows)
        .map_err(|e| SceError::ParseError(format!("Failed to create DataFrame: {}", e)))?;

    Ok(df)
}

/// Parse rowData (gene metadata) from SCE
fn parse_row_data(row_data: &RObject, expected_rows: usize, file: &RdsFile) -> Result<DataFrame> {
    // rowData is stored in elementMetadata slot, which is a DFrame
    parse_col_data(row_data, expected_rows, file)
}

/// Parse reducedDims from SCE (stored in int_colData)
fn parse_reduced_dims(
    int_col_data: &RObject,
    expected_rows: usize,
    file: &RdsFile,
) -> Result<Option<HashMap<String, Embedding>>> {
    // int_colData is a DFrame containing reducedDims
    let list_data = match int_col_data {
        RObject::S4Object(s4) => s4
            .attributes
            .get("listData")
            .or_else(|| s4.attributes.get("data")),
        RObject::GenericVector(_) => Some(int_col_data),
        _ => None,
    };

    let list_data = match list_data {
        Some(data) => data,
        None => return Ok(None),
    };

    // Look for reducedDims column
    let columns = get_named_list_items(list_data).unwrap_or_default();

    for (col_name, col_value) in &columns {
        if *col_name == "reducedDims" {
            return parse_reduced_dims_list(col_value, expected_rows, file);
        }
    }

    Ok(None)
}

/// Parse reducedDims SimpleList
fn parse_reduced_dims_list(
    reduced_dims: &RObject,
    expected_rows: usize,
    _file: &RdsFile,
) -> Result<Option<HashMap<String, Embedding>>> {
    // reducedDims is a SimpleList (S4) containing matrices
    let list_data = match reduced_dims {
        RObject::S4Object(s4) => s4
            .attributes
            .get("listData")
            .or_else(|| s4.attributes.get("data")),
        RObject::GenericVector(_) => Some(reduced_dims),
        _ => None,
    };

    let list_data = match list_data {
        Some(data) => data,
        None => return Ok(None),
    };

    let dim_list = get_named_list_items(list_data).unwrap_or_default();

    if dim_list.is_empty() {
        return Ok(None);
    }

    let mut embeddings = HashMap::new();

    for (name, matrix) in dim_list {
        match parse_dense_matrix(matrix, expected_rows) {
            Ok(embedding) => {
                embeddings.insert(name.to_string(), embedding);
            }
            Err(e) => {
                log::warn!("Failed to parse reducedDim '{}': {}", name, e);
            }
        }
    }

    if embeddings.is_empty() {
        Ok(None)
    } else {
        Ok(Some(embeddings))
    }
}

/// Parse dense matrix (R matrix) to Embedding
fn parse_dense_matrix(robj: &RObject, _expected_rows: usize) -> Result<Embedding> {
    // R matrix is stored as a real vector with dim attribute
    let (data, dim) = match robj {
        RObject::DoubleVector(dv) => {
            // Check for dim attribute
            if let Some(dim_obj) = dv.attributes.get("dim") {
                let dim = extract_integer_vector(dim_obj)?;
                (dv.data.clone(), dim)
            } else {
                // No dim attribute, assume it's a vector (n x 1 matrix)
                let n = dv.data.len();
                (dv.data.clone(), vec![n as i32, 1])
            }
        }
        _ => {
            return Err(SceError::ParseError(format!(
                "Expected real vector or matrix, got {}",
                robj.type_name()
            )))
        }
    };

    if dim.len() != 2 {
        return Err(SceError::ParseError(format!(
            "Expected 2D dimensions, got {}",
            dim.len()
        )));
    }

    let n_rows = dim[0] as usize;
    let n_cols = dim[1] as usize;

    if data.len() != n_rows * n_cols {
        return Err(SceError::ParseError(format!(
            "Matrix data length {} doesn't match dimensions {} x {}",
            data.len(),
            n_rows,
            n_cols
        )));
    }

    // R stores matrices in column-major order, convert to row-major
    let mut row_major_data = vec![0.0; n_rows * n_cols];
    for col in 0..n_cols {
        for row in 0..n_rows {
            row_major_data[row * n_cols + col] = data[col * n_rows + row];
        }
    }

    // Create Embedding
    let embedding = Embedding::new("embedding".to_string(), row_major_data, n_rows, n_cols)
        .map_err(|e| SceError::ValidationError(e))?;

    Ok(embedding)
}

/// Convert RObject to Arrow array
fn robject_to_arrow_array(
    robj: &RObject,
    col_name: &str,
    expected_len: usize,
) -> Result<Arc<dyn arrow::array::Array>> {
    use arrow::array::*;

    match robj {
        RObject::IntegerVector(iv) => Ok(Arc::new(Int32Array::from(iv.data.clone()))),
        RObject::DoubleVector(dv) => Ok(Arc::new(Float64Array::from(dv.data.clone()))),
        RObject::StringVector(sv) => Ok(Arc::new(StringArray::from(sv.data.clone()))),
        RObject::LogicalVector(lv) => {
            // Convert i32 to bool
            let bools: Vec<bool> = lv.data.iter().map(|&v| v != 0).collect();
            Ok(Arc::new(BooleanArray::from(bools)))
        }
        _ => {
            // Check if it has attributes (might be a factor)
            if let Some(attrs) = robj.attributes() {
                let has_levels = attrs.get("levels").is_some();
                let is_factor = attrs
                    .get_class()
                    .map(|c| c.iter().any(|s| s == "factor"))
                    .unwrap_or(false);

                if has_levels || is_factor {
                    return robject_factor_to_arrow(robj, col_name, expected_len);
                }
            }

            // Unsupported type - create string array with placeholder values
            log::warn!(
                "Unsupported RObject type for column {}: {:?}, using placeholder",
                col_name,
                robj.type_name()
            );
            let placeholder: Vec<String> =
                (0..expected_len).map(|i| format!("Value_{}", i)).collect();
            Ok(Arc::new(StringArray::from(placeholder)))
        }
    }
}

/// Convert R factor to Arrow dictionary array
fn robject_factor_to_arrow(
    robj: &RObject,
    col_name: &str,
    _expected_len: usize,
) -> Result<Arc<dyn arrow::array::Array>> {
    use arrow::array::*;

    // Get the integer codes (1-based in R)
    let codes = match robj {
        RObject::IntegerVector(iv) => &iv.data,
        _ => {
            return Err(SceError::ParseError(format!(
                "Expected integer codes for factor column {}",
                col_name
            )));
        }
    };

    // Get the levels from attributes
    let levels = if let Some(attrs) = robj.attributes() {
        if let Some(RObject::StringVector(sv)) = attrs.get("levels") {
            sv.data.clone()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // Convert 1-based R codes to 0-based, handling NA
    let keys: Vec<Option<i32>> = codes
        .iter()
        .map(|&v| {
            if v <= 0 || v == i32::MIN {
                None // NA
            } else {
                Some(v - 1) // Convert 1-based to 0-based
            }
        })
        .collect();

    // Create dictionary array
    let keys_array = Int32Array::from(keys);
    let values_array = StringArray::from(levels);

    let dict_array = DictionaryArray::<Int32Type>::try_new(keys_array, Arc::new(values_array))
        .map_err(|e| SceError::ParseError(format!("Failed to create dictionary array: {}", e)))?;

    Ok(Arc::new(dict_array))
}
