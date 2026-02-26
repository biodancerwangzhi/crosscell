//! Seurat to IR conversion
//!
//! This module provides functions to convert simplified Seurat objects to CrossCell IR format.

use crate::ir::{
    DataFrame, DatasetMetadata, Embedding, ExpressionMatrix, PairwiseMatrix, SingleCellData,
    SpatialData,
};
use crate::rds::{parse_rds, RObject, RdsFile};
use crate::seurat::error::SeuratError;
use crate::seurat::extract::get_named_list_items;
use arrow::datatypes::Int32Type;
use std::collections::HashMap;
use std::sync::Arc;

/// Read simplified Seurat RDS file and convert to IR
///
/// This function reads a simplified Seurat object (created by simplify_seurat.R)
/// and converts it to CrossCell IR format.
pub fn seurat_rds_to_ir(path: &str) -> Result<SingleCellData, SeuratError> {
    use std::path::Path;

    let file = parse_rds(Path::new(path))
        .map_err(|e| SeuratError::ParseError(format!("Failed to read RDS: {}", e)))?;

    // Parse simplified Seurat structure
    parse_simplified_seurat(&file.object, &file)
}

/// Parse simplified Seurat object to IR
fn parse_simplified_seurat(robj: &RObject, file: &RdsFile) -> Result<SingleCellData, SeuratError> {
    // Get named list items from the object
    let fields = get_fields_from_object(robj, file)?;

    // Extract fields
    let mut _project_name = None;
    let mut active_assay = None;
    let mut assays = None;
    let mut meta_data = None;
    let mut reductions = None;
    let mut images = None;
    let mut graphs = None;

    for (name, value) in &fields {
        match *name {
            "project_name" | "project.name" => {
                _project_name = Some(extract_string_scalar(value)?);
            }
            "active_assay" | "active.assay" => {
                active_assay = Some(extract_string_scalar(value)?);
            }
            "assays" => {
                assays = Some(*value);
            }
            "meta_data" | "meta.data" => {
                meta_data = Some(*value);
            }
            "reductions" => {
                reductions = Some(*value);
            }
            "images" => {
                images = Some(*value);
            }
            "graphs" => {
                graphs = Some(*value);
            }
            "class" | "active_ident" | "active.ident" | "neighbors" | "misc" | "version"
            | "commands" | "tools" => {
                // Skip these fields
            }
            _ => {
                // Unknown field, skip
            }
        }
    }

    // Validate required fields
    let assays = assays
        .ok_or_else(|| SeuratError::ParseError("Missing required field: assays".to_string()))?;
    let meta_data = meta_data
        .ok_or_else(|| SeuratError::ParseError("Missing required field: meta_data".to_string()))?;

    // Parse assays (extract first assay as main expression matrix)
    let (expression, gene_metadata, layers) = parse_assays(assays, file)?;
    let (n_cells, n_genes) = expression.shape();

    // Parse metadata
    let cell_metadata = parse_metadata(meta_data, n_cells, file)?;

    // Parse reductions (if present) — also extracts gene loadings
    let (embeddings, gene_loadings) = if let Some(red) = reductions {
        let embs = parse_reductions(red, n_cells, file)?;
        let loads = parse_reduction_loadings(red, file);
        (
            Some(embs),
            if loads.is_empty() { None } else { Some(loads) },
        )
    } else {
        (None, None)
    };

    // Parse graphs → cell_pairwise (if present)
    let cell_pairwise = if let Some(g) = graphs {
        parse_graphs(g, n_cells, file).unwrap_or(None)
    } else {
        None
    };

    // Parse spatial data (if present and non-empty)
    let spatial = if let Some(imgs) = images {
        // Check if images is an empty list/S4 before trying to parse
        let is_empty = match imgs {
            RObject::GenericVector(gv) => gv.data.is_empty(),
            RObject::S4Object(s4) => s4.attributes.is_empty(),
            RObject::Null => true,
            _ => false,
        };
        if is_empty {
            None
        } else {
            Some(parse_spatial_data(imgs, n_cells, file)?)
        }
    } else {
        None
    };

    // Create dataset metadata
    let mut metadata = DatasetMetadata::new(n_cells, n_genes, "seurat".to_string());
    metadata.source_version = Some("5.0".to_string());
    metadata.active_assay = active_assay;

    // Create SingleCellData
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .map_err(|e| SeuratError::ValidationError(e))?;

    // Add optional fields
    data.layers = layers;
    data.embeddings = embeddings;
    data.spatial = spatial;
    data.gene_loadings = gene_loadings;
    data.cell_pairwise = cell_pairwise;

    Ok(data)
}

/// Get fields from an RObject (handles List with names, PairList, etc.)
fn get_fields_from_object<'a>(
    robj: &'a RObject,
    file: &'a RdsFile,
) -> Result<Vec<(&'a str, &'a RObject)>, SeuratError> {
    match robj {
        RObject::GenericVector(gv) => {
            if let Some(names) = gv.attributes.get_names() {
                Ok(names
                    .iter()
                    .zip(gv.data.iter())
                    .map(|(n, v)| (n.as_str(), v))
                    .collect())
            } else if gv.data.len() == 1 {
                // Single element list without names - recurse
                get_fields_from_object(&gv.data[0], file)
            } else {
                Err(SeuratError::ParseError(
                    "Expected named list for SimplifiedSeurat".to_string(),
                ))
            }
        }
        RObject::PairList(pl) => Ok(pl
            .tag_names
            .iter()
            .zip(pl.data.iter())
            .map(|(n, v)| (n.as_str(), v))
            .collect()),
        RObject::S4Object(s4) => {
            // S4 objects store slots as attributes
            Ok(s4
                .attributes
                .names
                .iter()
                .zip(s4.attributes.values.iter())
                .map(|(n, v)| (n.as_str(), v))
                .collect())
        }
        _ => Err(SeuratError::ParseError(format!(
            "Expected named list, got {}",
            robj.type_name()
        ))),
    }
}

/// Extract string scalar from RObject
fn extract_string_scalar(robj: &RObject) -> Result<String, SeuratError> {
    match robj {
        RObject::StringVector(sv) => {
            if sv.data.is_empty() {
                Err(SeuratError::ParseError("Empty string vector".to_string()))
            } else {
                Ok(sv.data[0].clone())
            }
        }
        _ => Err(SeuratError::ParseError(format!(
            "Expected string, got {}",
            robj.type_name()
        ))),
    }
}

/// Parse assays field
fn parse_assays(
    robj: &RObject,
    file: &RdsFile,
) -> Result<
    (
        ExpressionMatrix,
        DataFrame,
        Option<HashMap<String, ExpressionMatrix>>,
    ),
    SeuratError,
> {
    let assay_list = get_fields_from_object(robj, file)?;

    if assay_list.is_empty() {
        return Err(SeuratError::ParseError(
            "No assays found in simplified Seurat".to_string(),
        ));
    }

    // Extract first assay as main expression matrix
    let (_first_name, first_assay) = &assay_list[0];
    let (expression, gene_metadata) = parse_single_assay(first_assay, file)?;

    // If there are multiple assays, store others as layers
    let layers = if assay_list.len() > 1 {
        let mut layer_map = HashMap::new();
        for (name, assay) in assay_list.iter().skip(1) {
            let (layer_expr, _) = parse_single_assay(assay, file)?;
            layer_map.insert(name.to_string(), layer_expr);
        }
        Some(layer_map)
    } else {
        None
    };

    Ok((expression, gene_metadata, layers))
}

/// Parse single assay
fn parse_single_assay(
    robj: &RObject,
    file: &RdsFile,
) -> Result<(ExpressionMatrix, DataFrame), SeuratError> {
    use crate::seurat::assay5::{
        detect_assay_type, extract_assay5_gene_metadata, extract_assay5_main_matrix, AssayType,
    };

    // 检Assay 类型
    let assay_type = detect_assay_type(robj);
    log::debug!("Detected assay type: {:?}", assay_type);

    match assay_type {
        AssayType::Assay5 => {
            // 使用 Assay5 提取逻辑
            log::debug!("Using Assay5 extraction logic");
            let expression = extract_assay5_main_matrix(robj, file)?;
            let (n_cells, n_genes) = expression.shape();
            log::debug!("Assay5 matrix: {} cells × {} genes", n_cells, n_genes);

            let gene_metadata = extract_assay5_gene_metadata(robj, n_genes, file)?;
            Ok((expression, gene_metadata))
        }
        AssayType::Legacy | AssayType::Unknown => {
            // 使用传统 Assay 提取逻辑
            log::debug!("Using legacy Assay extraction logic");
            parse_legacy_assay(robj, file)
        }
    }
}

/// Parse legacy (V3/V4) assay
fn parse_legacy_assay(
    robj: &RObject,
    file: &RdsFile,
) -> Result<(ExpressionMatrix, DataFrame), SeuratError> {
    let assay_fields = get_fields_from_object(robj, file)?;

    let mut counts = None;
    let mut features = None;
    let mut _cells = None;
    let mut _assay_class = None;

    for (name, value) in &assay_fields {
        match *name {
            "class" => {
                _assay_class = Some(extract_string_scalar(value)?);
            }
            "counts" => {
                counts = Some(*value);
            }
            "layers" => {
                // For Assay5, counts is in layers$counts
                if let Some(layer_fields) = get_named_list_items(value) {
                    for (layer_name, layer_value) in layer_fields {
                        if layer_name == "counts" {
                            counts = Some(layer_value);
                            break;
                        }
                    }
                }
            }
            "features" => {
                features = Some(*value);
            }
            "cells" => {
                _cells = Some(*value);
            }
            _ => {}
        }
    }

    // Extract counts matrix
    let counts = counts
        .ok_or_else(|| SeuratError::ParseError("Missing counts matrix in assay".to_string()))?;
    let expression = parse_dgcmatrix(counts)?;

    let (n_cells, n_genes) = expression.shape();
    log::debug!("After transpose: {} cells × {} genes", n_cells, n_genes);

    // Extract feature names (if present)
    let _feature_names = if let Some(features) = features {
        let names = extract_string_vector(features)?;
        log::debug!("Feature names count: {}", names.len());
        Some(names)
    } else {
        log::debug!("No feature names provided");
        None
    };

    // Create gene metadata DataFrame with correct number of rows
    let gene_metadata = DataFrame::empty(n_genes);
    log::debug!("Created gene_metadata with {} rows", n_genes);

    Ok((expression, gene_metadata))
}

/// Parse dgCMatrix (R sparse matrix) to ExpressionMatrix
pub fn parse_dgcmatrix(robj: &RObject) -> Result<ExpressionMatrix, SeuratError> {
    // dgCMatrix is an S4 object with slots: i, p, x, Dim
    let (i, p, x, dim) = match robj {
        RObject::S4Object(s4) => {
            let i = s4.attributes.get("i").ok_or_else(|| {
                SeuratError::ParseError("Missing 'i' slot in dgCMatrix".to_string())
            })?;
            let p = s4.attributes.get("p").ok_or_else(|| {
                SeuratError::ParseError("Missing 'p' slot in dgCMatrix".to_string())
            })?;
            let x = s4.attributes.get("x").ok_or_else(|| {
                SeuratError::ParseError("Missing 'x' slot in dgCMatrix".to_string())
            })?;
            let dim = s4.attributes.get("Dim").ok_or_else(|| {
                SeuratError::ParseError("Missing 'Dim' slot in dgCMatrix".to_string())
            })?;
            (i, p, x, dim)
        }
        _ => {
            return Err(SeuratError::ParseError(format!(
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
        return Err(SeuratError::ParseError(format!(
            "Expected 2D dimensions, got {}",
            dimensions.len()
        )));
    }

    // R dgCMatrix dimensions: [n_genes, n_cells]
    let n_genes = dimensions[0] as usize;
    let n_cells = dimensions[1] as usize;

    log::debug!(
        "dgCMatrix dimensions from R: {} genes × {} cells",
        n_genes,
        n_cells
    );
    log::debug!(
        "R col_ptrs length: {} (expected {})",
        r_col_ptrs.len(),
        n_cells + 1
    );

    // Convert to usize
    let r_row_indices: Vec<usize> = r_row_indices.iter().map(|&x| x as usize).collect();
    let r_col_ptrs: Vec<usize> = r_col_ptrs.iter().map(|&x| x as usize).collect();

    use crate::ir::SparseMatrixCSR;
    use crate::sparse::convert::csr_to_csc;

    // Interpret R dgCMatrix as CSR (cells × genes)
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
    log::debug!(
        "CSC indptr length: {} (expected {})",
        csc.indptr.len(),
        csc.n_cols + 1
    );

    Ok(ExpressionMatrix::SparseCSC(csc))
}

/// Extract integer vector from RObject
fn extract_integer_vector(robj: &RObject) -> Result<Vec<i32>, SeuratError> {
    match robj {
        RObject::IntegerVector(iv) => Ok(iv.data.clone()),
        _ => Err(SeuratError::ParseError(format!(
            "Expected integer vector, got {}",
            robj.type_name()
        ))),
    }
}

/// Extract real vector from RObject
fn extract_real_vector(robj: &RObject) -> Result<Vec<f64>, SeuratError> {
    match robj {
        RObject::DoubleVector(dv) => Ok(dv.data.clone()),
        _ => Err(SeuratError::ParseError(format!(
            "Expected real vector, got {}",
            robj.type_name()
        ))),
    }
}

/// Extract string vector from RObject
fn extract_string_vector(robj: &RObject) -> Result<Vec<String>, SeuratError> {
    match robj {
        RObject::StringVector(sv) => Ok(sv.data.clone()),
        _ => Err(SeuratError::ParseError(format!(
            "Expected string vector, got {}",
            robj.type_name()
        ))),
    }
}

/// Parse metadata (data.frame)
fn parse_metadata(
    robj: &RObject,
    expected_rows: usize,
    file: &RdsFile,
) -> Result<DataFrame, SeuratError> {
    // data.frame is a named list with row.names attribute
    let columns = get_fields_from_object(robj, file).unwrap_or_default();

    // If no columns, return empty DataFrame
    if columns.is_empty() {
        return Ok(DataFrame::empty(expected_rows));
    }

    // Convert R columns to Arrow arrays
    use arrow::array::ArrayRef;

    let mut col_names: Vec<String> = Vec::new();
    let mut col_arrays: Vec<ArrayRef> = Vec::new();

    for (col_name, col_value) in columns {
        // Skip internal columns like _row_id
        if col_name.starts_with('_') {
            continue;
        }

        let array = robject_to_arrow_array(col_value, col_name, expected_rows)?;
        col_names.push(col_name.to_string());
        col_arrays.push(array);
    }

    // Build DataFrame from Arrow arrays
    if col_names.is_empty() {
        return Ok(DataFrame::empty(expected_rows));
    }

    let df = DataFrame::new(col_names, col_arrays, expected_rows)
        .map_err(|e| SeuratError::ParseError(format!("Failed to create DataFrame: {}", e)))?;

    Ok(df)
}

/// Convert RObject to Arrow array
fn robject_to_arrow_array(
    robj: &RObject,
    col_name: &str,
    expected_len: usize,
) -> Result<Arc<dyn arrow::array::Array>, SeuratError> {
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
            eprintln!(
                "Warning: Unsupported RObject type for column {}: {:?}, using placeholder",
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
) -> Result<Arc<dyn arrow::array::Array>, SeuratError> {
    use arrow::array::*;

    // Get the integer codes (1-based in R)
    let codes = match robj {
        RObject::IntegerVector(iv) => &iv.data,
        _ => {
            return Err(SeuratError::ParseError(format!(
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
        .map_err(|e| {
            SeuratError::ParseError(format!("Failed to create dictionary array: {}", e))
        })?;

    Ok(Arc::new(dict_array))
}

/// Parse reductions (dimensionality reductions)
fn parse_reductions(
    robj: &RObject,
    expected_rows: usize,
    file: &RdsFile,
) -> Result<HashMap<String, Embedding>, SeuratError> {
    // Handle empty or null cases
    if matches!(robj, RObject::Null) {
        return Ok(HashMap::new());
    }

    let reduction_list = match get_fields_from_object(robj, file) {
        Ok(list) if list.is_empty() => return Ok(HashMap::new()),
        Ok(list) => list,
        Err(_) => return Ok(HashMap::new()),
    };

    let mut embeddings = HashMap::new();

    for (name, reduction) in reduction_list {
        let embedding = parse_single_reduction(reduction, expected_rows, file)?;
        embeddings.insert(name.to_string(), embedding);
    }

    Ok(embeddings)
}

/// Parse single reduction
fn parse_single_reduction(
    robj: &RObject,
    _expected_rows: usize,
    file: &RdsFile,
) -> Result<Embedding, SeuratError> {
    let reduction_fields = get_fields_from_object(robj, file)?;

    let mut cell_embeddings = None;
    let mut reduction_name = "unknown".to_string();

    for (name, value) in &reduction_fields {
        match *name {
            "cell_embeddings" | "cell.embeddings" => {
                cell_embeddings = Some(*value);
            }
            "key" => {
                if let Ok(key) = extract_string_scalar(value) {
                    reduction_name = key;
                }
            }
            _ => {}
        }
    }

    let cell_embeddings = cell_embeddings.ok_or_else(|| {
        SeuratError::ParseError("Missing cell_embeddings in reduction".to_string())
    })?;

    // Parse matrix
    let matrix = parse_dense_matrix(cell_embeddings)?;

    // Flatten matrix to Vec<f64> (row-major order)
    let n_rows = matrix.len();
    let n_cols = if n_rows > 0 { matrix[0].len() } else { 0 };
    let mut data = Vec::with_capacity(n_rows * n_cols);
    for row in matrix {
        data.extend(row);
    }

    // Create Embedding
    let embedding = Embedding::new(reduction_name, data, n_rows, n_cols)
        .map_err(|e| SeuratError::ValidationError(e))?;

    Ok(embedding)
}

/// Extract feature.loadings from all reductions → gene_loadings (varm)
///
/// For each DimReduc S4 object, extract the feature.loadings matrix.
/// Returns a HashMap mapping reduction key (e.g., "PCs") to Embedding.
fn parse_reduction_loadings(robj: &RObject, file: &RdsFile) -> HashMap<String, Embedding> {
    let mut loadings = HashMap::new();

    if matches!(robj, RObject::Null) {
        return loadings;
    }

    let reduction_list = match get_fields_from_object(robj, file) {
        Ok(list) if !list.is_empty() => list,
        _ => return loadings,
    };

    for (name, reduction) in reduction_list {
        if let Ok(fields) = get_fields_from_object(reduction, file) {
            let mut feature_loadings = None;
            let mut key = name.to_string();

            for (fname, fvalue) in &fields {
                match *fname {
                    "feature.loadings" | "feature_loadings" => {
                        feature_loadings = Some(*fvalue);
                    }
                    "key" => {
                        if let Ok(k) = extract_string_scalar(fvalue) {
                            key = k;
                        }
                    }
                    _ => {}
                }
            }

            if let Some(fl) = feature_loadings {
                if let Ok(matrix) = parse_dense_matrix(fl) {
                    let n_rows = matrix.len();
                    let n_cols = if n_rows > 0 { matrix[0].len() } else { 0 };
                    if n_rows > 0 && n_cols > 0 {
                        let mut data = Vec::with_capacity(n_rows * n_cols);
                        for row in matrix {
                            data.extend(row);
                        }
                        // Use "PCs" as key for PCA loadings (AnnData convention)
                        let varm_key =
                            if name == "pca" || key.starts_with("pca") || key.starts_with("PC") {
                                "PCs".to_string()
                            } else {
                                name.to_string()
                            };
                        if let Ok(emb) = Embedding::new(varm_key.clone(), data, n_rows, n_cols) {
                            loadings.insert(varm_key, emb);
                        }
                    }
                }
            }
        }
    }

    loadings
}

/// Parse graphs slot → cell_pairwise (HashMap<String, PairwiseMatrix>)
///
/// Each graph in Seurat is a dgCMatrix (sparse matrix) with class "Graph".
fn parse_graphs(
    robj: &RObject,
    expected_cells: usize,
    file: &RdsFile,
) -> Result<Option<HashMap<String, PairwiseMatrix>>, SeuratError> {
    if matches!(robj, RObject::Null) {
        return Ok(None);
    }

    let graph_list = match get_fields_from_object(robj, file) {
        Ok(list) if !list.is_empty() => list,
        _ => return Ok(None),
    };

    let mut pairwise = HashMap::new();

    for (name, graph_obj) in graph_list {
        // Try to parse as a sparse matrix (dgCMatrix)
        if let Ok(matrix) = parse_sparse_matrix_from_s4(graph_obj) {
            let (n_rows, n_cols) = matrix.shape();
            if n_rows == expected_cells && n_cols == expected_cells {
                if let Ok(pw) = PairwiseMatrix::new(name.to_string(), matrix) {
                    pairwise.insert(name.to_string(), pw);
                }
            }
        }
    }

    if pairwise.is_empty() {
        Ok(None)
    } else {
        Ok(Some(pairwise))
    }
}

/// Parse a sparse matrix from an S4 object (dgCMatrix format)
fn parse_sparse_matrix_from_s4(robj: &RObject) -> Result<ExpressionMatrix, SeuratError> {
    use crate::ir::SparseMatrixCSC;

    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => {
            return Err(SeuratError::ParseError(
                "Expected S4 for sparse matrix".to_string(),
            ))
        }
    };

    // Extract i (row indices), p (column pointers), x (values), Dim
    let i = s4
        .attributes
        .get("i")
        .ok_or_else(|| SeuratError::ParseError("Missing 'i' in dgCMatrix".to_string()))?;
    let p = s4
        .attributes
        .get("p")
        .ok_or_else(|| SeuratError::ParseError("Missing 'p' in dgCMatrix".to_string()))?;
    let x = s4
        .attributes
        .get("x")
        .ok_or_else(|| SeuratError::ParseError("Missing 'x' in dgCMatrix".to_string()))?;
    let dim = s4
        .attributes
        .get("Dim")
        .ok_or_else(|| SeuratError::ParseError("Missing 'Dim' in dgCMatrix".to_string()))?;

    let dim_vec = extract_integer_vector(dim)?;
    if dim_vec.len() < 2 {
        return Err(SeuratError::ParseError(
            "Dim must have 2 elements".to_string(),
        ));
    }
    let n_rows = dim_vec[0] as usize;
    let n_cols = dim_vec[1] as usize;

    let indices: Vec<usize> = extract_integer_vector(i)?
        .into_iter()
        .map(|v| v as usize)
        .collect();
    let indptr: Vec<usize> = extract_integer_vector(p)?
        .into_iter()
        .map(|v| v as usize)
        .collect();
    let values: Vec<f64> = extract_real_vector(x)?;

    let csc = SparseMatrixCSC::new(values, indices, indptr, n_rows, n_cols)
        .map_err(|e| SeuratError::ParseError(format!("Invalid CSC matrix: {}", e)))?;

    Ok(ExpressionMatrix::SparseCSC(csc))
}

/// Parse dense matrix (R matrix)
fn parse_dense_matrix(robj: &RObject) -> Result<Vec<Vec<f64>>, SeuratError> {
    // R matrix is stored as a real vector with dim attribute
    let (data, dim) = match robj {
        RObject::DoubleVector(dv) => {
            // Check for dim attribute
            if let Some(dim_obj) = dv.attributes.get("dim") {
                let dim = extract_integer_vector(dim_obj)?;
                (dv.data.clone(), dim)
            } else {
                // No dim attribute, assume it's a vector (n x 1 matrix)
                return Ok(dv.data.iter().map(|&x| vec![x]).collect());
            }
        }
        _ => {
            return Err(SeuratError::ParseError(format!(
                "Expected real vector or matrix, got {}",
                robj.type_name()
            )))
        }
    };

    if dim.len() != 2 {
        return Err(SeuratError::ParseError(format!(
            "Expected 2D dimensions, got {}",
            dim.len()
        )));
    }

    let n_rows = dim[0] as usize;
    let n_cols = dim[1] as usize;

    if data.len() != n_rows * n_cols {
        return Err(SeuratError::ParseError(format!(
            "Matrix data length {} doesn't match dimensions {} x {}",
            data.len(),
            n_rows,
            n_cols
        )));
    }

    // R stores matrices in column-major order
    let mut matrix = vec![vec![0.0; n_cols]; n_rows];
    for col in 0..n_cols {
        for row in 0..n_rows {
            matrix[row][col] = data[col * n_rows + row];
        }
    }

    Ok(matrix)
}

/// Parse spatial data from Seurat images slot
fn parse_spatial_data(
    robj: &RObject,
    expected_cells: usize,
    file: &RdsFile,
) -> Result<SpatialData, SeuratError> {
    let image_list = get_fields_from_object(robj, file)?;

    if image_list.is_empty() {
        return Err(SeuratError::ParseError(
            "No images found in Seurat object".to_string(),
        ));
    }

    // Extract first image
    let (_image_name, image_obj) = &image_list[0];
    parse_single_image(image_obj, expected_cells, file)
}

/// Parse single Seurat image object
fn parse_single_image(
    robj: &RObject,
    _expected_cells: usize,
    file: &RdsFile,
) -> Result<SpatialData, SeuratError> {
    let image_fields = get_fields_from_object(robj, file)?;

    let mut coordinates = None;
    let mut scale_factors = None;
    let mut image_data = None;

    for (name, value) in &image_fields {
        match *name {
            "coordinates" => {
                coordinates = Some(*value);
            }
            "scale.factors" => {
                scale_factors = Some(*value);
            }
            "image" => {
                image_data = Some(*value);
            }
            _ => {}
        }
    }

    // Parse coordinates
    let coords = coordinates
        .ok_or_else(|| SeuratError::ParseError("Missing coordinates in image".to_string()))?;

    let coord_matrix = parse_dense_matrix(coords)?;
    let n_cells = coord_matrix.len();
    let n_dims = if n_cells > 0 {
        coord_matrix[0].len()
    } else {
        2
    };

    // Flatten coordinates
    let mut flat_coords = Vec::with_capacity(n_cells * n_dims);
    for row in coord_matrix {
        flat_coords.extend(row);
    }

    // Parse scale factors if present
    let mut parsed_scale_factors = None;
    if let Some(sf) = scale_factors {
        if let Ok(sf_fields) = get_fields_from_object(sf, file) {
            let mut factors = HashMap::new();
            for (name, value) in sf_fields {
                if let RObject::DoubleVector(dv) = value {
                    if !dv.data.is_empty() {
                        factors.insert(name.to_string(), dv.data[0]);
                    }
                }
            }
            if !factors.is_empty() {
                parsed_scale_factors = Some(factors);
            }
        }
    }

    // Create SpatialData
    let spatial = SpatialData::new(flat_coords, n_dims, None, parsed_scale_factors)
        .map_err(|e| SeuratError::ValidationError(e))?;

    // Parse image data if present
    if let Some(_img) = image_data {
        // TODO: Implement image parsing
        // For now, skip image data
    }

    Ok(spatial)
}
