//! IR → R object converters
//!
//! Converts CrossCell IR types to R objects (Seurat, SCE, dgCMatrix, etc.)

use extendr_api::prelude::*;

use crosscell::ir::unstructured::UnstructuredValue;
use crosscell::ir::{
    DataFrame, DenseMatrix, Embedding, ExpressionMatrix, SingleCellData, SparseMatrixCSC,
    SparseMatrixCSR,
};

use arrow::array::{
    Array, ArrayRef, BooleanArray, DictionaryArray, Float32Array, Float64Array, Int32Array,
    Int64Array, LargeStringArray, StringArray,
};
use arrow::datatypes::{DataType, Int32Type};

// ============================================================================
// 10.1: ExpressionMatrix → dgCMatrix / R matrix
// ============================================================================

/// Convert ExpressionMatrix to R dgCMatrix or matrix
pub fn expression_to_r(matrix: &ExpressionMatrix) -> Result<Robj> {
    match matrix {
        ExpressionMatrix::SparseCSR(m) => csr_to_dgcmatrix(m),
        ExpressionMatrix::SparseCSC(m) => csc_to_dgcmatrix(m),
        ExpressionMatrix::Dense(m) => dense_to_matrix(m),
        ExpressionMatrix::Lazy(_) => Err("Lazy matrices must be loaded before conversion".into()),
    }
}

/// Convert CSR sparse matrix to R dgCMatrix (CSC format)
fn csr_to_dgcmatrix(m: &SparseMatrixCSR) -> Result<Robj> {
    // dgCMatrix is CSC format, so we need to convert CSR to CSC
    let (n_rows, n_cols) = (m.n_rows, m.n_cols);
    let nnz = m.data.len();

    // Count elements per column
    let mut col_counts = vec![0usize; n_cols];
    for &col in &m.indices {
        col_counts[col] += 1;
    }

    // Build column pointers
    let mut col_ptr = vec![0i32; n_cols + 1];
    for i in 0..n_cols {
        col_ptr[i + 1] = col_ptr[i] + col_counts[i] as i32;
    }

    // Build row indices and data in CSC order
    let mut row_indices = vec![0i32; nnz];
    let mut csc_data = vec![0.0f64; nnz];
    let mut col_pos = col_ptr[..n_cols].to_vec();

    for row in 0..n_rows {
        let start = m.indptr[row];
        let end = m.indptr[row + 1];
        for idx in start..end {
            let col = m.indices[idx];
            let pos = col_pos[col] as usize;
            row_indices[pos] = row as i32;
            csc_data[pos] = m.data[idx];
            col_pos[col] += 1;
        }
    }

    create_dgcmatrix(&csc_data, &row_indices, &col_ptr, n_rows, n_cols)
}

/// Convert CSC sparse matrix to R dgCMatrix
fn csc_to_dgcmatrix(m: &SparseMatrixCSC) -> Result<Robj> {
    let col_ptr: Vec<i32> = m.indptr.iter().map(|&x| x as i32).collect();
    let row_indices: Vec<i32> = m.indices.iter().map(|&x| x as i32).collect();
    create_dgcmatrix(&m.data, &row_indices, &col_ptr, m.n_rows, m.n_cols)
}

/// Create dgCMatrix from CSC components
fn create_dgcmatrix(
    data: &[f64],
    row_indices: &[i32],
    col_ptr: &[i32],
    n_rows: usize,
    n_cols: usize,
) -> Result<Robj> {
    let r_i = row_indices.to_vec();
    let r_p = col_ptr.to_vec();
    let r_x = data.to_vec();
    let r_nrow = n_rows as i32;
    let r_ncol = n_cols as i32;

    // Ensure Matrix package is loaded and create dgCMatrix using triplet format
    // Convert CSC (i, p, x) to triplet (i, j, x) format
    R!("
        if (!requireNamespace('Matrix', quietly = TRUE)) {
            stop('Matrix package is required')
        }
        p <- {{r_p}}
        j <- rep(seq_along(diff(p)), diff(p))
        Matrix::sparseMatrix(
            i = {{r_i}} + 1L,
            j = j,
            x = {{r_x}},
            dims = c({{r_nrow}}, {{r_ncol}}),
            repr = 'C'
        )
    ")
}

/// Convert dense matrix to R matrix
fn dense_to_matrix(m: &DenseMatrix) -> Result<Robj> {
    // R matrices are column-major, our DenseMatrix is row-major
    let mut col_major = vec![0.0f64; m.n_rows * m.n_cols];
    for row in 0..m.n_rows {
        for col in 0..m.n_cols {
            col_major[col * m.n_rows + row] = m.data[row * m.n_cols + col];
        }
    }

    let r_data = col_major;
    let r_nrow = m.n_rows as i32;
    let r_ncol = m.n_cols as i32;

    R!("matrix({{r_data}}, nrow = {{r_nrow}}, ncol = {{r_ncol}})")
}

// ============================================================================
// 10.3: DataFrame_IR → R data.frame
// ============================================================================

/// Convert DataFrame to R data.frame
pub fn dataframe_to_r(df: &DataFrame) -> Result<Robj> {
    let mut columns: Vec<(&str, Robj)> = Vec::new();

    for (col_name, col_data) in df.columns.iter().zip(df.data.iter()) {
        let r_col = arrow_array_to_r(col_data)?;
        columns.push((col_name.as_str(), r_col));
    }

    let list = List::from_pairs(columns);
    let n_rows = df.n_rows as i32;

    R!("
        df <- as.data.frame({{list}}, stringsAsFactors = FALSE)
        if (nrow(df) == 0 && {{n_rows}} > 0) {
            df <- data.frame(row.names = seq_len({{n_rows}}))
        }
        df
    ")
}

/// Convert Arrow array to R vector
fn arrow_array_to_r(array: &ArrayRef) -> Result<Robj> {
    match array.data_type() {
        DataType::Int64 => {
            let arr = array.as_any().downcast_ref::<Int64Array>().unwrap();
            let values: Vec<i32> = (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        i32::MIN
                    } else {
                        arr.value(i) as i32
                    }
                })
                .collect();
            Ok(values.into_robj())
        }
        DataType::Int32 => {
            let arr = array.as_any().downcast_ref::<Int32Array>().unwrap();
            let values: Vec<i32> = (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        i32::MIN
                    } else {
                        arr.value(i)
                    }
                })
                .collect();
            Ok(values.into_robj())
        }
        DataType::Float64 => {
            let arr = array.as_any().downcast_ref::<Float64Array>().unwrap();
            let values: Vec<f64> = arr.values().iter().copied().collect();
            Ok(values.into_robj())
        }
        DataType::Float32 => {
            let arr = array.as_any().downcast_ref::<Float32Array>().unwrap();
            let values: Vec<f64> = arr.values().iter().map(|&v| v as f64).collect();
            Ok(values.into_robj())
        }
        DataType::Boolean => {
            let arr = array.as_any().downcast_ref::<BooleanArray>().unwrap();
            let values: Vec<Rbool> = (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        Rbool::na_value()
                    } else if arr.value(i) {
                        Rbool::true_value()
                    } else {
                        Rbool::false_value()
                    }
                })
                .collect();
            Ok(values.into_robj())
        }
        DataType::Utf8 => {
            let arr = array.as_any().downcast_ref::<StringArray>().unwrap();
            let values: Vec<Option<&str>> = (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i))
                    }
                })
                .collect();
            Ok(values.into_robj())
        }
        DataType::LargeUtf8 => {
            let arr = array.as_any().downcast_ref::<LargeStringArray>().unwrap();
            let values: Vec<Option<&str>> = (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i))
                    }
                })
                .collect();
            Ok(values.into_robj())
        }
        DataType::Dictionary(_, _) => arrow_dict_to_factor(array),
        _ => {
            let values: Vec<String> = (0..array.len()).map(|_| "NA".to_string()).collect();
            Ok(values.into_robj())
        }
    }
}

/// Convert Arrow DictionaryArray to R factor
fn arrow_dict_to_factor(array: &ArrayRef) -> Result<Robj> {
    if let Some(dict) = array.as_any().downcast_ref::<DictionaryArray<Int32Type>>() {
        let keys = dict.keys();
        let values = dict.values();
        let values_str = values.as_any().downcast_ref::<StringArray>().unwrap();

        let levels: Vec<&str> = (0..values_str.len()).map(|i| values_str.value(i)).collect();

        let codes: Vec<i32> = (0..keys.len())
            .map(|i| {
                if keys.is_null(i) {
                    i32::MIN
                } else {
                    keys.value(i) + 1
                }
            })
            .collect();

        let r_codes = codes;
        let r_levels = levels;

        R!("
            f <- {{r_codes}}
            attr(f, 'levels') <- {{r_levels}}
            class(f) <- 'factor'
            f
        ")
    } else {
        Err("Unsupported dictionary key type".into())
    }
}

// ============================================================================
// 10.5: Embedding → DimReduc / R matrix
// ============================================================================

/// Convert Embedding to R matrix (for SCE reducedDims)
pub fn embedding_to_matrix(emb: &Embedding) -> Result<Robj> {
    let mut col_major = vec![0.0f64; emb.n_rows * emb.n_cols];
    for row in 0..emb.n_rows {
        for col in 0..emb.n_cols {
            col_major[col * emb.n_rows + row] = emb.data[row * emb.n_cols + col];
        }
    }

    let r_data = col_major;
    let r_nrow = emb.n_rows as i32;
    let r_ncol = emb.n_cols as i32;
    let r_name = emb.name.clone();

    R!("
        m <- matrix({{r_data}}, nrow = {{r_nrow}}, ncol = {{r_ncol}})
        colnames(m) <- paste0({{r_name}}, '_', seq_len(ncol(m)))
        m
    ")
}

/// Convert Embedding to Seurat DimReduc S4 object
pub fn embedding_to_dimreduc(emb: &Embedding, cell_names: &[String]) -> Result<Robj> {
    let mat = embedding_to_matrix(emb)?;
    let r_key = format!("{}_", emb.name.to_lowercase());
    let r_cell_names: Vec<String> = cell_names.to_vec();

    R!("
        m <- {{mat}}
        rownames(m) <- {{r_cell_names}}
        Seurat::CreateDimReducObject(
            embeddings = m,
            key = {{r_key}},
            assay = 'RNA'
        )
    ")
}

// ============================================================================
// 10.7: ir_to_seurat and ir_to_sce
// ============================================================================

/// Convert SingleCellData to Seurat object
pub fn ir_to_seurat(data: &SingleCellData) -> Result<Robj> {
    // IR stores matrix as (n_cells, n_genes), but Seurat needs (n_genes, n_cells)
    // So we need to transpose the matrix
    let counts = expression_to_r(&data.expression)?;

    let cell_names = get_row_names(data.metadata.n_cells);
    let gene_names = get_row_names(data.metadata.n_genes);

    let r_cell_names: Vec<String> = cell_names.clone();
    let r_gene_names: Vec<String> = gene_names.clone();

    // Transpose the matrix: IR is (cells x genes), Seurat needs (genes x cells)
    let mut seurat = R!("
        counts <- Matrix::t({{counts}})
        rownames(counts) <- {{r_gene_names}}
        colnames(counts) <- {{r_cell_names}}
        Seurat::CreateSeuratObject(counts = counts, project = 'CrossCell')
    ")?;

    // Add metadata if present
    if !data.cell_metadata.columns.is_empty() {
        let meta = dataframe_to_r(&data.cell_metadata)?;
        let cell_names_for_meta: Vec<String> = cell_names.clone();
        seurat = R!("
            obj <- {{seurat}}
            meta <- {{meta}}
            rownames(meta) <- {{cell_names_for_meta}}
            obj@meta.data <- cbind(obj@meta.data, meta)
            obj
        ")?;
    }

    // Add embeddings
    if let Some(ref embeddings) = data.embeddings {
        for (name, emb) in embeddings {
            let dimreduc = embedding_to_dimreduc(emb, &cell_names)?;
            let r_name = name.to_lowercase();
            seurat = R!("
                obj <- {{seurat}}
                obj@reductions[[{{r_name}}]] <- {{dimreduc}}
                obj
            ")?;
        }
    }

    Ok(seurat)
}

/// Convert SingleCellData to SingleCellExperiment object
pub fn ir_to_sce(data: &SingleCellData) -> Result<Robj> {
    // IR stores matrix as (n_cells, n_genes), but SCE needs (n_genes, n_cells)
    let counts = expression_to_r(&data.expression)?;

    let cell_names = get_row_names(data.metadata.n_cells);
    let gene_names = get_row_names(data.metadata.n_genes);

    let col_data = dataframe_to_r(&data.cell_metadata)?;
    let row_data = dataframe_to_r(&data.gene_metadata)?;

    // Clone for R! macro which moves values
    let r_cell_names1: Vec<String> = cell_names.clone();
    let r_gene_names1: Vec<String> = gene_names.clone();
    let r_cell_names2: Vec<String> = cell_names.clone();
    let r_gene_names2: Vec<String> = gene_names.clone();

    // Transpose the matrix: IR is (cells x genes), SCE needs (genes x cells)
    let mut sce = R!("
        counts <- Matrix::t({{counts}})
        rownames(counts) <- {{r_gene_names1}}
        colnames(counts) <- {{r_cell_names1}}
        
        col_data <- {{col_data}}
        if (nrow(col_data) > 0) rownames(col_data) <- {{r_cell_names2}}
        
        row_data <- {{row_data}}
        if (nrow(row_data) > 0) rownames(row_data) <- {{r_gene_names2}}
        
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts),
            colData = col_data,
            rowData = row_data
        )
    ")?;

    // Add reducedDims
    if let Some(ref embeddings) = data.embeddings {
        for (name, emb) in embeddings {
            let mat = embedding_to_matrix(emb)?;
            let r_name = name.clone();
            let cell_names_for_emb: Vec<String> = cell_names.clone();
            sce = R!("
                obj <- {{sce}}
                m <- {{mat}}
                rownames(m) <- {{cell_names_for_emb}}
                SingleCellExperiment::reducedDim(obj, {{r_name}}) <- m
                obj
            ")?;
        }
    }

    Ok(sce)
}

/// Generate default row names
fn get_row_names(n: usize) -> Vec<String> {
    (1..=n).map(|i| format!("{}", i)).collect()
}

// ============================================================================
// Unstructured data conversion
// ============================================================================

/// Convert UnstructuredValue to R object
pub fn unstructured_to_r(value: &UnstructuredValue) -> Result<Robj> {
    match value {
        UnstructuredValue::Null => Ok(().into_robj()),
        UnstructuredValue::Boolean(b) => Ok((*b).into_robj()),
        UnstructuredValue::Integer(i) => Ok((*i as i32).into_robj()),
        UnstructuredValue::Float(f) => Ok((*f).into_robj()),
        UnstructuredValue::String(s) => Ok(s.as_str().into_robj()),
        UnstructuredValue::Array(arr) => {
            let items: Result<Vec<Robj>> = arr.iter().map(|v| unstructured_to_r(v)).collect();
            Ok(List::from_values(items?).into_robj())
        }
        UnstructuredValue::Dict(map) => {
            let pairs: Result<Vec<(&str, Robj)>> = map
                .iter()
                .map(|(k, v)| Ok((k.as_str(), unstructured_to_r(v)?)))
                .collect();
            Ok(List::from_pairs(pairs?).into_robj())
        }
    }
}
