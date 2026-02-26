//! IR to Seurat conversion
//!
//! This module converts CrossCell IR to native Seurat S4 objects
//! that can be read directly by R's readRDS() + SeuratObject package.

use crate::ir::{
    DataFrame, Embedding, ExpressionMatrix, PairwiseMatrix, SingleCellData, SparseMatrixCSC,
};
use crate::rds::StringEncoding;
use crate::rds::{
    Attributes, DoubleVector, GenericVector, IntegerVector, LogicalVector, RObject, RdsFile,
    S4Object, StringVector,
};
use crate::seurat::error::SeuratError;
use arrow::array::Array as ArrowArray;
use arrow::datatypes::Int32Type;
use std::collections::HashMap;

// ─── helpers ────────────────────────────────────────────────────────────────

/// Create a string scalar RObject
fn create_string_scalar(s: &str) -> RObject {
    let mut sv = StringVector::default();
    sv.add(s.to_string(), StringEncoding::Utf8);
    RObject::StringVector(sv)
}

/// Create an integer vector RObject
fn create_integer_vector(data: Vec<i32>) -> RObject {
    RObject::IntegerVector(IntegerVector {
        data,
        attributes: Attributes::new(),
    })
}

/// Create a real vector RObject
fn create_real_vector(data: Vec<f64>) -> RObject {
    RObject::DoubleVector(DoubleVector {
        data,
        attributes: Attributes::new(),
    })
}

/// Create a string vector RObject
fn create_string_vector(data: Vec<String>) -> RObject {
    let mut sv = StringVector::default();
    for s in data {
        sv.add(s, StringEncoding::Utf8);
    }
    RObject::StringVector(sv)
}

/// Create a named list RObject
fn create_named_list(items: Vec<(String, RObject)>) -> RObject {
    let mut data = Vec::new();
    let mut names = Vec::new();
    for (name, value) in items {
        names.push(name);
        data.push(value);
    }
    let mut attrs = Attributes::new();
    let mut names_vec = StringVector::default();
    for name in names {
        names_vec.add(name, StringEncoding::Utf8);
    }
    attrs.add(
        "names".to_string(),
        RObject::StringVector(names_vec),
        StringEncoding::Utf8,
    );
    RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    })
}

/// Extract names from a DataFrame's `_index` column, falling back to generated names.
fn extract_names_from_df(df: &DataFrame, n: usize, prefix: &str) -> Vec<String> {
    if let Some(array) = df.column("_index") {
        if let Some(str_arr) = array.as_any().downcast_ref::<arrow::array::StringArray>() {
            if str_arr.len() == n {
                return (0..n).map(|i| str_arr.value(i).to_string()).collect();
            }
        }
        if let Some(str_arr) = array
            .as_any()
            .downcast_ref::<arrow::array::LargeStringArray>()
        {
            if str_arr.len() == n {
                return (0..n).map(|i| str_arr.value(i).to_string()).collect();
            }
        }
    }
    // Fallback: generate names
    (0..n).map(|i| format!("{}_{}", prefix, i + 1)).collect()
}

/// Create an R `package_version` object (S3 class, list containing integer vector).
fn create_package_version(major: i32, minor: i32, patch: i32) -> RObject {
    let inner = create_integer_vector(vec![major, minor, patch]);
    let mut gv = GenericVector {
        data: vec![inner],
        attributes: Attributes::new(),
    };
    let mut class_vec = StringVector::default();
    class_vec.add("package_version".to_string(), StringEncoding::Utf8);
    class_vec.add("numeric_version".to_string(), StringEncoding::Utf8);
    gv.attributes.add(
        "class".to_string(),
        RObject::StringVector(class_vec),
        StringEncoding::Utf8,
    );
    RObject::GenericVector(gv)
}

/// Create a LogMap S4 object (logical matrix with row/col names).
/// LogMap inherits from matrix. In R's RDS format, S4 objects that extend
/// basic types are serialized as the base type with class attribute (not as S4 SEXP).
/// So we create a logical vector with dim, dimnames, and class attributes.
fn create_logmap(row_names: &[String], col_names: &[String]) -> RObject {
    let n_rows = row_names.len();
    let n_cols = col_names.len();
    // All TRUE — every cell/feature is present in every layer
    let data: Vec<i32> = vec![1; n_rows * n_cols]; // R logical TRUE = 1

    let mut attrs = Attributes::new();
    // dim
    attrs.add(
        "dim".to_string(),
        create_integer_vector(vec![n_rows as i32, n_cols as i32]),
        StringEncoding::Utf8,
    );
    // dimnames (list of two character vectors)
    let dimnames = RObject::GenericVector(GenericVector {
        data: vec![
            create_string_vector(row_names.to_vec()),
            create_string_vector(col_names.to_vec()),
        ],
        attributes: Attributes::new(),
    });
    attrs.add("dimnames".to_string(), dimnames, StringEncoding::Utf8);
    // class attribute with package
    let mut class_vec = StringVector::default();
    class_vec.add("LogMap".to_string(), StringEncoding::Utf8);
    let mut class_attrs = Attributes::new();
    let mut pkg_vec = StringVector::default();
    pkg_vec.add("SeuratObject".to_string(), StringEncoding::Utf8);
    class_attrs.add(
        "package".to_string(),
        RObject::StringVector(pkg_vec),
        StringEncoding::Utf8,
    );
    class_vec.attributes = class_attrs;
    attrs.add(
        "class".to_string(),
        RObject::StringVector(class_vec),
        StringEncoding::Utf8,
    );

    RObject::LogicalVector(LogicalVector {
        data,
        attributes: attrs,
    })
}

// ─── main entry point ───────────────────────────────────────────────────────

/// Convert IR to a native Seurat S4 RObject that R can read directly.
pub fn ir_to_seurat_rds(ir: &SingleCellData) -> Result<(RObject, RdsFile), SeuratError> {
    let file = RdsFile::default();

    let active_assay = ir
        .metadata
        .active_assay
        .clone()
        .unwrap_or_else(|| "RNA".to_string());
    let n_cells = ir.metadata.n_cells;
    let n_genes = ir.metadata.n_genes;

    let cell_names = extract_names_from_df(&ir.cell_metadata, n_cells, "Cell");
    let gene_names = extract_names_from_df(&ir.gene_metadata, n_genes, "Gene");

    // Build Seurat S4 object
    let mut s4 = S4Object::default();
    s4.class_name = "Seurat".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "SeuratObject".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // --- assays slot (named list of Assay5 S4 objects) ---
    let assays = construct_assays(ir, &active_assay, &cell_names, &gene_names)?;
    s4.attributes
        .add("assays".to_string(), assays, StringEncoding::Utf8);

    // --- meta.data slot (data.frame) ---
    let meta_data = construct_metadata(&ir.cell_metadata, &cell_names)?;
    s4.attributes
        .add("meta.data".to_string(), meta_data, StringEncoding::Utf8);

    // --- active.assay slot ---
    s4.attributes.add(
        "active.assay".to_string(),
        create_string_scalar(&active_assay),
        StringEncoding::Utf8,
    );

    // --- active.ident slot (factor with cell names) ---
    let ident = create_active_ident(&cell_names, "CrossCell");
    s4.attributes
        .add("active.ident".to_string(), ident, StringEncoding::Utf8);

    // --- graphs slot (from cell_pairwise / obsp) ---
    let graphs = if let Some(ref cell_pairwise) = ir.cell_pairwise {
        construct_graphs(cell_pairwise, &cell_names)?
    } else {
        create_named_list(vec![])
    };
    s4.attributes
        .add("graphs".to_string(), graphs, StringEncoding::Utf8);

    // --- neighbors slot (empty named list) ---
    s4.attributes.add(
        "neighbors".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    // --- reductions slot ---
    let reductions = if let Some(embeddings) = &ir.embeddings {
        construct_reductions(
            embeddings,
            &ir.gene_loadings,
            &cell_names,
            &gene_names,
            &active_assay,
        )?
    } else {
        create_named_list(vec![])
    };
    s4.attributes
        .add("reductions".to_string(), reductions, StringEncoding::Utf8);

    // --- images slot ---
    let images = if let Some(spatial) = &ir.spatial {
        construct_spatial_data(spatial, n_cells)?
    } else {
        create_named_list(vec![])
    };
    s4.attributes
        .add("images".to_string(), images, StringEncoding::Utf8);

    // --- project.name slot ---
    s4.attributes.add(
        "project.name".to_string(),
        create_string_scalar("CrossCell"),
        StringEncoding::Utf8,
    );

    // --- misc slot (empty named list) ---
    s4.attributes.add(
        "misc".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    // --- version slot (package_version "5.0.0") ---
    s4.attributes.add(
        "version".to_string(),
        create_package_version(5, 0, 0),
        StringEncoding::Utf8,
    );

    // --- commands slot (empty named list) ---
    s4.attributes.add(
        "commands".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    // --- tools slot (empty named list) ---
    s4.attributes.add(
        "tools".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    Ok((RObject::S4Object(s4), file))
}

/// Create active.ident factor (named factor with one level)
fn create_active_ident(cell_names: &[String], project_name: &str) -> RObject {
    let n = cell_names.len();
    // All cells belong to the same identity class
    let codes: Vec<i32> = vec![1; n]; // 1-based index into levels

    let mut int_vec = IntegerVector {
        data: codes,
        attributes: Attributes::new(),
    };
    // levels
    int_vec.attributes.add(
        "levels".to_string(),
        create_string_vector(vec![project_name.to_string()]),
        StringEncoding::Utf8,
    );
    // class = "factor"
    let mut class_vec = StringVector::default();
    class_vec.add("factor".to_string(), StringEncoding::Utf8);
    int_vec.attributes.add(
        "class".to_string(),
        RObject::StringVector(class_vec),
        StringEncoding::Utf8,
    );
    // names = cell names
    int_vec.attributes.add(
        "names".to_string(),
        create_string_vector(cell_names.to_vec()),
        StringEncoding::Utf8,
    );

    RObject::IntegerVector(int_vec)
}

// ─── assay construction ─────────────────────────────────────────────────────

/// Construct assays slot (named list of Assay5 S4 objects)
fn construct_assays(
    ir: &SingleCellData,
    active_assay: &str,
    cell_names: &[String],
    gene_names: &[String],
) -> Result<RObject, SeuratError> {
    let mut assay_list = Vec::new();

    // Main assay
    let main_assay = construct_single_assay(
        &ir.expression,
        &ir.gene_metadata,
        cell_names,
        gene_names,
        active_assay,
    )?;
    assay_list.push((active_assay.to_string(), main_assay));

    // Additional assays from layers
    if let Some(layers) = &ir.layers {
        for (name, layer_expr) in layers {
            let layer_assay = construct_single_assay(
                layer_expr,
                &ir.gene_metadata,
                cell_names,
                gene_names,
                name,
            )?;
            assay_list.push((name.clone(), layer_assay));
        }
    }

    Ok(create_named_list(assay_list))
}

/// Construct a single Assay5 S4 object
fn construct_single_assay(
    expression: &ExpressionMatrix,
    _gene_metadata: &DataFrame,
    cell_names: &[String],
    gene_names: &[String],
    assay_name: &str,
) -> Result<RObject, SeuratError> {
    let mut s4 = S4Object::default();
    s4.class_name = "Assay5".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "SeuratObject".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // layers slot: named list with "counts" → dgCMatrix
    let dgc = expression_to_dgcmatrix(expression, gene_names, cell_names)?;
    let layers = create_named_list(vec![("counts".to_string(), dgc)]);
    s4.attributes
        .add("layers".to_string(), layers, StringEncoding::Utf8);

    // cells slot: LogMap (cell_names × layer_names)
    let layer_names = vec!["counts".to_string()];
    let cells_logmap = create_logmap(cell_names, &layer_names);
    s4.attributes
        .add("cells".to_string(), cells_logmap, StringEncoding::Utf8);

    // features slot: LogMap (gene_names × layer_names)
    let features_logmap = create_logmap(gene_names, &layer_names);
    s4.attributes.add(
        "features".to_string(),
        features_logmap,
        StringEncoding::Utf8,
    );

    // default slot: integer(1) = 1
    s4.attributes.add(
        "default".to_string(),
        create_integer_vector(vec![1]),
        StringEncoding::Utf8,
    );

    // assay.orig slot: character("")
    s4.attributes.add(
        "assay.orig".to_string(),
        create_string_scalar(""),
        StringEncoding::Utf8,
    );

    // meta.data slot: empty data.frame with n_genes rows
    let assay_meta = construct_empty_dataframe(gene_names.len())?;
    s4.attributes
        .add("meta.data".to_string(), assay_meta, StringEncoding::Utf8);

    // misc slot: empty named list
    s4.attributes.add(
        "misc".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    // key slot: lowercase assay name + "_"
    let key = format!("{}_", assay_name.to_lowercase());
    s4.attributes.add(
        "key".to_string(),
        create_string_scalar(&key),
        StringEncoding::Utf8,
    );

    Ok(RObject::S4Object(s4))
}

/// Construct an empty data.frame with given number of rows
fn construct_empty_dataframe(n_rows: usize) -> Result<RObject, SeuratError> {
    let mut attrs = Attributes::new();
    // names: character(0)
    attrs.add(
        "names".to_string(),
        create_string_vector(vec![]),
        StringEncoding::Utf8,
    );
    // class
    let mut class_vec = StringVector::default();
    class_vec.add("data.frame".to_string(), StringEncoding::Utf8);
    attrs.add(
        "class".to_string(),
        RObject::StringVector(class_vec),
        StringEncoding::Utf8,
    );
    // row.names: compact form c(NA, -n_rows)
    let row_names = create_integer_vector(vec![i32::MIN, -(n_rows as i32)]);
    attrs.add("row.names".to_string(), row_names, StringEncoding::Utf8);

    Ok(RObject::GenericVector(GenericVector {
        data: vec![],
        attributes: attrs,
    }))
}

// ─── sparse matrix conversion ───────────────────────────────────────────────

/// Convert ExpressionMatrix to dgCMatrix (R sparse matrix)
pub fn expression_to_dgcmatrix(
    expression: &ExpressionMatrix,
    gene_names: &[String],
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    match expression {
        ExpressionMatrix::SparseCSC(sparse) => {
            sparse_csc_to_dgcmatrix(sparse, gene_names, cell_names)
        }
        ExpressionMatrix::SparseCSR(sparse) => {
            use crate::sparse::convert::csr_to_csc;
            let csc = csr_to_csc(sparse);
            sparse_csc_to_dgcmatrix(&csc, gene_names, cell_names)
        }
        ExpressionMatrix::Dense(dense) => dense_to_dgcmatrix(dense, gene_names, cell_names),
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                expression_to_dgcmatrix(&cached, gene_names, cell_names)
            } else {
                Err(SeuratError::ConversionError(
                    "Cannot convert LazyMatrix without loading data first.".to_string(),
                ))
            }
        }
    }
}

/// Convert DenseMatrix to dgCMatrix
fn dense_to_dgcmatrix(
    dense: &crate::ir::DenseMatrix,
    gene_names: &[String],
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    let n_rows = dense.n_rows;
    let n_cols = dense.n_cols;
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0usize];
    for col in 0..n_cols {
        for row in 0..n_rows {
            let val = dense.data[row * n_cols + col];
            if val != 0.0 {
                data.push(val);
                indices.push(row);
            }
        }
        indptr.push(data.len());
    }
    let csc = SparseMatrixCSC {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    };
    sparse_csc_to_dgcmatrix(&csc, gene_names, cell_names)
}

/// Convert SparseMatrixCSC to dgCMatrix S4 object
fn sparse_csc_to_dgcmatrix(
    sparse: &SparseMatrixCSC,
    gene_names: &[String],
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    log::debug!(
        "sparse_csc_to_dgcmatrix: IR input: {} cells x {} genes",
        sparse.n_rows,
        sparse.n_cols
    );

    // dgCMatrix in R stores genes × cells (transposed from our cells × genes).
    // Convert our CSC (cells × genes) → CSR, then treat CSR as CSC for (genes × cells).
    use crate::sparse::convert::csc_to_csr;
    let csr = csc_to_csr(sparse);

    let mut s4 = S4Object::default();
    s4.class_name = "dgCMatrix".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "Matrix".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // i: row indices (0-based)
    let i_vec: Vec<i32> = csr.indices.iter().map(|&x| x as i32).collect();
    s4.attributes.add(
        "i".to_string(),
        create_integer_vector(i_vec),
        StringEncoding::Utf8,
    );

    // p: column pointers
    let p_vec: Vec<i32> = csr.indptr.iter().map(|&x| x as i32).collect();
    s4.attributes.add(
        "p".to_string(),
        create_integer_vector(p_vec),
        StringEncoding::Utf8,
    );

    // x: non-zero values
    s4.attributes.add(
        "x".to_string(),
        create_real_vector(csr.data.clone()),
        StringEncoding::Utf8,
    );

    // Dim: [n_genes, n_cells]
    let dim_vec = vec![sparse.n_cols as i32, sparse.n_rows as i32];
    s4.attributes.add(
        "Dim".to_string(),
        create_integer_vector(dim_vec),
        StringEncoding::Utf8,
    );

    // Dimnames: list(gene_names_or_NULL, cell_names_or_NULL)
    let row_dimnames = if gene_names.is_empty() {
        RObject::Null
    } else {
        create_string_vector(gene_names.to_vec())
    };
    let col_dimnames = if cell_names.is_empty() {
        RObject::Null
    } else {
        create_string_vector(cell_names.to_vec())
    };
    let dimnames = RObject::GenericVector(GenericVector {
        data: vec![row_dimnames, col_dimnames],
        attributes: Attributes::new(),
    });
    s4.attributes
        .add("Dimnames".to_string(), dimnames, StringEncoding::Utf8);

    // factors: empty list
    s4.attributes.add(
        "factors".to_string(),
        RObject::GenericVector(GenericVector::default()),
        StringEncoding::Utf8,
    );

    Ok(RObject::S4Object(s4))
}

// ─── metadata construction ──────────────────────────────────────────────────

/// Construct metadata (data.frame) with proper row.names
fn construct_metadata(df: &DataFrame, cell_names: &[String]) -> Result<RObject, SeuratError> {
    if df.n_rows == 0 {
        let mut attrs = Attributes::new();
        attrs.add(
            "names".to_string(),
            create_string_vector(vec![]),
            StringEncoding::Utf8,
        );
        let mut class_vec = StringVector::default();
        class_vec.add("data.frame".to_string(), StringEncoding::Utf8);
        attrs.add(
            "class".to_string(),
            RObject::StringVector(class_vec),
            StringEncoding::Utf8,
        );
        attrs.add(
            "row.names".to_string(),
            create_integer_vector(vec![]),
            StringEncoding::Utf8,
        );
        return Ok(RObject::GenericVector(GenericVector {
            data: vec![],
            attributes: attrs,
        }));
    }

    // Convert DataFrame columns to R vectors, skipping _index
    let mut columns = Vec::new();
    for col_name in &df.columns {
        if col_name == "_index" {
            continue;
        }
        if let Some(array) = df.column(col_name) {
            let r_col = arrow_array_to_robject(array, col_name)?;
            columns.push((col_name.clone(), r_col));
        }
    }

    // If no columns, add a dummy
    if columns.is_empty() {
        let dummy_col: Vec<i32> = (0..df.n_rows as i32).collect();
        columns.push(("_row_id".to_string(), create_integer_vector(dummy_col)));
    }

    let col_names: Vec<String> = columns.iter().map(|(n, _)| n.clone()).collect();
    let data: Vec<RObject> = columns.into_iter().map(|(_, v)| v).collect();

    let mut attrs = Attributes::new();
    // names (column names) — must come before row.names for R data.frame compatibility
    attrs.add(
        "names".to_string(),
        create_string_vector(col_names),
        StringEncoding::Utf8,
    );
    // class
    let mut class_vec = StringVector::default();
    class_vec.add("data.frame".to_string(), StringEncoding::Utf8);
    attrs.add(
        "class".to_string(),
        RObject::StringVector(class_vec),
        StringEncoding::Utf8,
    );
    // row.names = cell names (character vector)
    attrs.add(
        "row.names".to_string(),
        create_string_vector(cell_names.to_vec()),
        StringEncoding::Utf8,
    );

    Ok(RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    }))
}

/// Convert Arrow array to RObject
fn arrow_array_to_robject(
    array: &dyn arrow::array::Array,
    col_name: &str,
) -> Result<RObject, SeuratError> {
    use arrow::array::*;
    use arrow::datatypes::DataType;

    match array.data_type() {
        DataType::Int64 => {
            let arr = array.as_any().downcast_ref::<Int64Array>().ok_or_else(|| {
                SeuratError::ConversionError(format!(
                    "Failed to downcast {} to Int64Array",
                    col_name
                ))
            })?;
            let values: Vec<i32> = arr
                .iter()
                .map(|v| v.map(|x| x as i32).unwrap_or(i32::MIN))
                .collect();
            Ok(create_integer_vector(values))
        }
        DataType::Int32 => {
            let arr = array.as_any().downcast_ref::<Int32Array>().ok_or_else(|| {
                SeuratError::ConversionError(format!(
                    "Failed to downcast {} to Int32Array",
                    col_name
                ))
            })?;
            let values: Vec<i32> = arr.iter().map(|v| v.unwrap_or(i32::MIN)).collect();
            Ok(create_integer_vector(values))
        }
        DataType::Float64 => {
            let arr = array
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| {
                    SeuratError::ConversionError(format!(
                        "Failed to downcast {} to Float64Array",
                        col_name
                    ))
                })?;
            let values: Vec<f64> = arr.iter().map(|v| v.unwrap_or(f64::NAN)).collect();
            Ok(create_real_vector(values))
        }
        DataType::Float32 => {
            let arr = array
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| {
                    SeuratError::ConversionError(format!(
                        "Failed to downcast {} to Float32Array",
                        col_name
                    ))
                })?;
            let values: Vec<f64> = arr
                .iter()
                .map(|v| v.map(|x| x as f64).unwrap_or(f64::NAN))
                .collect();
            Ok(create_real_vector(values))
        }
        DataType::Utf8 => {
            let arr = array
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    SeuratError::ConversionError(format!(
                        "Failed to downcast {} to StringArray",
                        col_name
                    ))
                })?;
            let values: Vec<String> = arr
                .iter()
                .map(|v| v.map(|s| s.to_string()).unwrap_or_else(|| "NA".to_string()))
                .collect();
            Ok(create_string_vector(values))
        }
        DataType::LargeUtf8 => {
            let arr = array
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .ok_or_else(|| {
                    SeuratError::ConversionError(format!(
                        "Failed to downcast {} to LargeStringArray",
                        col_name
                    ))
                })?;
            let values: Vec<String> = arr
                .iter()
                .map(|v| v.map(|s| s.to_string()).unwrap_or_else(|| "NA".to_string()))
                .collect();
            Ok(create_string_vector(values))
        }
        DataType::Boolean => {
            let arr = array
                .as_any()
                .downcast_ref::<BooleanArray>()
                .ok_or_else(|| {
                    SeuratError::ConversionError(format!(
                        "Failed to downcast {} to BooleanArray",
                        col_name
                    ))
                })?;
            let values: Vec<i32> = arr
                .iter()
                .map(|v| match v {
                    Some(true) => 1,
                    Some(false) => 0,
                    None => i32::MIN,
                })
                .collect();
            Ok(create_integer_vector(values))
        }
        DataType::Dictionary(key_type, _value_type) => {
            match key_type.as_ref() {
                DataType::Int32 => {
                    let dict_arr = array
                        .as_any()
                        .downcast_ref::<DictionaryArray<Int32Type>>()
                        .ok_or_else(|| {
                            SeuratError::ConversionError(format!(
                                "Failed to downcast {} to DictionaryArray",
                                col_name
                            ))
                        })?;
                    let keys: Vec<i32> = dict_arr
                        .keys()
                        .iter()
                        .map(|v: Option<i32>| v.map(|k| k + 1).unwrap_or(i32::MIN))
                        .collect();
                    let values = dict_arr.values();
                    let levels: Vec<String> =
                        if let Some(str_arr) = values.as_any().downcast_ref::<StringArray>() {
                            str_arr
                                .iter()
                                .map(|v: Option<&str>| v.unwrap_or("").to_string())
                                .collect()
                        } else {
                            (0..values.len()).map(|i| format!("Level_{}", i)).collect()
                        };
                    let mut int_vec = IntegerVector {
                        data: keys,
                        attributes: Attributes::new(),
                    };
                    int_vec.attributes.add(
                        "levels".to_string(),
                        create_string_vector(levels),
                        StringEncoding::Utf8,
                    );
                    let mut class_vec = StringVector::default();
                    class_vec.add("factor".to_string(), StringEncoding::Utf8);
                    int_vec.attributes.add(
                        "class".to_string(),
                        RObject::StringVector(class_vec),
                        StringEncoding::Utf8,
                    );
                    Ok(RObject::IntegerVector(int_vec))
                }
                _ => {
                    eprintln!("Warning: Unsupported dictionary key type for column {}, converting to string", col_name);
                    let values: Vec<String> =
                        (0..array.len()).map(|i| format!("Value_{}", i)).collect();
                    Ok(create_string_vector(values))
                }
            }
        }
        _ => {
            eprintln!(
                "Warning: Unsupported data type {:?} for column {}, converting to string",
                array.data_type(),
                col_name
            );
            let values: Vec<String> = (0..array.len()).map(|i| format!("Value_{}", i)).collect();
            Ok(create_string_vector(values))
        }
    }
}

// ─── reductions ─────────────────────────────────────────────────────────────

/// Construct reductions slot (named list of DimReduc S4 objects)
fn construct_reductions(
    embeddings: &HashMap<String, Embedding>,
    gene_loadings: &Option<HashMap<String, Embedding>>,
    cell_names: &[String],
    gene_names: &[String],
    assay_name: &str,
) -> Result<RObject, SeuratError> {
    let mut reduction_list = Vec::new();
    for (name, embedding) in embeddings {
        // Find matching gene loadings for this reduction
        // AnnData convention: obsm key "X_pca" → varm key "PCs"
        // We try multiple matching strategies
        let loading = gene_loadings.as_ref().and_then(|gl| {
            // Direct name match (e.g., "PCs" → "PCs")
            gl.get(name)
                // AnnData convention: "X_pca" in obsm → "PCs" in varm
                .or_else(|| {
                    if name == "X_pca" || name == "pca" {
                        gl.get("PCs")
                    } else {
                        None
                    }
                })
        });
        let reduction =
            construct_single_reduction(embedding, loading, cell_names, gene_names, assay_name)?;
        // Strip AnnData "X_" prefix for Seurat reduction names
        let seurat_name = if name.starts_with("X_") || name.starts_with("x_") {
            name[2..].to_string()
        } else {
            name.clone()
        };
        reduction_list.push((seurat_name, reduction));
    }
    Ok(create_named_list(reduction_list))
}

/// Construct a single DimReduc S4 object
fn construct_single_reduction(
    embedding: &Embedding,
    loading: Option<&Embedding>,
    cell_names: &[String],
    gene_names: &[String],
    assay_name: &str,
) -> Result<RObject, SeuratError> {
    let mut s4 = S4Object::default();
    s4.class_name = "DimReduc".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "SeuratObject".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // cell.embeddings slot: matrix with dimnames
    let cell_embeddings = embedding_to_matrix(embedding, cell_names)?;
    s4.attributes.add(
        "cell.embeddings".to_string(),
        cell_embeddings,
        StringEncoding::Utf8,
    );

    // feature.loadings: gene loadings matrix or empty
    let feature_loadings = if let Some(load) = loading {
        loading_to_matrix(load, gene_names, embedding.n_cols)?
    } else {
        RObject::DoubleVector(DoubleVector {
            data: vec![],
            attributes: Attributes::new(),
        })
    };
    s4.attributes.add(
        "feature.loadings".to_string(),
        feature_loadings,
        StringEncoding::Utf8,
    );

    // feature.loadings.projected: empty matrix
    let empty_mat = DoubleVector {
        data: vec![],
        attributes: Attributes::new(),
    };
    s4.attributes.add(
        "feature.loadings.projected".to_string(),
        RObject::DoubleVector(empty_mat),
        StringEncoding::Utf8,
    );

    // assay.used
    s4.attributes.add(
        "assay.used".to_string(),
        create_string_scalar(assay_name),
        StringEncoding::Utf8,
    );

    // global: FALSE
    let lv = LogicalVector {
        data: vec![0],
        attributes: Attributes::new(),
    };
    s4.attributes.add(
        "global".to_string(),
        RObject::LogicalVector(lv),
        StringEncoding::Utf8,
    );

    // stdev: numeric(0)
    s4.attributes.add(
        "stdev".to_string(),
        create_real_vector(vec![]),
        StringEncoding::Utf8,
    );

    // jackstraw: JackStrawData S4 with empty matrix slots
    let jackstraw = create_jackstraw_data();
    s4.attributes
        .add("jackstraw".to_string(), jackstraw, StringEncoding::Utf8);

    // misc: empty list
    s4.attributes.add(
        "misc".to_string(),
        create_named_list(vec![]),
        StringEncoding::Utf8,
    );

    // key: e.g. "pca_", "umap_"
    // Strip AnnData "X_" prefix for Seurat compatibility
    let base_name = embedding.name.to_lowercase();
    let clean_name = if base_name.starts_with("x_") {
        &base_name[2..]
    } else {
        &base_name
    };
    let key = format!("{}_", clean_name);
    s4.attributes.add(
        "key".to_string(),
        create_string_scalar(&key),
        StringEncoding::Utf8,
    );

    Ok(RObject::S4Object(s4))
}

/// Create an empty JackStrawData S4 object
fn create_jackstraw_data() -> RObject {
    let mut s4 = S4Object::default();
    s4.class_name = "JackStrawData".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "SeuratObject".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    let empty_mat = RObject::DoubleVector(DoubleVector {
        data: vec![],
        attributes: Attributes::new(),
    });
    s4.attributes.add(
        "empirical.p.values".to_string(),
        empty_mat.clone(),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "fake.reduction.scores".to_string(),
        empty_mat.clone(),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "empirical.p.values.full".to_string(),
        empty_mat.clone(),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "overall.p.values".to_string(),
        empty_mat,
        StringEncoding::Utf8,
    );

    RObject::S4Object(s4)
}

/// Convert Embedding to R matrix with dimnames
fn embedding_to_matrix(
    embedding: &Embedding,
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    let n_rows = embedding.n_rows;
    let n_cols = embedding.n_cols;

    // Convert row-major → column-major
    let mut col_major_data = Vec::with_capacity(n_rows * n_cols);
    for col in 0..n_cols {
        for row in 0..n_rows {
            col_major_data.push(embedding.data[row * n_cols + col]);
        }
    }

    let mut attrs = Attributes::new();
    // dim
    attrs.add(
        "dim".to_string(),
        create_integer_vector(vec![n_rows as i32, n_cols as i32]),
        StringEncoding::Utf8,
    );
    // dimnames: list(cell_names, component_names)
    // Strip AnnData "X_" prefix for Seurat compatibility
    let base = &embedding.name;
    let clean = if base.starts_with("X_") || base.starts_with("x_") {
        &base[2..]
    } else {
        base.as_str()
    };
    let comp_names: Vec<String> = (1..=n_cols).map(|i| format!("{}_{}", clean, i)).collect();
    let dimnames = RObject::GenericVector(GenericVector {
        data: vec![
            create_string_vector(cell_names.to_vec()),
            create_string_vector(comp_names),
        ],
        attributes: Attributes::new(),
    });
    attrs.add("dimnames".to_string(), dimnames, StringEncoding::Utf8);

    Ok(RObject::DoubleVector(DoubleVector {
        data: col_major_data,
        attributes: attrs,
    }))
}

/// Convert gene loadings (varm) to R matrix for feature.loadings slot
///
/// The loading matrix is n_genes × n_components (row-major in IR).
/// R expects column-major storage.
fn loading_to_matrix(
    loading: &Embedding,
    gene_names: &[String],
    expected_cols: usize,
) -> Result<RObject, SeuratError> {
    let n_rows = loading.n_rows; // n_genes
    let n_cols = std::cmp::min(loading.n_cols, expected_cols); // n_components

    // Convert row-major → column-major
    let mut col_major_data = Vec::with_capacity(n_rows * n_cols);
    for col in 0..n_cols {
        for row in 0..n_rows {
            col_major_data.push(loading.data[row * loading.n_cols + col]);
        }
    }

    let mut attrs = Attributes::new();
    // dim
    attrs.add(
        "dim".to_string(),
        create_integer_vector(vec![n_rows as i32, n_cols as i32]),
        StringEncoding::Utf8,
    );
    // dimnames: list(gene_names, component_names)
    let comp_names: Vec<String> = (1..=n_cols).map(|i| format!("PC_{}", i)).collect();
    let row_names = if gene_names.len() == n_rows {
        gene_names.to_vec()
    } else {
        (1..=n_rows).map(|i| format!("Gene_{}", i)).collect()
    };
    let dimnames = RObject::GenericVector(GenericVector {
        data: vec![
            create_string_vector(row_names),
            create_string_vector(comp_names),
        ],
        attributes: Attributes::new(),
    });
    attrs.add("dimnames".to_string(), dimnames, StringEncoding::Utf8);

    Ok(RObject::DoubleVector(DoubleVector {
        data: col_major_data,
        attributes: attrs,
    }))
}

// ─── graphs (obsp → Seurat@graphs) ─────────────────────────────────────────

/// Construct graphs slot from cell_pairwise matrices
///
/// Each pairwise matrix becomes a Graph S4 object (inherits from dgCMatrix).
/// Seurat stores neighbor graphs (connectivities, distances) in @graphs.
fn construct_graphs(
    cell_pairwise: &HashMap<String, PairwiseMatrix>,
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    let mut graph_list = Vec::new();
    for (name, pw) in cell_pairwise {
        let graph = construct_single_graph(pw, cell_names)?;
        graph_list.push((name.clone(), graph));
    }
    Ok(create_named_list(graph_list))
}

/// Construct a single Graph S4 object from a PairwiseMatrix
///
/// A Seurat Graph is essentially a dgCMatrix with class = c("Graph", "dgCMatrix", ...)
fn construct_single_graph(
    pw: &PairwiseMatrix,
    cell_names: &[String],
) -> Result<RObject, SeuratError> {
    // Convert the pairwise matrix to dgCMatrix format
    // For graphs, both row and col names are cell names (square matrix)
    let dgc = expression_to_dgcmatrix(&pw.matrix, cell_names, cell_names)?;

    // Add dimnames (row = cell_names, col = cell_names)
    match dgc {
        RObject::S4Object(mut s4) => {
            let dimnames = RObject::GenericVector(GenericVector {
                data: vec![
                    create_string_vector(cell_names.to_vec()),
                    create_string_vector(cell_names.to_vec()),
                ],
                attributes: Attributes::new(),
            });
            s4.attributes
                .add("Dimnames".to_string(), dimnames, StringEncoding::Utf8);

            // Override class to Graph (which inherits from dgCMatrix)
            let mut class_vec = StringVector::default();
            class_vec.add("Graph".to_string(), StringEncoding::Utf8);
            class_vec.add("dgCMatrix".to_string(), StringEncoding::Utf8);
            class_vec.add("CsparseMatrix".to_string(), StringEncoding::Utf8);
            class_vec.add("sparseMatrix".to_string(), StringEncoding::Utf8);
            s4.class_name = "Graph".to_string();
            s4.package_name = "SeuratObject".to_string();

            Ok(RObject::S4Object(s4))
        }
        _ => Err(SeuratError::ConversionError(
            "Expected S4 dgCMatrix from expression_to_dgcmatrix".to_string(),
        )),
    }
}

// ─── spatial data ───────────────────────────────────────────────────────────

/// Construct spatial data (images slot)
fn construct_spatial_data(
    spatial: &crate::ir::SpatialData,
    expected_cells: usize,
) -> Result<RObject, SeuratError> {
    if spatial.n_cells() != expected_cells {
        return Err(SeuratError::ValidationError(format!(
            "Spatial data cell count {} doesn't match expected {}",
            spatial.n_cells(),
            expected_cells
        )));
    }
    let image_obj = construct_image_object(spatial)?;
    Ok(create_named_list(vec![("slice1".to_string(), image_obj)]))
}

fn construct_image_object(spatial: &crate::ir::SpatialData) -> Result<RObject, SeuratError> {
    let mut image_fields = Vec::new();
    let coordinates = construct_spatial_coordinates(spatial)?;
    image_fields.push(("coordinates".to_string(), coordinates));
    if let Some(scale_factors) = &spatial.scale_factors {
        image_fields.push((
            "scale.factors".to_string(),
            construct_scale_factors(scale_factors)?,
        ));
    }
    if let Some(images) = &spatial.images {
        if !images.is_empty() {
            image_fields.push(("image".to_string(), construct_image_data(&images[0])?));
        }
    }
    Ok(create_named_list(image_fields))
}

fn construct_spatial_coordinates(spatial: &crate::ir::SpatialData) -> Result<RObject, SeuratError> {
    let n_cells = spatial.n_cells();
    let n_dims = spatial.n_dims;
    let mut col_major_data = Vec::with_capacity(n_cells * n_dims);
    for col in 0..n_dims {
        for row in 0..n_cells {
            col_major_data.push(spatial.coordinates[row * n_dims + col]);
        }
    }
    let mut double_vec = DoubleVector {
        data: col_major_data,
        attributes: Attributes::new(),
    };
    double_vec.attributes.add(
        "dim".to_string(),
        create_integer_vector(vec![n_cells as i32, n_dims as i32]),
        StringEncoding::Utf8,
    );
    Ok(RObject::DoubleVector(double_vec))
}

fn construct_scale_factors(
    scale_factors: &std::collections::HashMap<String, f64>,
) -> Result<RObject, SeuratError> {
    let factors: Vec<(String, RObject)> = scale_factors
        .iter()
        .map(|(name, value)| (name.clone(), create_real_vector(vec![*value])))
        .collect();
    Ok(create_named_list(factors))
}

fn construct_image_data(image: &crate::ir::SpatialImage) -> Result<RObject, SeuratError> {
    let mut double_vec = DoubleVector {
        data: Vec::new(),
        attributes: Attributes::new(),
    };
    let dim_vec = vec![image.height as i32, image.width as i32, 3];
    double_vec.attributes.add(
        "dim".to_string(),
        create_integer_vector(dim_vec),
        StringEncoding::Utf8,
    );
    Ok(RObject::DoubleVector(double_vec))
}

// ─── write convenience ──────────────────────────────────────────────────────

/// Write IR to Seurat RDS file
pub fn write_seurat_rds(ir: &SingleCellData, path: &str) -> Result<(), SeuratError> {
    use crate::rds::write_rds;
    use std::path::Path;

    let (seurat_robj, mut file) = ir_to_seurat_rds(ir)?;
    file.object = seurat_robj;

    write_rds(&file, Path::new(path))
        .map_err(|e| SeuratError::IoError(format!("Failed to write RDS: {}", e)))?;

    Ok(())
}

// ─── tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{DatasetMetadata, ExpressionMatrix, SparseMatrixCSC};

    #[test]
    fn test_ir_to_seurat_minimal() {
        let n_cells = 10;
        let n_genes = 5;
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 2, 4, 6, 8, 10],
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],
            data: vec![1.0; 10],
        };
        let expression = ExpressionMatrix::SparseCSC(sparse);
        let cell_metadata = DataFrame::empty(n_cells);
        let gene_metadata = DataFrame::empty(n_genes);
        let metadata = DatasetMetadata::new(n_cells, n_genes, "test".to_string());

        let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
            .expect("Failed to create IR");

        let result = ir_to_seurat_rds(&ir);
        assert!(
            result.is_ok(),
            "Failed to convert IR to Seurat: {:?}",
            result.err()
        );

        let (seurat, _file) = result.unwrap();
        match &seurat {
            RObject::S4Object(s4) => {
                assert_eq!(s4.class_name, "Seurat");
                assert_eq!(s4.package_name, "SeuratObject");
                // Check required slots exist
                assert!(s4.attributes.get("assays").is_some());
                assert!(s4.attributes.get("meta.data").is_some());
                assert!(s4.attributes.get("active.assay").is_some());
                assert!(s4.attributes.get("active.ident").is_some());
                assert!(s4.attributes.get("version").is_some());
                assert!(s4.attributes.get("project.name").is_some());
            }
            _ => panic!("Expected S4 object, got {:?}", seurat.type_name()),
        }
    }

    #[test]
    fn test_write_seurat_rds() {
        let n_cells = 5;
        let n_genes = 10;
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            indices: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
        };
        let expression = ExpressionMatrix::SparseCSC(sparse);
        let cell_metadata = DataFrame::empty(n_cells);
        let gene_metadata = DataFrame::empty(n_genes);
        let metadata = DatasetMetadata::new(n_cells, n_genes, "test".to_string());

        let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
            .expect("Failed to create IR");

        let output_path = "tests/data/rust_generated_seurat.rds";
        let result = write_seurat_rds(&ir, output_path);
        assert!(
            result.is_ok(),
            "Failed to write Seurat RDS: {:?}",
            result.err()
        );

        use std::path::Path;
        assert!(Path::new(output_path).exists(), "Output file not created");
    }

    #[test]
    fn test_sparse_csc_to_dgcmatrix() {
        let n_cells = 3;
        let n_genes = 5;
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 1, 2, 3, 4, 5],
            indices: vec![0, 1, 2, 0, 1],
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        };
        let gene_names: Vec<String> = (1..=5).map(|i| format!("Gene_{}", i)).collect();
        let cell_names: Vec<String> = (1..=3).map(|i| format!("Cell_{}", i)).collect();

        let result = sparse_csc_to_dgcmatrix(&sparse, &gene_names, &cell_names);
        assert!(
            result.is_ok(),
            "Failed to convert to dgCMatrix: {:?}",
            result.err()
        );

        match result.unwrap() {
            RObject::S4Object(s4) => {
                assert_eq!(s4.class_name, "dgCMatrix");
                assert_eq!(s4.package_name, "Matrix");
                // Check Dim
                if let Some(RObject::IntegerVector(int_vec)) = s4.attributes.get("Dim") {
                    assert_eq!(int_vec.data[0], 5); // genes
                    assert_eq!(int_vec.data[1], 3); // cells
                }
                // Check Dimnames is a list, not NULL
                match s4.attributes.get("Dimnames") {
                    Some(RObject::GenericVector(gv)) => {
                        assert_eq!(gv.data.len(), 2);
                    }
                    other => panic!(
                        "Expected Dimnames to be a list, got {:?}",
                        other.map(|o| o.type_name())
                    ),
                }
            }
            _ => panic!("Expected S4 object"),
        }
    }

    #[test]
    fn test_embedding_to_matrix() {
        let embedding = Embedding::new("pca".to_string(), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 3, 2)
            .expect("Failed to create embedding");
        let cell_names: Vec<String> = (1..=3).map(|i| format!("Cell_{}", i)).collect();

        let result = embedding_to_matrix(&embedding, &cell_names);
        assert!(result.is_ok());

        match result.unwrap() {
            RObject::DoubleVector(dv) => {
                // Column-major: [1, 3, 5, 2, 4, 6]
                assert_eq!(dv.data, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
                if let Some(RObject::IntegerVector(int_vec)) = dv.attributes.get("dim") {
                    assert_eq!(int_vec.data, vec![3, 2]);
                }
            }
            _ => panic!("Expected Double vector"),
        }
    }

    #[test]
    fn test_construct_spatial_data() {
        use crate::ir::SpatialData;
        let coordinates = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut scale_factors = HashMap::new();
        scale_factors.insert("tissue_hires_scalef".to_string(), 0.08);
        scale_factors.insert("spot_diameter_fullres".to_string(), 89.43);
        let spatial = SpatialData::new(coordinates, 2, None, Some(scale_factors))
            .expect("Failed to create spatial data");

        let result = construct_spatial_data(&spatial, 3);
        assert!(result.is_ok());
        match result.unwrap() {
            RObject::GenericVector(gv) => {
                let names = gv.attributes.get_names();
                assert!(names.is_some());
                assert_eq!(names.unwrap()[0], "slice1");
            }
            _ => panic!("Expected List for images"),
        }
    }

    #[test]
    fn test_construct_spatial_coordinates() {
        use crate::ir::SpatialData;
        let coordinates = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let spatial =
            SpatialData::new(coordinates, 2, None, None).expect("Failed to create spatial data");

        let result = construct_spatial_coordinates(&spatial);
        assert!(result.is_ok());
        match result.unwrap() {
            RObject::DoubleVector(dv) => {
                assert_eq!(dv.data, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
                if let Some(RObject::IntegerVector(int_vec)) = dv.attributes.get("dim") {
                    assert_eq!(int_vec.data, vec![3, 2]);
                }
            }
            _ => panic!("Expected Double vector"),
        }
    }

    #[test]
    fn test_construct_scale_factors() {
        let mut scale_factors = HashMap::new();
        scale_factors.insert("tissue_hires_scalef".to_string(), 0.08);
        scale_factors.insert("spot_diameter_fullres".to_string(), 89.43);

        let result = construct_scale_factors(&scale_factors);
        assert!(result.is_ok());
        match result.unwrap() {
            RObject::GenericVector(gv) => {
                assert_eq!(gv.data.len(), 2);
            }
            _ => panic!("Expected List"),
        }
    }
}
