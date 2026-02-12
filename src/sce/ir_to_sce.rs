//! IR to SingleCellExperiment conversion
//!
//! This module provides functions to convert CrossCell IR to SingleCellExperiment RDS format.

use crate::ir::{ExpressionMatrix, SingleCellData, DataFrame as IrDataFrame, Embedding};
use crate::rds::{
    RObject, RdsFile, Attributes, GenericVector, S4Object,
    IntegerVector, DoubleVector, StringVector,
    write_rds,
};
use crate::rds::StringEncoding;
use crate::sce::error::{SceError, Result};
use std::collections::HashMap;
use std::path::Path;

/// Convert IR to SingleCellExperiment RDS structure
pub fn ir_to_sce_rds(ir: &SingleCellData) -> Result<RdsFile> {
    log::debug!("Converting IR to SingleCellExperiment");
    
    let assays = build_assays(ir)?;
    let col_data = build_col_data(&ir.cell_metadata)?;
    let element_metadata = build_row_data(&ir.gene_metadata)?;
    let int_col_data = build_int_col_data(&ir.embeddings, ir.metadata.n_cells)?;
    let sce = build_sce_object(assays, col_data, element_metadata, int_col_data, ir)?;
    
    let mut file = RdsFile::default();
    file.object = sce;
    Ok(file)
}

/// Write IR to SingleCellExperiment RDS file
pub fn write_sce_rds(ir: &SingleCellData, path: &str) -> Result<()> {
    let rds_file = ir_to_sce_rds(ir)?;
    write_rds(&rds_file, Path::new(path))
        .map_err(|e| SceError::IoError(format!("Failed to write RDS: {}", e)))?;
    Ok(())
}

fn attrs_with_names(names: Vec<String>) -> Attributes {
    let mut attrs = Attributes::new();
    let len = names.len();
    attrs.add(
        "names".to_string(),
        RObject::StringVector(StringVector {
            data: names,
            encodings: vec![StringEncoding::Utf8; len],
            missing: vec![false; len],
            attributes: Attributes::new(),
        }),
        StringEncoding::Utf8,
    );
    attrs
}

fn set_attr(attrs: &mut Attributes, name: &str, value: RObject) {
    attrs.add(name.to_string(), value, StringEncoding::Utf8);
}

fn build_assays(ir: &SingleCellData) -> Result<RObject> {
    let mut assay_list = Vec::new();
    let mut assay_names = Vec::new();
    
    let counts_matrix = build_dgc_matrix(&ir.expression)?;
    assay_list.push(counts_matrix);
    assay_names.push("counts".to_string());
    
    if let Some(ref layers) = ir.layers {
        for (name, matrix) in layers {
            let layer_matrix = build_dgc_matrix(matrix)?;
            assay_list.push(layer_matrix);
            assay_names.push(name.clone());
        }
    }
    
    // Build the inner list with assay data
    let list_data = RObject::GenericVector(GenericVector {
        data: assay_list,
        attributes: attrs_with_names(assay_names),
    });
    
    // Build SimpleList (data slot of SimpleAssays)
    let mut simple_list_attrs = Attributes::new();
    set_attr(&mut simple_list_attrs, "listData", list_data);
    
    let simple_list = RObject::S4Object(S4Object {
        class_name: "SimpleList".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "S4Vectors".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: simple_list_attrs,
    });
    
    // Build SimpleAssays with data slot containing SimpleList
    let mut assays_attrs = Attributes::new();
    set_attr(&mut assays_attrs, "data", simple_list);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "SimpleAssays".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "SummarizedExperiment".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: assays_attrs,
    }))
}

fn build_dgc_matrix(expr: &ExpressionMatrix) -> Result<RObject> {
    let (n_cells, n_genes) = expr.shape();
    
    let (col_ptrs, row_indices, data): (Vec<usize>, Vec<usize>, Vec<f64>) = match expr {
        ExpressionMatrix::SparseCSC(csc) => transpose_csc_to_dgc(csc),
        ExpressionMatrix::SparseCSR(csr) => {
            let csc = crate::sparse::convert::csr_to_csc(csr);
            transpose_csc_to_dgc(&csc)
        }
        ExpressionMatrix::Dense(dense) => {
            let csc = dense_to_csc(dense, n_cells, n_genes);
            transpose_csc_to_dgc(&csc)
        }
        ExpressionMatrix::Lazy(_) => {
            return Err(SceError::ConversionError(
                "Lazy matrices must be loaded before conversion to SCE".to_string()
            ));
        }
    };
    
    let mut attrs = Attributes::new();
    
    set_attr(&mut attrs, "i", RObject::IntegerVector(IntegerVector {
        data: row_indices.iter().map(|&x| x as i32).collect(),
        attributes: Attributes::new(),
    }));
    
    set_attr(&mut attrs, "p", RObject::IntegerVector(IntegerVector {
        data: col_ptrs.iter().map(|&x| x as i32).collect(),
        attributes: Attributes::new(),
    }));
    
    set_attr(&mut attrs, "x", RObject::DoubleVector(DoubleVector {
        data,
        attributes: Attributes::new(),
    }));
    
    set_attr(&mut attrs, "Dim", RObject::IntegerVector(IntegerVector {
        data: vec![n_genes as i32, n_cells as i32],
        attributes: Attributes::new(),
    }));
    
    let dimnames = RObject::GenericVector(GenericVector {
        data: vec![RObject::Null, RObject::Null],
        attributes: Attributes::new(),
    });
    set_attr(&mut attrs, "Dimnames", dimnames);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "dgCMatrix".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "Matrix".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: attrs,
    }))
}

fn transpose_csc_to_dgc(csc: &crate::ir::SparseMatrixCSC) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let n_cells = csc.n_rows;
    let n_genes = csc.n_cols;
    let nnz = csc.data.len();
    
    let mut row_counts = vec![0usize; n_cells];
    for &row_idx in &csc.indices {
        row_counts[row_idx] += 1;
    }
    
    let mut col_ptrs = vec![0usize; n_cells + 1];
    for i in 0..n_cells {
        col_ptrs[i + 1] = col_ptrs[i] + row_counts[i];
    }
    
    let mut row_indices = vec![0usize; nnz];
    let mut data = vec![0.0; nnz];
    let mut current_pos = col_ptrs.clone();
    
    for col in 0..n_genes {
        let start = csc.indptr[col];
        let end = csc.indptr[col + 1];
        for idx in start..end {
            let row = csc.indices[idx];
            let val = csc.data[idx];
            let pos = current_pos[row];
            row_indices[pos] = col;
            data[pos] = val;
            current_pos[row] += 1;
        }
    }
    
    (col_ptrs, row_indices, data)
}

fn dense_to_csc(dense: &crate::ir::DenseMatrix, n_rows: usize, n_cols: usize) -> crate::ir::SparseMatrixCSC {
    let mut indptr = vec![0usize; n_cols + 1];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    
    for col in 0..n_cols {
        for row in 0..n_rows {
            let val = dense.data[row * n_cols + col];
            if val != 0.0 {
                indices.push(row);
                data.push(val);
            }
        }
        indptr[col + 1] = indices.len();
    }
    
    crate::ir::SparseMatrixCSC { n_rows, n_cols, indptr, indices, data }
}

fn build_col_data(df: &IrDataFrame) -> Result<RObject> {
    build_dframe(df)
}

fn build_row_data(df: &IrDataFrame) -> Result<RObject> {
    build_dframe(df)
}

fn build_dframe(df: &IrDataFrame) -> Result<RObject> {
    let n_rows = df.n_rows;
    let mut col_list = Vec::new();
    let mut col_names = Vec::new();
    
    for (name, array) in df.columns() {
        let robj = arrow_array_to_robject(array)?;
        col_list.push(robj);
        col_names.push(name.clone());
    }
    
    let list_data = RObject::GenericVector(GenericVector {
        data: col_list,
        attributes: attrs_with_names(col_names),
    });
    
    let mut attrs = Attributes::new();
    set_attr(&mut attrs, "listData", list_data);
    set_attr(&mut attrs, "nrows", RObject::IntegerVector(IntegerVector {
        data: vec![n_rows as i32],
        attributes: Attributes::new(),
    }));
    set_attr(&mut attrs, "rownames", RObject::Null);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "DFrame".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "S4Vectors".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: attrs,
    }))
}


fn arrow_array_to_robject(array: &arrow::array::ArrayRef) -> Result<RObject> {
    use arrow::array::*;
    use arrow::datatypes::DataType;
    
    match array.data_type() {
        DataType::Int32 => {
            let arr = array.as_any().downcast_ref::<Int32Array>().unwrap();
            let data: Vec<i32> = arr.iter().map(|v| v.unwrap_or(i32::MIN)).collect();
            Ok(RObject::IntegerVector(IntegerVector { data, attributes: Attributes::new() }))
        }
        DataType::Int64 => {
            let arr = array.as_any().downcast_ref::<Int64Array>().unwrap();
            let data: Vec<f64> = arr.iter().map(|v| v.unwrap_or(0) as f64).collect();
            Ok(RObject::DoubleVector(DoubleVector { data, attributes: Attributes::new() }))
        }
        DataType::Float64 => {
            let arr = array.as_any().downcast_ref::<Float64Array>().unwrap();
            let data: Vec<f64> = arr.iter().map(|v| v.unwrap_or(f64::NAN)).collect();
            Ok(RObject::DoubleVector(DoubleVector { data, attributes: Attributes::new() }))
        }
        DataType::Utf8 => {
            let arr = array.as_any().downcast_ref::<StringArray>().unwrap();
            let data: Vec<String> = arr.iter().map(|v| v.unwrap_or("").to_string()).collect();
            let len = data.len();
            Ok(RObject::StringVector(StringVector { 
                data, 
                encodings: vec![StringEncoding::Utf8; len], 
                missing: vec![false; len], 
                attributes: Attributes::new() 
            }))
        }
        DataType::Boolean => {
            let arr = array.as_any().downcast_ref::<BooleanArray>().unwrap();
            let data: Vec<i32> = arr.iter().map(|v| if v.unwrap_or(false) { 1 } else { 0 }).collect();
            Ok(RObject::IntegerVector(IntegerVector { data, attributes: Attributes::new() }))
        }
        DataType::Dictionary(_, _) => arrow_dict_to_factor(array),
        _ => Err(SceError::ConversionError(format!("Unsupported Arrow type: {:?}", array.data_type())))
    }
}

fn arrow_dict_to_factor(array: &arrow::array::ArrayRef) -> Result<RObject> {
    use arrow::array::*;
    use arrow::datatypes::Int32Type;
    
    let dict_array = array.as_any()
        .downcast_ref::<DictionaryArray<Int32Type>>()
        .ok_or_else(|| SceError::ConversionError("Expected Int32 dictionary".to_string()))?;
    
    let keys: Vec<i32> = dict_array.keys().iter()
        .map(|v| v.map(|k| k + 1).unwrap_or(i32::MIN))
        .collect();
    
    let values = dict_array.values();
    let levels: Vec<String> = if let Some(str_arr) = values.as_any().downcast_ref::<StringArray>() {
        str_arr.iter().map(|v| v.unwrap_or("").to_string()).collect()
    } else {
        return Err(SceError::ConversionError("Dictionary values must be strings".to_string()));
    };
    
    let levels_len = levels.len();
    let mut attrs = Attributes::new();
    set_attr(&mut attrs, "levels", RObject::StringVector(StringVector { 
        data: levels, 
        encodings: vec![StringEncoding::Utf8; levels_len], 
        missing: vec![false; levels_len], 
        attributes: Attributes::new() 
    }));
    set_attr(&mut attrs, "class", RObject::StringVector(StringVector { 
        data: vec!["factor".to_string()], 
        encodings: vec![StringEncoding::Utf8], 
        missing: vec![false], 
        attributes: Attributes::new() 
    }));
    
    Ok(RObject::IntegerVector(IntegerVector { data: keys, attributes: attrs }))
}

fn build_int_col_data(embeddings: &Option<HashMap<String, Embedding>>, n_cells: usize) -> Result<RObject> {
    let reduced_dims = if let Some(embs) = embeddings {
        build_reduced_dims(embs)?
    } else {
        let mut attrs = Attributes::new();
        set_attr(&mut attrs, "listData", RObject::GenericVector(GenericVector {
            data: vec![], attributes: Attributes::new(),
        }));
        RObject::S4Object(S4Object {
            class_name: "SimpleList".to_string(),
            class_encoding: StringEncoding::Utf8,
            package_name: "S4Vectors".to_string(),
            package_encoding: StringEncoding::Utf8,
            attributes: attrs,
        })
    };
    
    let mut attrs = Attributes::new();
    set_attr(&mut attrs, "listData", RObject::GenericVector(GenericVector {
        data: vec![reduced_dims],
        attributes: attrs_with_names(vec!["reducedDims".to_string()]),
    }));
    set_attr(&mut attrs, "nrows", RObject::IntegerVector(IntegerVector {
        data: vec![n_cells as i32], attributes: Attributes::new(),
    }));
    set_attr(&mut attrs, "rownames", RObject::Null);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "DFrame".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "S4Vectors".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: attrs,
    }))
}

fn build_reduced_dims(embeddings: &HashMap<String, Embedding>) -> Result<RObject> {
    let mut dim_list = Vec::new();
    let mut dim_names = Vec::new();
    
    for (name, emb) in embeddings {
        let matrix = build_matrix(emb)?;
        dim_list.push(matrix);
        dim_names.push(name.clone());
    }
    
    let list_data = RObject::GenericVector(GenericVector {
        data: dim_list,
        attributes: attrs_with_names(dim_names),
    });
    
    let mut attrs = Attributes::new();
    set_attr(&mut attrs, "listData", list_data);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "SimpleList".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "S4Vectors".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: attrs,
    }))
}

fn build_matrix(emb: &Embedding) -> Result<RObject> {
    let n_rows = emb.n_rows;
    let n_cols = emb.n_cols;
    
    let mut col_major_data = vec![0.0; n_rows * n_cols];
    for row in 0..n_rows {
        for col in 0..n_cols {
            col_major_data[col * n_rows + row] = emb.data[row * n_cols + col];
        }
    }
    
    let mut attrs = Attributes::new();
    set_attr(&mut attrs, "dim", RObject::IntegerVector(IntegerVector {
        data: vec![n_rows as i32, n_cols as i32], attributes: Attributes::new(),
    }));
    
    Ok(RObject::DoubleVector(DoubleVector { data: col_major_data, attributes: attrs }))
}

fn build_sce_object(
    assays: RObject,
    col_data: RObject,
    element_metadata: RObject,
    int_col_data: RObject,
    ir: &SingleCellData,
) -> Result<RObject> {
    let mut attrs = Attributes::new();
    
    set_attr(&mut attrs, "assays", assays);
    set_attr(&mut attrs, "colData", col_data);
    set_attr(&mut attrs, "elementMetadata", element_metadata);
    set_attr(&mut attrs, "int_colData", int_col_data);
    
    let mut int_elem_attrs = Attributes::new();
    set_attr(&mut int_elem_attrs, "listData", RObject::GenericVector(GenericVector {
        data: vec![], attributes: Attributes::new(),
    }));
    set_attr(&mut int_elem_attrs, "nrows", RObject::IntegerVector(IntegerVector {
        data: vec![ir.metadata.n_genes as i32], attributes: Attributes::new(),
    }));
    set_attr(&mut int_elem_attrs, "rownames", RObject::Null);
    set_attr(&mut attrs, "int_elementMetadata", RObject::S4Object(S4Object {
        class_name: "DFrame".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "S4Vectors".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: int_elem_attrs,
    }));
    
    set_attr(&mut attrs, "int_metadata", RObject::GenericVector(GenericVector {
        data: vec![], attributes: Attributes::new(),
    }));
    set_attr(&mut attrs, "metadata", RObject::GenericVector(GenericVector {
        data: vec![], attributes: Attributes::new(),
    }));
    set_attr(&mut attrs, "NAMES", RObject::Null);
    set_attr(&mut attrs, "rowRanges", RObject::Null);
    
    Ok(RObject::S4Object(S4Object {
        class_name: "SingleCellExperiment".to_string(),
        class_encoding: StringEncoding::Utf8,
        package_name: "SingleCellExperiment".to_string(),
        package_encoding: StringEncoding::Utf8,
        attributes: attrs,
    }))
}
