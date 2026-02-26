//! R object → IR converters
use extendr_api::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;
use crosscell::ir::{DataFrame, DatasetMetadata, DenseMatrix, Embedding, ExpressionMatrix, SingleCellData, SparseMatrixCSC};
use crosscell::ir::unstructured::UnstructuredValue;
use arrow::array::{ArrayRef, BooleanArray, Float64Array, Int64Array, StringArray, DictionaryArray, Int32Array};
use arrow::datatypes::Int32Type;

pub fn r_to_expression(obj: Robj) -> Result<ExpressionMatrix> {
    if obj.inherits("dgCMatrix") { dgcmatrix_to_expression(obj) }
    else if obj.is_matrix() { rmatrix_to_expression(obj) }
    else { Err("Expected dgCMatrix or matrix".into()) }
}

fn dgcmatrix_to_expression(obj: Robj) -> Result<ExpressionMatrix> {
    let i_slot: Vec<i32> = obj.dollar("i")?.as_integer_vector().unwrap();
    let p_slot: Vec<i32> = obj.dollar("p")?.as_integer_vector().unwrap();
    let x_slot: Vec<f64> = obj.dollar("x")?.as_real_vector().unwrap();
    let dim: Vec<i32> = obj.dollar("Dim")?.as_integer_vector().unwrap();
    let n_rows = dim[0] as usize;
    let n_cols = dim[1] as usize;
    let indices: Vec<usize> = i_slot.iter().map(|&x| x as usize).collect();
    let indptr: Vec<usize> = p_slot.iter().map(|&x| x as usize).collect();
    let csc = SparseMatrixCSC::new(x_slot, indices, indptr, n_rows, n_cols).map_err(|e| Error::Other(e))?;
    Ok(ExpressionMatrix::SparseCSC(csc))
}

fn rmatrix_to_expression(obj: Robj) -> Result<ExpressionMatrix> {
    let data: Vec<f64> = obj.as_real_vector().ok_or("Cannot convert matrix to real vector")?;
    let dim: Vec<i32> = obj.dim().ok_or("Matrix has no dimensions")?.iter().map(|x| x.inner()).collect();
    let n_rows = dim[0] as usize;
    let n_cols = dim[1] as usize;
    let mut row_major = vec![0.0f64; n_rows * n_cols];
    for row in 0..n_rows { for col in 0..n_cols { row_major[row * n_cols + col] = data[col * n_rows + row]; } }
    let dense = DenseMatrix::new(row_major, n_rows, n_cols).map_err(|e| Error::Other(e))?;
    Ok(ExpressionMatrix::Dense(dense))
}

pub fn r_to_dataframe(obj: Robj) -> Result<DataFrame> {
    if !obj.inherits("data.frame") { return Err("Expected data.frame".into()); }
    let obj_clone = obj.clone();
    let n_rows: i32 = R!("nrow({{obj_clone}})")?.as_integer().unwrap_or(0);
    let n_rows = n_rows as usize;
    let obj_clone2 = obj.clone();
    let col_names: Vec<String> = R!("colnames({{obj_clone2}})")?.as_string_vector().unwrap_or_default();
    let mut columns = Vec::new();
    let mut arrow_data: Vec<ArrayRef> = Vec::new();
    for col_name in &col_names {
        let col = obj.dollar(col_name.as_str())?;
        columns.push(col_name.clone());
        arrow_data.push(r_vector_to_arrow(&col)?);
    }
    DataFrame::new(columns, arrow_data, n_rows).map_err(|e| Error::Other(e))
}


fn r_vector_to_arrow(obj: &Robj) -> Result<ArrayRef> {
    if obj.inherits("factor") { return r_factor_to_arrow(obj); }
    match obj.rtype() {
        Rtype::Integers => {
            let values: Vec<i32> = obj.as_integer_vector().unwrap();
            let av: Vec<Option<i64>> = values.iter().map(|&v| if v == i32::MIN { None } else { Some(v as i64) }).collect();
            Ok(Arc::new(Int64Array::from(av)))
        }
        Rtype::Doubles => {
            let values: Vec<f64> = obj.as_real_vector().unwrap();
            let av: Vec<Option<f64>> = values.iter().map(|&v| if v.is_nan() { None } else { Some(v) }).collect();
            Ok(Arc::new(Float64Array::from(av)))
        }
        Rtype::Logicals => {
            let values: Vec<Rbool> = obj.as_logical_vector().unwrap();
            let av: Vec<Option<bool>> = values.iter().map(|&v| if v.is_na() { None } else { Some(v.is_true()) }).collect();
            Ok(Arc::new(BooleanArray::from(av)))
        }
        Rtype::Strings => {
            let values: Vec<String> = obj.as_string_vector().unwrap_or_default();
            let av: Vec<Option<&str>> = values.iter().map(|s| Some(s.as_str())).collect();
            Ok(Arc::new(StringArray::from(av)))
        }
        _ => {
            let obj_clone = obj.clone();
            let sv: Vec<String> = R!("as.character({{obj_clone}})")?.as_string_vector().unwrap_or_default();
            let av: Vec<Option<&str>> = sv.iter().map(|s| Some(s.as_str())).collect();
            Ok(Arc::new(StringArray::from(av)))
        }
    }
}

fn r_factor_to_arrow(obj: &Robj) -> Result<ArrayRef> {
    let obj_clone = obj.clone();
    let levels: Vec<String> = R!("levels({{obj_clone}})")?.as_string_vector().unwrap_or_default();
    let codes: Vec<i32> = obj.as_integer_vector().unwrap();
    let ac: Vec<Option<i32>> = codes.iter().map(|&v| if v == i32::MIN { None } else { Some(v - 1) }).collect();
    let keys = Int32Array::from(ac);
    let values = StringArray::from(levels.iter().map(|s| Some(s.as_str())).collect::<Vec<_>>());
    let da = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values)).map_err(|e| Error::Other(format!("{}", e)))?;
    Ok(Arc::new(da))
}

pub fn r_to_embedding(name: &str, obj: Robj) -> Result<Embedding> {
    let mat = if obj.inherits("DimReduc") { let o = obj.clone(); R!("{{o}}@cell.embeddings")? } else { obj };
    let data: Vec<f64> = mat.as_real_vector().ok_or("Cannot convert to real vector")?;
    let dim: Vec<i32> = mat.dim().ok_or("Matrix has no dimensions")?.iter().map(|x| x.inner()).collect();
    let n_rows = dim[0] as usize;
    let n_cols = dim[1] as usize;
    let mut rm = vec![0.0f64; n_rows * n_cols];
    for row in 0..n_rows { for col in 0..n_cols { rm[row * n_cols + col] = data[col * n_rows + row]; } }
    Embedding::new(name.to_string(), rm, n_rows, n_cols).map_err(|e| Error::Other(e))
}


pub fn seurat_to_ir(obj: Robj) -> Result<SingleCellData> {
    if !obj.inherits("Seurat") { return Err("Expected Seurat object".into()); }
    let o1 = obj.clone();
    let aa: String = R!("Seurat::DefaultAssay({{o1}})")?.as_str().unwrap_or("RNA").to_string();
    let aa1 = aa.clone();
    let o2 = obj.clone();
    let counts = R!("Seurat::GetAssayData({{o2}}, assay = {{aa1}}, layer = 'counts')")?;
    let expression = r_to_expression(counts)?;
    let o3 = obj.clone();
    let nc: i32 = R!("ncol({{o3}})")?.as_integer().unwrap_or(0);
    let o4 = obj.clone();
    let ng: i32 = R!("nrow({{o4}})")?.as_integer().unwrap_or(0);
    let (n_cells, n_genes) = (nc as usize, ng as usize);
    let o5 = obj.clone();
    let md = R!("{{o5}}@meta.data")?;
    let cell_metadata = r_to_dataframe(md)?;
    let aa2 = aa.clone();
    let o6 = obj.clone();
    let gm = R!("{{o6}}[[{{aa2}}]]@meta.data")?;
    let gene_metadata = if gm.is_null() { DataFrame::empty(n_genes) } else { r_to_dataframe(gm).unwrap_or_else(|_| DataFrame::empty(n_genes)) };
    let metadata = DatasetMetadata::new(n_cells, n_genes, "seurat".to_string());
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).map_err(|e| Error::Other(e))?;
    let o7 = obj.clone();
    let rn: Vec<String> = R!("names({{o7}}@reductions)")?.as_string_vector().unwrap_or_default();
    if !rn.is_empty() {
        let mut embs = HashMap::new();
        for name in &rn {
            let o8 = obj.clone();
            let dr = R!("{{o8}}@reductions[[{{name}}]]")?;
            if let Ok(emb) = r_to_embedding(name, dr) { embs.insert(name.clone(), emb); }
        }
        if !embs.is_empty() { data.embeddings = Some(embs); }
    }
    let aa3 = aa.clone();
    let o9 = obj.clone();
    let dm = R!("tryCatch(Seurat::GetAssayData({{o9}}, assay = {{aa3}}, layer = 'data'), error = function(e) NULL)")?;
    if !dm.is_null() && dm.len() > 0 {
        let mut layers = data.layers.clone().unwrap_or_default();
        if let Ok(expr) = r_to_expression(dm) { layers.insert("data".to_string(), expr); }
        if !layers.is_empty() { data.layers = Some(layers); }
    }
    Ok(data)
}


pub fn sce_to_ir(obj: Robj) -> Result<SingleCellData> {
    if !obj.inherits("SingleCellExperiment") { return Err("Expected SingleCellExperiment object".into()); }
    let o1 = obj.clone();
    let counts = R!("SummarizedExperiment::assay({{o1}}, 'counts')")?;
    let expression = r_to_expression(counts)?;
    let o2 = obj.clone();
    let nc: i32 = R!("ncol({{o2}})")?.as_integer().unwrap_or(0);
    let o3 = obj.clone();
    let ng: i32 = R!("nrow({{o3}})")?.as_integer().unwrap_or(0);
    let (n_cells, n_genes) = (nc as usize, ng as usize);
    let o4 = obj.clone();
    let cd = R!("as.data.frame(SummarizedExperiment::colData({{o4}}))")?;
    let cell_metadata = r_to_dataframe(cd)?;
    let o5 = obj.clone();
    let rd = R!("as.data.frame(SummarizedExperiment::rowData({{o5}}))")?;
    let gene_metadata = r_to_dataframe(rd)?;
    let metadata = DatasetMetadata::new(n_cells, n_genes, "sce".to_string());
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).map_err(|e| Error::Other(e))?;
    let o6 = obj.clone();
    let rdn: Vec<String> = R!("SingleCellExperiment::reducedDimNames({{o6}})")?.as_string_vector().unwrap_or_default();
    if !rdn.is_empty() {
        let mut embs = HashMap::new();
        for name in &rdn {
            let o7 = obj.clone();
            let mat = R!("SingleCellExperiment::reducedDim({{o7}}, {{name}})")?;
            if let Ok(emb) = r_to_embedding(name, mat) { embs.insert(name.clone(), emb); }
        }
        if !embs.is_empty() { data.embeddings = Some(embs); }
    }
    let o8 = obj.clone();
    let an: Vec<String> = R!("SummarizedExperiment::assayNames({{o8}})")?.as_string_vector().unwrap_or_default();
    let oa: Vec<&String> = an.iter().filter(|&n| n != "counts").collect();
    if !oa.is_empty() {
        let mut layers = HashMap::new();
        for name in oa {
            let o9 = obj.clone();
            let assay = R!("SummarizedExperiment::assay({{o9}}, {{name}})")?;
            if let Ok(expr) = r_to_expression(assay) { layers.insert(name.clone(), expr); }
        }
        if !layers.is_empty() { data.layers = Some(layers); }
    }
    Ok(data)
}


pub fn r_to_unstructured(obj: Robj) -> Result<UnstructuredValue> {
    if obj.is_null() || obj.is_na() { return Ok(UnstructuredValue::Null); }
    if obj.is_list() {
        let o1 = obj.clone();
        let names: Option<Vec<String>> = R!("names({{o1}})")?.as_string_vector();
        if let Some(names) = names {
            if !names.is_empty() {
                let mut map = HashMap::new();
                for name in &names {
                    let item = obj.dollar(name.as_str())?;
                    map.insert(name.clone(), r_to_unstructured(item)?);
                }
                return Ok(UnstructuredValue::Dict(map));
            }
        }
        let o2 = obj.clone();
        let len: i32 = R!("length({{o2}})")?.as_integer().unwrap_or(0);
        let mut arr = Vec::new();
        for i in 1..=len { let o3 = obj.clone(); let item = R!("{{o3}}[[{{i}}]]")?; arr.push(r_to_unstructured(item)?); }
        return Ok(UnstructuredValue::Array(arr));
    }
    match obj.rtype() {
        Rtype::Logicals => {
            let v: Vec<Rbool> = obj.as_logical_vector().unwrap();
            if v.len() == 1 { if v[0].is_na() { Ok(UnstructuredValue::Null) } else { Ok(UnstructuredValue::Boolean(v[0].is_true())) } }
            else { Ok(UnstructuredValue::Array(v.iter().map(|&x| if x.is_na() { UnstructuredValue::Null } else { UnstructuredValue::Boolean(x.is_true()) }).collect())) }
        }
        Rtype::Integers => {
            let v: Vec<i32> = obj.as_integer_vector().unwrap();
            if v.len() == 1 { if v[0] == i32::MIN { Ok(UnstructuredValue::Null) } else { Ok(UnstructuredValue::Integer(v[0] as i64)) } }
            else { Ok(UnstructuredValue::Array(v.iter().map(|&x| if x == i32::MIN { UnstructuredValue::Null } else { UnstructuredValue::Integer(x as i64) }).collect())) }
        }
        Rtype::Doubles => {
            let v: Vec<f64> = obj.as_real_vector().unwrap();
            if v.len() == 1 { if v[0].is_nan() { Ok(UnstructuredValue::Null) } else { Ok(UnstructuredValue::Float(v[0])) } }
            else { Ok(UnstructuredValue::Array(v.iter().map(|&x| if x.is_nan() { UnstructuredValue::Null } else { UnstructuredValue::Float(x) }).collect())) }
        }
        Rtype::Strings => {
            let v: Vec<String> = obj.as_string_vector().unwrap_or_default();
            if v.len() == 1 { Ok(UnstructuredValue::String(v[0].clone())) }
            else { Ok(UnstructuredValue::Array(v.iter().map(|s| UnstructuredValue::String(s.clone())).collect())) }
        }
        _ => { let o = obj.clone(); let s: String = R!("as.character({{o}})[[1]]")?.as_str().unwrap_or("").to_string(); Ok(UnstructuredValue::String(s)) }
    }
}
