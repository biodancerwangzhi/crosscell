//! AnnData -> IR converter
//!
//! Converts Python AnnData objects to CrossCell IR (SingleCellData).

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

use crosscell::ir::{
    DataFrame, DatasetMetadata, DenseMatrix, Embedding, ExpressionMatrix,
    PairwiseMatrix, SingleCellData, SparseMatrixCSC, SparseMatrixCSR,
};
use crosscell::ir::unstructured::UnstructuredValue;

/// Convert a Python AnnData object to SingleCellData
pub fn anndata_to_ir(py: Python<'_>, adata: &Bound<'_, PyAny>) -> PyResult<SingleCellData> {
    let shape = adata.getattr("shape")?;
    let n_cells: usize = shape.get_item(0)?.extract()?;
    let n_genes: usize = shape.get_item(1)?.extract()?;

    let x_obj = adata.getattr("X")?;
    let expression = python_to_expression(py, &x_obj)?;

    let obs_obj = adata.getattr("obs")?;
    let cell_metadata = pandas_to_dataframe(py, &obs_obj, n_cells)?;

    let var_obj = adata.getattr("var")?;
    let gene_metadata = pandas_to_dataframe(py, &var_obj, n_genes)?;

    let metadata = DatasetMetadata::new(n_cells, n_genes, "anndata".to_string());
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    // Helper to convert keys() to list
    let builtins = py.import_bound("builtins")?;

    // layers
    let layers_obj = adata.getattr("layers")?;
    let layers_keys: Vec<String> = builtins
        .call_method1("list", (layers_obj.call_method0("keys")?,))?
        .extract()?;
    if !layers_keys.is_empty() {
        let mut layers = HashMap::new();
        for key in &layers_keys {
            let layer = layers_obj.get_item(key)?;
            layers.insert(key.clone(), python_to_expression(py, &layer)?);
        }
        data.layers = Some(layers);
    }

    // obsm
    let obsm_obj = adata.getattr("obsm")?;
    let obsm_keys: Vec<String> = builtins
        .call_method1("list", (obsm_obj.call_method0("keys")?,))?
        .extract()?;
    if !obsm_keys.is_empty() {
        let mut embeddings = HashMap::new();
        for key in &obsm_keys {
            let arr = obsm_obj.get_item(key)?;
            embeddings.insert(key.clone(), numpy_to_embedding(py, key, &arr)?);
        }
        data.embeddings = Some(embeddings);
    }

    // obsp
    let obsp_obj = adata.getattr("obsp")?;
    let obsp_keys: Vec<String> = builtins
        .call_method1("list", (obsp_obj.call_method0("keys")?,))?
        .extract()?;
    if !obsp_keys.is_empty() {
        let mut cell_pairwise = HashMap::new();
        for key in &obsp_keys {
            let mat = obsp_obj.get_item(key)?;
            let matrix = python_to_expression(py, &mat)?;
            cell_pairwise.insert(
                key.clone(),
                PairwiseMatrix { name: key.clone(), matrix },
            );
        }
        data.cell_pairwise = Some(cell_pairwise);
    }

    // varm
    let varm_obj = adata.getattr("varm")?;
    let varm_keys: Vec<String> = builtins
        .call_method1("list", (varm_obj.call_method0("keys")?,))?
        .extract()?;
    if !varm_keys.is_empty() {
        let mut gene_loadings = HashMap::new();
        for key in &varm_keys {
            let arr = varm_obj.get_item(key)?;
            gene_loadings.insert(key.clone(), numpy_to_embedding(py, key, &arr)?);
        }
        data.gene_loadings = Some(gene_loadings);
    }

    // uns
    let uns_obj = adata.getattr("uns")?;
    let uns_keys: Vec<String> = builtins
        .call_method1("list", (uns_obj.call_method0("keys")?,))?
        .extract()?;
    if !uns_keys.is_empty() {
        let mut uns = HashMap::new();
        for key in &uns_keys {
            let val = uns_obj.get_item(key)?;
            uns.insert(key.clone(), python_to_unstructured(py, &val)?);
        }
        data.unstructured = Some(uns);
    }

    Ok(data)
}

/// Convert scipy.sparse / numpy.ndarray to ExpressionMatrix
pub fn python_to_expression(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
) -> PyResult<ExpressionMatrix> {
    let scipy_sparse = py.import_bound("scipy.sparse")?;
    let is_sparse: bool = scipy_sparse
        .call_method1("issparse", (obj,))?
        .extract()?;

    if is_sparse {
        sparse_to_expression(py, obj)
    } else {
        dense_to_expression(py, obj)
    }
}

fn sparse_to_expression(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
) -> PyResult<ExpressionMatrix> {
    let format: String = obj.getattr("format")?.extract()?;
    match format.as_str() {
        "csr" => {
            let data: Vec<f64> = obj.getattr("data")?.call_method0("tolist")?.extract()?;
            let indices: Vec<i64> = obj.getattr("indices")?.call_method0("tolist")?.extract()?;
            let indptr: Vec<i64> = obj.getattr("indptr")?.call_method0("tolist")?.extract()?;
            let shape = obj.getattr("shape")?;
            let n_rows: usize = shape.get_item(0)?.extract()?;
            let n_cols: usize = shape.get_item(1)?.extract()?;
            let m = SparseMatrixCSR::new(
                data,
                indices.into_iter().map(|i| i as usize).collect(),
                indptr.into_iter().map(|i| i as usize).collect(),
                n_rows,
                n_cols,
            )
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
            Ok(ExpressionMatrix::SparseCSR(m))
        }
        "csc" => {
            let data: Vec<f64> = obj.getattr("data")?.call_method0("tolist")?.extract()?;
            let indices: Vec<i64> = obj.getattr("indices")?.call_method0("tolist")?.extract()?;
            let indptr: Vec<i64> = obj.getattr("indptr")?.call_method0("tolist")?.extract()?;
            let shape = obj.getattr("shape")?;
            let n_rows: usize = shape.get_item(0)?.extract()?;
            let n_cols: usize = shape.get_item(1)?.extract()?;
            let m = SparseMatrixCSC::new(
                data,
                indices.into_iter().map(|i| i as usize).collect(),
                indptr.into_iter().map(|i| i as usize).collect(),
                n_rows,
                n_cols,
            )
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
            Ok(ExpressionMatrix::SparseCSC(m))
        }
        _ => {
            let csr = obj.call_method0("tocsr")?;
            sparse_to_expression(py, &csr)
        }
    }
}

fn dense_to_expression(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
) -> PyResult<ExpressionMatrix> {
    let numpy = py.import_bound("numpy")?;
    let arr = numpy.call_method1("asarray", (obj,))?;
    let flat = arr.call_method0("flatten")?.call_method0("tolist")?;
    let data: Vec<f64> = flat.extract()?;
    let shape = arr.getattr("shape")?;
    let n_rows: usize = shape.get_item(0)?.extract()?;
    let n_cols: usize = shape.get_item(1)?.extract()?;
    let m = DenseMatrix::new(data, n_rows, n_cols)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
    Ok(ExpressionMatrix::Dense(m))
}

/// Convert pandas.DataFrame to DataFrame_IR
pub fn pandas_to_dataframe(
    py: Python<'_>,
    df: &Bound<'_, PyAny>,
    n_rows: usize,
) -> PyResult<DataFrame> {
    use arrow::array::{BooleanArray, Float64Array, Int64Array, StringArray};
    use std::sync::Arc;

    let col_names: Vec<String> = df
        .getattr("columns")?
        .call_method0("tolist")?
        .extract()?;
    let mut columns = Vec::new();
    let mut arrow_data = Vec::new();

    for col_name in &col_names {
        let series = df.get_item(col_name)?;
        let dtype_str: String = series
            .getattr("dtype")?
            .call_method0("__str__")?
            .extract()?;

        let arrow_array: arrow::array::ArrayRef = if dtype_str == "category" {
            pandas_categorical_to_arrow(py, &series)?
        } else if dtype_str.starts_with("int") {
            let values: Vec<Option<i64>> = series.call_method0("tolist")?.extract()?;
            Arc::new(Int64Array::from(values))
        } else if dtype_str.starts_with("float") {
            let values: Vec<f64> = series
                .call_method1("fillna", (0.0,))?
                .call_method0("tolist")?
                .extract()?;
            Arc::new(Float64Array::from(values))
        } else if dtype_str == "bool" {
            let values: Vec<Option<bool>> = series.call_method0("tolist")?.extract()?;
            Arc::new(BooleanArray::from(values))
        } else {
            let str_series = series.call_method1("astype", ("str",))?;
            let values: Vec<String> = str_series.call_method0("tolist")?.extract()?;
            let str_values: Vec<Option<&str>> =
                values.iter().map(|s| Some(s.as_str())).collect();
            Arc::new(StringArray::from(str_values))
        };

        columns.push(col_name.clone());
        arrow_data.push(arrow_array);
    }

    DataFrame::new(columns, arrow_data, n_rows)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))
}

fn pandas_categorical_to_arrow(
    _py: Python<'_>,
    series: &Bound<'_, PyAny>,
) -> PyResult<arrow::array::ArrayRef> {
    use arrow::array::{DictionaryArray, Int32Array, StringArray};
    use std::sync::Arc;

    let cat = series.getattr("cat")?;
    let codes: Vec<i32> = cat.getattr("codes")?.call_method0("tolist")?.extract()?;
    let categories: Vec<String> = cat
        .getattr("categories")?
        .call_method0("tolist")?
        .extract()?;

    let cat_array = StringArray::from(
        categories
            .iter()
            .map(|s| Some(s.as_str()))
            .collect::<Vec<_>>(),
    );
    let keys = Int32Array::from(codes);
    let dict_array = DictionaryArray::try_new(keys, Arc::new(cat_array))
        .map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Arrow error: {}", e))
        })?;
    Ok(Arc::new(dict_array))
}

/// Convert numpy.ndarray to Embedding
pub fn numpy_to_embedding(
    py: Python<'_>,
    name: &str,
    arr: &Bound<'_, PyAny>,
) -> PyResult<Embedding> {
    let numpy = py.import_bound("numpy")?;
    let arr = numpy.call_method1("asarray", (arr,))?;
    let shape = arr.getattr("shape")?;
    let n_rows: usize = shape.get_item(0)?.extract()?;
    let n_cols: usize = shape.get_item(1)?.extract()?;
    let flat = arr.call_method0("flatten")?.call_method0("tolist")?;
    let data: Vec<f64> = flat.extract()?;
    Embedding::new(name.to_string(), data, n_rows, n_cols)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))
}

/// Convert Python object to UnstructuredValue
pub fn python_to_unstructured(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
) -> PyResult<UnstructuredValue> {
    if obj.is_none() {
        return Ok(UnstructuredValue::Null);
    }
    if let Ok(b) = obj.extract::<bool>() {
        return Ok(UnstructuredValue::Boolean(b));
    }
    if let Ok(i) = obj.extract::<i64>() {
        return Ok(UnstructuredValue::Integer(i));
    }
    if let Ok(f) = obj.extract::<f64>() {
        return Ok(UnstructuredValue::Float(f));
    }
    if let Ok(s) = obj.extract::<String>() {
        return Ok(UnstructuredValue::String(s));
    }
    if let Ok(dict) = obj.downcast::<PyDict>() {
        let mut map = HashMap::new();
        for (k, v) in dict.iter() {
            let key: String = k.extract()?;
            let val = python_to_unstructured(py, &v)?;
            map.insert(key, val);
        }
        return Ok(UnstructuredValue::Dict(map));
    }
    if let Ok(items) = obj.call_method0("__iter__") {
        let mut arr = Vec::new();
        loop {
            match items.call_method0("__next__") {
                Ok(item) => arr.push(python_to_unstructured(py, &item)?),
                Err(_) => break,
            }
        }
        return Ok(UnstructuredValue::Array(arr));
    }
    let s: String = obj.str()?.extract()?;
    Ok(UnstructuredValue::String(s))
}
