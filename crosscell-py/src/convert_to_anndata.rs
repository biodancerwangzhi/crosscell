//! IR → AnnData 转换器
//!
//! 将 CrossCell IR (SingleCellData) 转换为 Python AnnData 对象。

use arrow::array::{
    Array, ArrayRef, BooleanArray, DictionaryArray, Float32Array, Float64Array,
    Int32Array, Int64Array, LargeStringArray, StringArray,
};
use arrow::datatypes::{DataType, Int32Type, Int8Type};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::collections::HashMap;

use crosscell::ir::{
    DataFrame, DenseMatrix, Embedding, ExpressionMatrix,
    SingleCellData, SparseMatrixCSC, SparseMatrixCSR,
};
use crosscell::ir::unstructured::UnstructuredValue;

/// 将 SingleCellData 转换为 Python AnnData 对象
pub fn ir_to_anndata(py: Python<'_>, data: &SingleCellData) -> PyResult<PyObject> {
    let anndata_mod = py.import_bound("anndata")?;
    let kwargs = PyDict::new_bound(py);

    kwargs.set_item("X", expression_to_python(py, &data.expression)?)?;
    kwargs.set_item("obs", dataframe_to_pandas(py, &data.cell_metadata)?)?;
    kwargs.set_item("var", dataframe_to_pandas(py, &data.gene_metadata)?)?;

    if let Some(ref layers) = data.layers {
        let d = PyDict::new_bound(py);
        for (name, matrix) in layers {
            d.set_item(name, expression_to_python(py, matrix)?)?;
        }
        kwargs.set_item("layers", d)?;
    }

    if let Some(ref embeddings) = data.embeddings {
        let d = PyDict::new_bound(py);
        for (name, emb) in embeddings {
            d.set_item(name, embedding_to_numpy(py, emb)?)?;
        }
        kwargs.set_item("obsm", d)?;
    }

    if let Some(ref cell_pairwise) = data.cell_pairwise {
        let d = PyDict::new_bound(py);
        for (name, pw) in cell_pairwise {
            d.set_item(name, expression_to_python(py, &pw.matrix)?)?;
        }
        kwargs.set_item("obsp", d)?;
    }

    if let Some(ref gene_pairwise) = data.gene_pairwise {
        let d = PyDict::new_bound(py);
        for (name, pw) in gene_pairwise {
            d.set_item(name, expression_to_python(py, &pw.matrix)?)?;
        }
        kwargs.set_item("varp", d)?;
    }

    if let Some(ref gene_loadings) = data.gene_loadings {
        let d = PyDict::new_bound(py);
        for (name, loading) in gene_loadings {
            let numpy = py.import_bound("numpy")?;
            let arr = numpy.call_method1(
                "array",
                (loading.data.clone(),),
            )?;
            let reshaped = arr.call_method1("reshape", ((loading.n_rows, loading.n_cols),))?;
            d.set_item(name, reshaped)?;
        }
        kwargs.set_item("varm", d)?;
    }

    if let Some(ref uns) = data.unstructured {
        kwargs.set_item("uns", unstructured_map_to_dict(py, uns)?)?;
    }

    let adata = anndata_mod.getattr("AnnData")?.call((), Some(&kwargs))?;
    Ok(adata.into())
}

/// 将 ExpressionMatrix 转换为 scipy.sparse 或 numpy.ndarray
pub fn expression_to_python(py: Python<'_>, matrix: &ExpressionMatrix) -> PyResult<PyObject> {
    match matrix {
        ExpressionMatrix::SparseCSR(m) => csr_to_scipy(py, m),
        ExpressionMatrix::SparseCSC(m) => csc_to_scipy(py, m),
        ExpressionMatrix::Dense(m) => dense_to_numpy(py, m),
        ExpressionMatrix::Lazy(_) => Err(pyo3::exceptions::PyValueError::new_err(
            "Cannot convert lazy matrix to Python; load data first",
        )),
    }
}

fn csr_to_scipy(py: Python<'_>, m: &SparseMatrixCSR) -> PyResult<PyObject> {
    let scipy_sparse = py.import_bound("scipy.sparse")?;
    let numpy = py.import_bound("numpy")?;

    let data = numpy.call_method1("array", (m.data.clone(),))?;
    let indices = numpy.call_method1(
        "array",
        (m.indices.iter().map(|&i| i as i64).collect::<Vec<_>>(),),
    )?;
    let indptr = numpy.call_method1(
        "array",
        (m.indptr.iter().map(|&i| i as i64).collect::<Vec<_>>(),),
    )?;

    let tuple = PyTuple::new_bound(py, &[data.as_ref(), indices.as_ref(), indptr.as_ref()]);
    let shape = PyTuple::new_bound(py, &[m.n_rows as i64, m.n_cols as i64]);
    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("shape", shape)?;

    let result = scipy_sparse.getattr("csr_matrix")?.call((tuple,), Some(&kwargs))?;
    Ok(result.into())
}

fn csc_to_scipy(py: Python<'_>, m: &SparseMatrixCSC) -> PyResult<PyObject> {
    let scipy_sparse = py.import_bound("scipy.sparse")?;
    let numpy = py.import_bound("numpy")?;

    let data = numpy.call_method1("array", (m.data.clone(),))?;
    let indices = numpy.call_method1(
        "array",
        (m.indices.iter().map(|&i| i as i64).collect::<Vec<_>>(),),
    )?;
    let indptr = numpy.call_method1(
        "array",
        (m.indptr.iter().map(|&i| i as i64).collect::<Vec<_>>(),),
    )?;

    let tuple = PyTuple::new_bound(py, &[data.as_ref(), indices.as_ref(), indptr.as_ref()]);
    let shape = PyTuple::new_bound(py, &[m.n_rows as i64, m.n_cols as i64]);
    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("shape", shape)?;

    let result = scipy_sparse.getattr("csc_matrix")?.call((tuple,), Some(&kwargs))?;
    Ok(result.into())
}

fn dense_to_numpy(py: Python<'_>, m: &DenseMatrix) -> PyResult<PyObject> {
    let numpy = py.import_bound("numpy")?;
    let flat = numpy.call_method1("array", (m.data.clone(),))?;
    let shape = PyTuple::new_bound(py, &[m.n_rows as i64, m.n_cols as i64]);
    let reshaped = flat.call_method1("reshape", (shape,))?;
    Ok(reshaped.into())
}

/// 将 DataFrame_IR 转换为 pandas.DataFrame
pub fn dataframe_to_pandas(py: Python<'_>, df: &DataFrame) -> PyResult<PyObject> {
    let pandas = py.import_bound("pandas")?;
    let numpy = py.import_bound("numpy")?;
    let data_dict = PyDict::new_bound(py);

    for (col_name, col_data) in df.columns.iter().zip(df.data.iter()) {
        let py_col = arrow_array_to_python(py, col_data, &numpy, &pandas)?;
        data_dict.set_item(col_name, py_col)?;
    }

    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("data", data_dict)?;
    // When DataFrame has no columns, we still need the correct number of rows
    // so AnnData's obs/var dimension checks pass.
    if df.columns.is_empty() && df.n_rows > 0 {
        let range = py.import_bound("builtins")?.call_method1("range", (df.n_rows,))?;
        kwargs.set_item("index", range)?;
    }
    let result = pandas.getattr("DataFrame")?.call((), Some(&kwargs))?;
    Ok(result.into())
}

/// 将 Arrow ArrayRef 转换为 Python 对象
fn arrow_array_to_python(
    py: Python<'_>,
    array: &ArrayRef,
    numpy: &Bound<'_, PyModule>,
    pandas: &Bound<'_, PyModule>,
) -> PyResult<PyObject> {
    match array.data_type() {
        DataType::Int64 => {
            let arr = array.as_any().downcast_ref::<Int64Array>().unwrap();
            let values: Vec<Option<i64>> = (0..arr.len())
                .map(|i| if arr.is_null(i) { None } else { Some(arr.value(i)) })
                .collect();
            Ok(values.into_py(py))
        }
        DataType::Int32 => {
            let arr = array.as_any().downcast_ref::<Int32Array>().unwrap();
            let values: Vec<Option<i64>> = (0..arr.len())
                .map(|i| if arr.is_null(i) { None } else { Some(arr.value(i) as i64) })
                .collect();
            Ok(values.into_py(py))
        }
        DataType::Float64 => {
            let arr = array.as_any().downcast_ref::<Float64Array>().unwrap();
            let values: Vec<f64> = arr.values().iter().copied().collect();
            let np_arr = numpy.call_method1("array", (values,))?;
            Ok(np_arr.into())
        }
        DataType::Float32 => {
            let arr = array.as_any().downcast_ref::<Float32Array>().unwrap();
            let values: Vec<f64> = arr.values().iter().map(|&v| v as f64).collect();
            let np_arr = numpy.call_method1("array", (values,))?;
            Ok(np_arr.into())
        }
        DataType::Boolean => {
            let arr = array.as_any().downcast_ref::<BooleanArray>().unwrap();
            let values: Vec<Option<bool>> = (0..arr.len())
                .map(|i| if arr.is_null(i) { None } else { Some(arr.value(i)) })
                .collect();
            Ok(values.into_py(py))
        }
        DataType::Utf8 => {
            let arr = array.as_any().downcast_ref::<StringArray>().unwrap();
            let values: Vec<Option<&str>> = (0..arr.len())
                .map(|i| if arr.is_null(i) { None } else { Some(arr.value(i)) })
                .collect();
            Ok(values.into_py(py))
        }
        DataType::LargeUtf8 => {
            let arr = array.as_any().downcast_ref::<LargeStringArray>().unwrap();
            let values: Vec<Option<&str>> = (0..arr.len())
                .map(|i| if arr.is_null(i) { None } else { Some(arr.value(i)) })
                .collect();
            Ok(values.into_py(py))
        }
        DataType::Dictionary(_, _) => {
            arrow_dict_to_categorical(py, array, pandas)
        }
        _ => {
            // Fallback: convert each element to string
            let values: Vec<String> = (0..array.len())
                .map(|_| "unsupported".to_string())
                .collect();
            Ok(values.into_py(py))
        }
    }
}

/// 将 Arrow DictionaryArray 转换为 pandas Categorical
fn arrow_dict_to_categorical(
    py: Python<'_>,
    array: &ArrayRef,
    pandas: &Bound<'_, PyModule>,
) -> PyResult<PyObject> {
    // Try Int32 key (most common)
    if let Some(dict_arr) = array.as_any().downcast_ref::<DictionaryArray<Int32Type>>() {
        return dict_array_to_categorical::<Int32Type>(py, dict_arr, pandas);
    }
    // Try Int8 key
    if let Some(dict_arr) = array.as_any().downcast_ref::<DictionaryArray<Int8Type>>() {
        return dict_array_to_categorical::<Int8Type>(py, dict_arr, pandas);
    }
    Err(pyo3::exceptions::PyValueError::new_err(
        "Unsupported dictionary key type",
    ))
}

fn dict_array_to_categorical<K: arrow::datatypes::ArrowDictionaryKeyType>(
    py: Python<'_>,
    dict_arr: &DictionaryArray<K>,
    pandas: &Bound<'_, PyModule>,
) -> PyResult<PyObject>
where
    K::Native: Into<i64>,
{
    let keys = dict_arr.keys();
    let values = dict_arr.values();

    let categories: Vec<String> = if let Some(str_arr) = values.as_any().downcast_ref::<StringArray>() {
        (0..str_arr.len())
            .map(|i| if str_arr.is_null(i) { String::new() } else { str_arr.value(i).to_string() })
            .collect()
    } else {
        (0..values.len()).map(|i| format!("{}", i)).collect()
    };

    let codes: Vec<i64> = (0..keys.len())
        .map(|i| {
            if keys.is_null(i) { -1i64 } else { keys.value(i).into() }
        })
        .collect();

    let py_codes = PyList::new_bound(py, &codes);
    let py_categories = PyList::new_bound(py, &categories);

    let categorical_type = pandas.getattr("Categorical")?;
    let result = categorical_type.call_method1("from_codes", (py_codes, py_categories))?;
    Ok(result.into())
}

/// 将 Embedding 转换为 numpy.ndarray
pub fn embedding_to_numpy(py: Python<'_>, emb: &Embedding) -> PyResult<PyObject> {
    let numpy = py.import_bound("numpy")?;
    let flat = numpy.call_method1("array", (emb.data.clone(),))?;
    let shape = PyTuple::new_bound(py, &[emb.n_rows as i64, emb.n_cols as i64]);
    let reshaped = flat.call_method1("reshape", (shape,))?;
    Ok(reshaped.into())
}

/// 将 HashMap<String, UnstructuredValue> 转换为 Python dict
pub fn unstructured_map_to_dict(
    py: Python<'_>,
    uns: &HashMap<String, UnstructuredValue>,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    for (key, value) in uns {
        dict.set_item(key, unstructured_to_python(py, value)?)?;
    }
    Ok(dict.into())
}

fn unstructured_to_python(py: Python<'_>, value: &UnstructuredValue) -> PyResult<PyObject> {
    match value {
        UnstructuredValue::String(s) => Ok(s.into_py(py)),
        UnstructuredValue::Integer(i) => Ok(i.into_py(py)),
        UnstructuredValue::Float(f) => Ok(f.into_py(py)),
        UnstructuredValue::Boolean(b) => Ok(b.into_py(py)),
        UnstructuredValue::Null => Ok(py.None()),
        UnstructuredValue::Array(arr) => {
            let py_list = PyList::empty_bound(py);
            for item in arr {
                py_list.append(unstructured_to_python(py, item)?)?;
            }
            Ok(py_list.into())
        }
        UnstructuredValue::Dict(dict) => {
            let py_dict = PyDict::new_bound(py);
            for (k, v) in dict {
                py_dict.set_item(k, unstructured_to_python(py, v)?)?;
            }
            Ok(py_dict.into())
        }
    }
}
