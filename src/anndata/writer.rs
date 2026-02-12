//! AnnData (.h5ad) 文件写入器
//!
//! 将 IR 数据写入 HDF5 格式的 AnnData 文件。

use super::Result;
use crate::ir::{SingleCellData, ExpressionMatrix, SparseMatrixCSR, DenseMatrix, DataFrame, Embedding};
use hdf5::File;
use std::path::Path;
use std::collections::HashMap;
use arrow::array::{Array, AsArray};
use arrow::datatypes::DataType;

/// 将 IR 数据写入 .h5ad 文件
///
/// # 参数
/// - `data`: IR 数据
/// - `path`: 输出 .h5ad 文件路径
///
/// # 返回
/// - `Ok(())`: 成功写入
/// - `Err(AnnDataError)`: 写入失败
///
/// # 示例
/// ```no_run
/// use crosscell::anndata::write_h5ad;
/// use crosscell::ir::SingleCellData;
///
/// // let data: SingleCellData = ...;
/// // write_h5ad(&data, "output.h5ad").unwrap();
/// ```
pub fn write_h5ad<P: AsRef<Path>>(data: &SingleCellData, path: P) -> Result<()> {
    let file = File::create(path.as_ref())?;
    
    write_expression_matrix(&file, &data.expression)?;
    write_metadata_dataframe(&file, "obs", &data.cell_metadata)?;
    write_metadata_dataframe(&file, "var", &data.gene_metadata)?;
    
    // 分离 layers：同维度 → /layers，不同维度（多 Assay）→ /obsm
    let (same_dim_layers, diff_dim_layers) = if let Some(ref layers) = data.layers {
        split_layers_by_dimension(layers, data.metadata.n_cells, data.metadata.n_genes)
    } else {
        (HashMap::new(), HashMap::new())
    };
    
    // 写入降维嵌入 /obsm
    if let Some(ref embeddings) = data.embeddings {
        write_embeddings(&file, embeddings, data.metadata.n_cells)?;
    }
    
    // 写入不同维度的 assay 到 /obsm（如 ADT）
    if !diff_dim_layers.is_empty() {
        write_multi_assay_to_obsm(&file, &diff_dim_layers, data.metadata.n_cells)?;
    }
    
    // 写入同维度 layers 到 /layers
    if !same_dim_layers.is_empty() {
        write_layers(&file, &same_dim_layers, data.metadata.n_cells, data.metadata.n_genes)?;
    }
    
    // 写入多 Assay 元信息到 /uns/crosscell_multi_assay
    if !diff_dim_layers.is_empty() {
        write_multi_assay_metadata(&file, &diff_dim_layers)?;
    }
    
    // 写入 obsp
    if let Some(ref cell_pairwise) = data.cell_pairwise {
        if !cell_pairwise.is_empty() {
            write_pairwise_matrices(&file, "obsp", cell_pairwise)?;
        }
    }
    
    // 写入 varp
    if let Some(ref gene_pairwise) = data.gene_pairwise {
        if !gene_pairwise.is_empty() {
            write_pairwise_matrices(&file, "varp", gene_pairwise)?;
        }
    }
    
    // 写入空间数据
    if let Some(ref spatial) = data.spatial {
        write_spatial_data(&file, spatial, data.metadata.n_cells)?;
    }
    
    Ok(())
}


/// 写入表达矩阵到 HDF5 文件 (/X)
///
/// 支持两种格式：
/// 1. 稀疏矩阵（CSR 格式）：/X 是一个 Group，包含 data, indices, indptr
/// 2. 稠密矩阵：/X 是一个 Dataset
fn write_expression_matrix(file: &File, matrix: &ExpressionMatrix) -> Result<()> {
    match matrix {
        ExpressionMatrix::SparseCSR(csr) => write_sparse_csr_matrix(file, "X", csr),
        ExpressionMatrix::SparseCSC(csc) => {
            // 将 CSC 转换为 CSR 后写入（AnnData 标准格式是 CSR）
            let csr = crate::sparse::convert::csc_to_csr(csc);
            write_sparse_csr_matrix(file, "X", &csr)
        }
        ExpressionMatrix::Dense(dense) => write_dense_matrix(file, "X", dense),
        ExpressionMatrix::Lazy(lazy) => {
            // 延迟加载矩阵：尝试从缓存获取，否则返回错误
            if let Some(cached) = lazy.get_cached() {
                write_expression_matrix(file, &cached)
            } else {
                Err(super::AnnDataError::InvalidFormat(
                    "Cannot write LazyMatrix without loading data first. Call load() before writing.".to_string()
                ))
            }
        }
    }
}

/// 写入稀疏 CSR 矩阵
///
/// HDF5 结构：
/// - /X: Group
/// - /X/data: 非零元素值 (float64)
/// - /X/indices: 列索引 (int32)
/// - /X/indptr: 行指针 (int32)
/// - /X 的属性：shape, encoding-type, encoding-version
fn write_sparse_csr_matrix(file: &File, name: &str, csr: &SparseMatrixCSR) -> Result<()> {
    // 创建 Group
    let group = file.create_group(name)?;
    
    // 写入 shape 属性
    let shape = vec![csr.n_rows as i64, csr.n_cols as i64];
    group.new_attr::<i64>()
        .shape([2])
        .create("shape")?
        .write(&shape)?;
    
    // 写入 encoding-type 属性
    group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-type")?
        .write_scalar(&"csr_matrix".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    // 写入 encoding-version 属性
    group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-version")?
        .write_scalar(&"0.1.0".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    // 写入 data (非零元素值)
    let data_dataset = group.new_dataset::<f64>()
        .shape([csr.data.len()])
        .create("data")?;
    data_dataset.write(&csr.data)?;
    
    // 写入 indices (列索引，转换为 i32)
    let indices_i32: Vec<i32> = csr.indices.iter().map(|&x| x as i32).collect();
    let indices_dataset = group.new_dataset::<i32>()
        .shape([indices_i32.len()])
        .create("indices")?;
    indices_dataset.write(&indices_i32)?;
    
    // 写入 indptr (行指针，转换为 i32)
    let indptr_i32: Vec<i32> = csr.indptr.iter().map(|&x| x as i32).collect();
    let indptr_dataset = group.new_dataset::<i32>()
        .shape([indptr_i32.len()])
        .create("indptr")?;
    indptr_dataset.write(&indptr_i32)?;
    
    Ok(())
}

/// 写入稠密矩阵
///
/// HDF5 结构：
/// - /X: Dataset (n_cells × n_genes)
fn write_dense_matrix(file: &File, name: &str, dense: &DenseMatrix) -> Result<()> {
    // 创建 Dataset
    let dataset = file.new_dataset::<f64>()
        .shape([dense.n_rows, dense.n_cols])
        .create(name)?;
    
    // 写入数据（行优先）
    // HDF5 expects data in the shape specified, so we need to write it as a 2D array
    // The data is already in row-major order, so we can write it directly
    // but we need to use write_raw to specify the shape
    dataset.write_raw(&dense.data)?;
    
    Ok(())
}

/// 写入元数据 DataFrame (/obs 或 /var)
///
/// HDF5 结构：
/// - /obs 或 /var: Group，包含多个列
/// - 每列是一个 Dataset（普通列）或 Group（Categorical 列）
/// - Categorical 列（新版本格式 0.2.0）：
///   - /obs/{column_name}: Group
///   - /obs/{column_name}/categories: 类别标签
///   - /obs/{column_name}/codes: 整数编码
fn write_metadata_dataframe(file: &File, group_name: &str, df: &DataFrame) -> Result<()> {
    // 如果 DataFrame 为空，创建空 Group
    if df.n_cols() == 0 {
        file.create_group(group_name)?;
        return Ok(());
    }
    
    // 创建主 Group
    let group = file.create_group(group_name)?;
    
    // 写入每一列
    for (col_name, array) in df.columns.iter().zip(df.data.iter()) {
        match array.data_type() {
            DataType::Float64 => {
                write_float64_column(&group, col_name, array)?;
            }
            DataType::Int64 => {
                write_int64_column(&group, col_name, array)?;
            }
            DataType::Int32 => {
                write_int32_column(&group, col_name, array)?;
            }
            DataType::Boolean => {
                write_boolean_column(&group, col_name, array)?;
            }
            DataType::Utf8 => {
                write_string_column(&group, col_name, array)?;
            }
            DataType::Dictionary(_, _) => {
                // Categorical 列（新版本格式，不需要 categories_group）
                write_categorical_column_v2(&group, col_name, array)?;
            }
            _ => {
                return Err(super::AnnDataError::UnsupportedType(format!(
                    "Unsupported column type: {:?}",
                    array.data_type()
                )));
            }
        }
    }
    
    Ok(())
}

/// 写入 Float64 列
fn write_float64_column(group: &hdf5::Group, col_name: &str, array: &arrow::array::ArrayRef) -> Result<()> {
    let float_array = array.as_primitive::<arrow::datatypes::Float64Type>();
    let values: Vec<f64> = float_array.values().to_vec();
    
    let dataset = group.new_dataset::<f64>()
        .shape([values.len()])
        .create(col_name)?;
    dataset.write(&values)?;
    
    Ok(())
}

/// 写入 Int64 列
fn write_int64_column(group: &hdf5::Group, col_name: &str, array: &arrow::array::ArrayRef) -> Result<()> {
    let int_array = array.as_primitive::<arrow::datatypes::Int64Type>();
    let values: Vec<i64> = int_array.values().to_vec();
    
    let dataset = group.new_dataset::<i64>()
        .shape([values.len()])
        .create(col_name)?;
    dataset.write(&values)?;
    
    Ok(())
}

/// 写入 Int32 列
fn write_int32_column(group: &hdf5::Group, col_name: &str, array: &arrow::array::ArrayRef) -> Result<()> {
    let int_array = array.as_primitive::<arrow::datatypes::Int32Type>();
    let values: Vec<i32> = int_array.values().to_vec();
    
    let dataset = group.new_dataset::<i32>()
        .shape([values.len()])
        .create(col_name)?;
    dataset.write(&values)?;
    
    Ok(())
}

/// 写入 Boolean 列
fn write_boolean_column(group: &hdf5::Group, col_name: &str, array: &arrow::array::ArrayRef) -> Result<()> {
    let bool_array = array.as_boolean();
    let values: Vec<bool> = (0..bool_array.len())
        .map(|i| bool_array.value(i))
        .collect();
    
    let dataset = group.new_dataset::<bool>()
        .shape([values.len()])
        .create(col_name)?;
    dataset.write(&values)?;
    
    Ok(())
}

/// 写入 String 列
fn write_string_column(group: &hdf5::Group, col_name: &str, array: &arrow::array::ArrayRef) -> Result<()> {
    let string_array = array.as_string::<i32>();
    
    // 将字符串转换为 VarLenUnicode
    let values: Vec<hdf5::types::VarLenUnicode> = (0..string_array.len())
        .map(|i| {
            let s = string_array.value(i);
            s.parse::<hdf5::types::VarLenUnicode>()
                .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    
    let dataset = group.new_dataset::<hdf5::types::VarLenUnicode>()
        .shape([values.len()])
        .create(col_name)?;
    dataset.write(&values)?;
    
    Ok(())
}

/// 写入 Categorical 列（新版本格式 0.2.0）
///
/// HDF5 结构：
/// - /obs/{column_name}: Group（不是 dataset）
/// - /obs/{column_name}/categories: 类别标签 (string array)
/// - /obs/{column_name}/codes: 整数编码 (int32)
/// - 属性：encoding-type="categorical", encoding-version="0.2.0", ordered=False
fn write_categorical_column_v2(
    group: &hdf5::Group,
    col_name: &str,
    array: &arrow::array::ArrayRef,
) -> Result<()> {
    use arrow::array::DictionaryArray;
    use arrow::datatypes::Int32Type;
    
    // 转换为 DictionaryArray
    let dict_array = array
        .as_any()
        .downcast_ref::<DictionaryArray<Int32Type>>()
        .ok_or_else(|| {
            super::AnnDataError::InvalidFormat("Expected DictionaryArray<Int32Type>".to_string())
        })?;
    
    // 提取整数编码
    let keys = dict_array.keys();
    let codes: Vec<i32> = keys.values().to_vec();
    
    // 提取类别标签
    let values = dict_array.values();
    let string_values = values
        .as_any()
        .downcast_ref::<arrow::array::StringArray>()
        .ok_or_else(|| {
            super::AnnDataError::InvalidFormat("Expected StringArray for dictionary values".to_string())
        })?;
    
    let categories: Vec<hdf5::types::VarLenUnicode> = (0..string_values.len())
        .map(|i| {
            let s = string_values.value(i);
            s.parse::<hdf5::types::VarLenUnicode>()
                .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    
    // 创建 categorical 列的 Group（新版本格式）
    let cat_col_group = group.create_group(col_name)?;
    
    // 写入属性
    cat_col_group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-type")?
        .write_scalar(&"categorical".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    cat_col_group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-version")?
        .write_scalar(&"0.2.0".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    // ordered 属性（默认 False）
    cat_col_group.new_attr::<bool>()
        .create("ordered")?
        .write_scalar(&false)?;
    
    // 写入 codes
    let codes_dataset = cat_col_group.new_dataset::<i32>()
        .shape([codes.len()])
        .create("codes")?;
    codes_dataset.write(&codes)?;
    
    // 写入 categories
    let cat_dataset = cat_col_group.new_dataset::<hdf5::types::VarLenUnicode>()
        .shape([categories.len()])
        .create("categories")?;
    cat_dataset.write(&categories)?;
    
    Ok(())
}

/// 写入降维嵌入 (/obsm)
///
/// HDF5 结构：
/// - /obsm: Group，包含多个嵌入
/// - 每个嵌入是一个 Dataset (n_cells × n_components)
fn write_embeddings(
    file: &File,
    embeddings: &HashMap<String, Embedding>,
    expected_cells: usize,
) -> Result<()> {
    // 创建 /obsm Group
    let obsm_group = file.create_group("obsm")?;
    
    // 写入每个嵌入
    for (name, embedding) in embeddings {
        // 验证行数
        if embedding.n_rows != expected_cells {
            return Err(super::AnnDataError::DimensionMismatch(format!(
                "Embedding '{}' has {} rows, expected {}",
                name, embedding.n_rows, expected_cells
            )));
        }
        
        // 创建 Dataset (n_rows × n_cols)
        let dataset = obsm_group.new_dataset::<f64>()
            .shape([embedding.n_rows, embedding.n_cols])
            .create(name.as_str())?;
        
        // 写入数据（行优先）
        dataset.write_raw(&embedding.data)?;
    }
    
    Ok(())
}

/// 写入 layers (/layers)
///
/// HDF5 结构：
/// - /layers: Group，包含多个层
/// - 每层是一个表达矩阵（稀疏或稠密）
fn write_layers(
    file: &File,
    layers: &HashMap<String, ExpressionMatrix>,
    expected_cells: usize,
    expected_genes: usize,
) -> Result<()> {
    // 创建 /layers Group
    let layers_group = file.create_group("layers")?;
    
    // 写入每个 layer
    for (name, matrix) in layers {
        // 验证维度
        let (n_rows, n_cols) = matrix.shape();
        if n_rows != expected_cells || n_cols != expected_genes {
            return Err(super::AnnDataError::DimensionMismatch(format!(
                "Layer '{}' has shape ({}, {}), expected ({}, {})",
                name, n_rows, n_cols, expected_cells, expected_genes
            )));
        }
        
        // 根据矩阵类型写入
        match matrix {
            ExpressionMatrix::SparseCSR(csr) => {
                // 创建子 Group 并写入稀疏矩阵
                let layer_group = layers_group.create_group(name.as_str())?;
                write_sparse_csr_to_group(&layer_group, csr)?;
            }
            ExpressionMatrix::SparseCSC(csc) => {
                // 转换为 CSR 后写入
                let csr = crate::sparse::convert::csc_to_csr(csc);
                let layer_group = layers_group.create_group(name.as_str())?;
                write_sparse_csr_to_group(&layer_group, &csr)?;
            }
            ExpressionMatrix::Dense(dense) => {
                // 创建 Dataset 并写入稠密矩阵
                let dataset = layers_group.new_dataset::<f64>()
                    .shape([dense.n_rows, dense.n_cols])
                    .create(name.as_str())?;
                dataset.write_raw(&dense.data)?;
            }
            ExpressionMatrix::Lazy(lazy) => {
                // 延迟加载矩阵：尝试从缓存获取
                if let Some(cached) = lazy.get_cached() {
                    // 递归处理缓存的矩阵
                    let mut temp_layers = HashMap::new();
                    temp_layers.insert(name.clone(), cached);
                    // 注意：这里我们直接处理缓存的矩阵
                    match temp_layers.get(name).unwrap() {
                        ExpressionMatrix::SparseCSR(csr) => {
                            let layer_group = layers_group.create_group(name.as_str())?;
                            write_sparse_csr_to_group(&layer_group, csr)?;
                        }
                        ExpressionMatrix::SparseCSC(csc) => {
                            let csr = crate::sparse::convert::csc_to_csr(csc);
                            let layer_group = layers_group.create_group(name.as_str())?;
                            write_sparse_csr_to_group(&layer_group, &csr)?;
                        }
                        ExpressionMatrix::Dense(dense) => {
                            let dataset = layers_group.new_dataset::<f64>()
                                .shape([dense.n_rows, dense.n_cols])
                                .create(name.as_str())?;
                            dataset.write_raw(&dense.data)?;
                        }
                        ExpressionMatrix::Lazy(_) => {
                            return Err(super::AnnDataError::DimensionMismatch(
                                "Nested LazyMatrix not supported".to_string()
                            ));
                        }
                    }
                } else {
                    return Err(super::AnnDataError::DimensionMismatch(format!(
                        "Layer '{}' is a LazyMatrix without cached data. Load data first.",
                        name
                    )));
                }
            }
        }
    }
    
    Ok(())
}

/// 按维度分组 layers：与主矩阵同维度的放入 /layers，不同维度的放入 /obsm
///
/// 多 Assay 数据集（如 CITE-seq 的 ADT assay）的基因数与主 RNA assay 不同，
/// 无法放入 /layers（要求与 /X 同维度），改为转成 dense 写入 /obsm。
fn split_layers_by_dimension(
    layers: &HashMap<String, ExpressionMatrix>,
    n_cells: usize,
    n_genes: usize,
) -> (HashMap<String, ExpressionMatrix>, HashMap<String, ExpressionMatrix>) {
    let mut same_dim = HashMap::new();
    let mut diff_dim = HashMap::new();
    
    for (name, matrix) in layers {
        let (rows, cols) = matrix.shape();
        if rows == n_cells && cols == n_genes {
            same_dim.insert(name.clone(), matrix.clone());
        } else if rows == n_cells {
            // 细胞数相同但基因数不同 → 多 Assay，放入 /obsm
            diff_dim.insert(name.clone(), matrix.clone());
        } else {
            // 维度完全不匹配，仍然尝试放入 diff_dim（后续写入时会报错）
            diff_dim.insert(name.clone(), matrix.clone());
        }
    }
    
    (same_dim, diff_dim)
}

/// 将不同维度的 ExpressionMatrix（多 Assay）转为 dense 写入 /obsm
///
/// 例如 CITE-seq 的 ADT assay (n_cells × n_adt_features) 写入 /obsm/X_ADT
fn write_multi_assay_to_obsm(
    file: &File,
    diff_dim_layers: &HashMap<String, ExpressionMatrix>,
    expected_cells: usize,
) -> Result<()> {
    // 确保 /obsm group 存在
    let obsm_group = if file.link_exists("obsm") {
        file.group("obsm")?
    } else {
        file.create_group("obsm")?
    };
    
    for (name, matrix) in diff_dim_layers {
        let (n_rows, _n_cols) = matrix.shape();
        if n_rows != expected_cells {
            return Err(super::AnnDataError::DimensionMismatch(format!(
                "Multi-assay layer '{}' has {} rows, expected {} cells",
                name, n_rows, expected_cells
            )));
        }
        
        // 命名约定：X_{assay_name}
        let obsm_key = format!("X_{}", name);
        
        // 转为 dense 写入
        let dense = expression_to_dense(matrix)?;
        let dataset = obsm_group.new_dataset::<f64>()
            .shape([dense.n_rows, dense.n_cols])
            .create(obsm_key.as_str())?;
        dataset.write_raw(&dense.data)?;
    }
    
    Ok(())
}

/// 将 ExpressionMatrix 转为 DenseMatrix
fn expression_to_dense(matrix: &ExpressionMatrix) -> Result<DenseMatrix> {
    match matrix {
        ExpressionMatrix::Dense(dense) => Ok(dense.clone()),
        ExpressionMatrix::SparseCSR(csr) => {
            Ok(crate::sparse::convert::csr_to_dense(csr))
        }
        ExpressionMatrix::SparseCSC(csc) => {
            Ok(crate::sparse::convert::csc_to_dense(csc))
        }
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                expression_to_dense(&cached)
            } else {
                Err(super::AnnDataError::InvalidFormat(
                    "Cannot convert LazyMatrix to dense without cached data".to_string()
                ))
            }
        }
    }
}

/// 在 /uns/crosscell_multi_assay 记录多 Assay 元信息
///
/// 记录每个额外 assay 的名称、维度、在 obsm 中的 key，
/// 便于反向转换时恢复为独立 Assay。
fn write_multi_assay_metadata(
    file: &File,
    diff_dim_layers: &HashMap<String, ExpressionMatrix>,
) -> Result<()> {
    // 确保 /uns group 存在
    let uns_group = if file.link_exists("uns") {
        file.group("uns")?
    } else {
        file.create_group("uns")?
    };
    
    // 创建 /uns/crosscell_multi_assay group
    let meta_group = uns_group.create_group("crosscell_multi_assay")?;
    
    // 记录 assay 名称列表
    let assay_names: Vec<hdf5::types::VarLenUnicode> = diff_dim_layers.keys()
        .map(|name| name.parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string())))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    
    let names_ds = meta_group.new_dataset::<hdf5::types::VarLenUnicode>()
        .shape([assay_names.len()])
        .create("assay_names")?;
    names_ds.write(&assay_names)?;
    
    // 为每个 assay 记录维度信息
    for (name, matrix) in diff_dim_layers {
        let assay_group = meta_group.create_group(name.as_str())?;
        let (n_rows, n_cols) = matrix.shape();
        
        // obsm_key
        let obsm_key = format!("X_{}", name);
        assay_group.new_attr::<hdf5::types::VarLenUnicode>()
            .create("obsm_key")?
            .write_scalar(&obsm_key.parse::<hdf5::types::VarLenUnicode>()
                .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
        
        // shape
        let shape = vec![n_rows as i64, n_cols as i64];
        assay_group.new_attr::<i64>()
            .shape([2])
            .create("shape")?
            .write(&shape)?;
    }
    
    Ok(())
}

/// 写入成对矩阵 (/obsp 或 /varp)
///
/// HDF5 结构：
/// - /obsp 或 /varp: Group，包含多个成对矩阵
/// - 每个矩阵是一个方阵（稀疏或稠密）
fn write_pairwise_matrices(
    file: &File,
    group_name: &str,
    pairwise_matrices: &HashMap<String, crate::ir::PairwiseMatrix>,
) -> Result<()> {
    // 创建 Group
    let pw_group = file.create_group(group_name)?;
    
    // 写入每个成对矩阵
    for (name, pw_matrix) in pairwise_matrices {
        // 验证是方阵
        let (n_rows, n_cols) = pw_matrix.matrix.shape();
        if n_rows != n_cols {
            return Err(super::AnnDataError::InvalidFormat(format!(
                "Pairwise matrix '{}' is not square: {} × {}",
                name, n_rows, n_cols
            )));
        }
        
        // 根据矩阵类型写入
        match &pw_matrix.matrix {
            ExpressionMatrix::SparseCSR(csr) => {
                // 创建子 Group 并写入稀疏矩阵
                let matrix_group = pw_group.create_group(name.as_str())?;
                write_sparse_csr_to_group(&matrix_group, csr)?;
            }
            ExpressionMatrix::SparseCSC(csc) => {
                // 转换为 CSR 后写入
                let csr = crate::sparse::convert::csc_to_csr(csc);
                let matrix_group = pw_group.create_group(name.as_str())?;
                write_sparse_csr_to_group(&matrix_group, &csr)?;
            }
            ExpressionMatrix::Dense(dense) => {
                // 创建 Dataset 并写入稠密矩阵
                let dataset = pw_group.new_dataset::<f64>()
                    .shape([dense.n_rows, dense.n_cols])
                    .create(name.as_str())?;
                dataset.write_raw(&dense.data)?;
            }
            ExpressionMatrix::Lazy(lazy) => {
                // 延迟加载矩阵：尝试从缓存获取
                if let Some(cached) = lazy.get_cached() {
                    match &cached {
                        ExpressionMatrix::SparseCSR(csr) => {
                            let matrix_group = pw_group.create_group(name.as_str())?;
                            write_sparse_csr_to_group(&matrix_group, csr)?;
                        }
                        ExpressionMatrix::SparseCSC(csc) => {
                            let csr = crate::sparse::convert::csc_to_csr(csc);
                            let matrix_group = pw_group.create_group(name.as_str())?;
                            write_sparse_csr_to_group(&matrix_group, &csr)?;
                        }
                        ExpressionMatrix::Dense(dense) => {
                            let dataset = pw_group.new_dataset::<f64>()
                                .shape([dense.n_rows, dense.n_cols])
                                .create(name.as_str())?;
                            dataset.write_raw(&dense.data)?;
                        }
                        ExpressionMatrix::Lazy(_) => {
                            return Err(super::AnnDataError::InvalidFormat(
                                "Nested LazyMatrix not supported".to_string()
                            ));
                        }
                    }
                } else {
                    return Err(super::AnnDataError::InvalidFormat(format!(
                        "Pairwise matrix '{}' is a LazyMatrix without cached data. Load data first.",
                        name
                    )));
                }
            }
        }
    }
    
    Ok(())
}

/// 辅助函数：将 CSR 矩阵写入已存在的 Group
///
/// 这个函数用于 layers 和 pairwise matrices，避免代码重复
fn write_sparse_csr_to_group(group: &hdf5::Group, csr: &SparseMatrixCSR) -> Result<()> {
    // 写入 shape 属性
    let shape = vec![csr.n_rows as i64, csr.n_cols as i64];
    group.new_attr::<i64>()
        .shape([2])
        .create("shape")?
        .write(&shape)?;
    
    // 写入 encoding-type 属性
    group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-type")?
        .write_scalar(&"csr_matrix".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    // 写入 encoding-version 属性
    group.new_attr::<hdf5::types::VarLenUnicode>()
        .create("encoding-version")?
        .write_scalar(&"0.1.0".parse::<hdf5::types::VarLenUnicode>()
            .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?)?;
    
    // 写入 data (非零元素值)
    let data_dataset = group.new_dataset::<f64>()
        .shape([csr.data.len()])
        .create("data")?;
    data_dataset.write(&csr.data)?;
    
    // 写入 indices (列索引，转换为 i32)
    let indices_i32: Vec<i32> = csr.indices.iter().map(|&x| x as i32).collect();
    let indices_dataset = group.new_dataset::<i32>()
        .shape([indices_i32.len()])
        .create("indices")?;
    indices_dataset.write(&indices_i32)?;
    
    // 写入 indptr (行指针，转换为 i32)
    let indptr_i32: Vec<i32> = csr.indptr.iter().map(|&x| x as i32).collect();
    let indptr_dataset = group.new_dataset::<i32>()
        .shape([indptr_i32.len()])
        .create("indptr")?;
    indptr_dataset.write(&indptr_i32)?;
    
    Ok(())
}

/// 写入空间数据 (/obsm['spatial'] 和 /uns['spatial'])
///
/// AnnData 空间数据结构：
/// - /obsm/spatial: 空间坐标矩阵 (n_cells × 2 或 n_cells × 3)
/// - /uns/spatial: Group，包含图像和缩放因子
///   - /uns/spatial/{library_id}/images: Group，包含多分辨率图像
///   - /uns/spatial/{library_id}/scalefactors: 缩放因子
fn write_spatial_data(
    file: &File,
    spatial: &crate::ir::SpatialData,
    expected_cells: usize,
) -> Result<()> {
    // 验证细胞数
    if spatial.n_cells() != expected_cells {
        return Err(super::AnnDataError::DimensionMismatch(format!(
            "Spatial data has {} cells, expected {}",
            spatial.n_cells(),
            expected_cells
        )));
    }
    
    // 1. 写入空间坐标到 /obsm/spatial
    // 确保 /obsm group 存在
    let obsm_group = if file.link_exists("obsm") {
        file.group("obsm")?
    } else {
        file.create_group("obsm")?
    };
    
    // 写入坐标数据 - 检查是否已存在
    if !obsm_group.link_exists("spatial") {
        let spatial_dataset = obsm_group.new_dataset::<f64>()
            .shape([spatial.n_cells(), spatial.n_dims])
            .create("spatial")?;
        spatial_dataset.write_raw(&spatial.coordinates)?;
    }
    
    // 2. 写入图像和缩放因子到 /uns/spatial（如果存在）
    if spatial.images.is_some() || spatial.scale_factors.is_some() {
        // 确保 /uns group 存在
        let uns_group = if file.link_exists("uns") {
            file.group("uns")?
        } else {
            file.create_group("uns")?
        };
        
        // 创建 /uns/spatial group
        let spatial_group = if uns_group.link_exists("spatial") {
            uns_group.group("spatial")?
        } else {
            uns_group.create_group("spatial")?
        };
        
        // 使用默认 library_id "library_1"
        let library_id = "library_1";
        let library_group = if spatial_group.link_exists(library_id) {
            spatial_group.group(library_id)?
        } else {
            spatial_group.create_group(library_id)?
        };
        
        // 写入图像
        if let Some(ref images) = spatial.images {
            let images_group = if library_group.link_exists("images") {
                library_group.group("images")?
            } else {
                library_group.create_group("images")?
            };
            
            for img in images {
                // 验证图像数据大小
                let expected_size = img.width * img.height * 3; // RGB
                if img.data.len() != expected_size {
                    return Err(super::AnnDataError::InvalidFormat(format!(
                        "Image '{}' data size {} doesn't match dimensions {}×{}×3 = {}",
                        img.name, img.data.len(), img.width, img.height, expected_size
                    )));
                }
                
                // 创建 3D dataset (height, width, channels)
                let img_dataset = images_group.new_dataset::<u8>()
                    .shape([img.height, img.width, 3])
                    .create(img.name.as_str())?;
                
                // 写入图像数据（已经是行优先的 RGB 格式）
                img_dataset.write_raw(&img.data)?;
            }
        }
        
        // 写入缩放因子
        if let Some(ref scale_factors) = spatial.scale_factors {
            let scalefactors_group = if library_group.link_exists("scalefactors") {
                library_group.group("scalefactors")?
            } else {
                library_group.create_group("scalefactors")?
            };
            
            for (factor_name, &factor_value) in scale_factors {
                // 创建标量 dataset
                let factor_dataset = scalefactors_group.new_dataset::<f64>()
                    .create(factor_name.as_str())?;
                factor_dataset.write_scalar(&factor_value)?;
            }
        }
    }
    
    Ok(())
}
