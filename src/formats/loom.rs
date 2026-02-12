//! Loom (.loom) 格式转换器
//!
//! Loom 是一种基于 HDF5 的单细胞数据格式，由 Linnarsson Lab 开发。
//!
//! ## Loom 文件结构
//!
//! ```text
//! /matrix          - 主表达矩阵 (genes × cells)，注意是转置的
//! /layers/         - 额外的矩阵层 (spliced, unspliced, ambiguous)
//! /row_attrs/      - 基因属性 (相当于 AnnData 的 var)
//! /col_attrs/      - 细胞属性 (相当于 AnnData 的 obs)
//! /col_graphs/     - 细胞图 (KNN 等)
//! /row_graphs/     - 基因图
//! ```
//!
//! ## 参考
//!
//! - Loom 文件规范: https://linnarssonlab.org/loompy/format/

use std::path::Path;
use std::collections::HashMap;
use std::sync::Arc;

use hdf5::{File as H5File, Group};
use arrow::array::{ArrayRef, Float64Array, Int64Array, StringArray, BooleanArray};

use crate::error::CrossCellError;
use crate::ir::{
    SingleCellData, ExpressionMatrix, DenseMatrix, SparseMatrixCSR,
    DataFrame, DatasetMetadata, Embedding, PairwiseMatrix,
};

use super::FormatConverter;

/// Loom 格式转换器
pub struct LoomConverter;

impl FormatConverter for LoomConverter {
    fn name(&self) -> &str {
        "loom"
    }
    
    fn display_name(&self) -> &str {
        "Loom"
    }
    
    fn extensions(&self) -> &[&str] {
        &[".loom"]
    }
    
    fn can_read(&self) -> bool {
        true
    }
    
    fn can_write(&self) -> bool {
        true
    }
    
    fn description(&self) -> &str {
        "Loom HDF5 format for scVelo/velocyto"
    }
    
    fn read(&self, path: &Path) -> Result<SingleCellData, CrossCellError> {
        read_loom(path)
    }
    
    fn write(&self, data: &SingleCellData, path: &Path) -> Result<(), CrossCellError> {
        write_loom(data, path)
    }
}

/// 读取 Loom 文件
fn read_loom(path: &Path) -> Result<SingleCellData, CrossCellError> {
    let file = H5File::open(path).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open loom file: {}", e))
    })?;
    
    // 读取主矩阵 (genes × cells)
    let matrix_ds = file.dataset("matrix").map_err(|e| {
        CrossCellError::MissingField { field: format!("/matrix: {}", e) }
    })?;
    
    let shape = matrix_ds.shape();
    if shape.len() != 2 {
        return Err(CrossCellError::InvalidDimensions(
            format!("Expected 2D matrix, got {}D", shape.len())
        ));
    }
    
    let n_genes = shape[0];
    let n_cells = shape[1];
    
    // 读取矩阵数据并转置 (genes × cells → cells × genes)
    let matrix_data: Vec<f64> = matrix_ds.read_raw().map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to read matrix data: {}", e))
    })?;
    
    // 转置矩阵
    let mut transposed = vec![0.0; n_cells * n_genes];
    for i in 0..n_genes {
        for j in 0..n_cells {
            transposed[j * n_genes + i] = matrix_data[i * n_cells + j];
        }
    }
    
    let expression = ExpressionMatrix::Dense(DenseMatrix {
        n_rows: n_cells,
        n_cols: n_genes,
        data: transposed,
    });
    
    // 读取细胞属性 (/col_attrs)
    let (cell_metadata, embeddings) = read_col_attrs(&file, n_cells)?;
    
    // 读取基因属性 (/row_attrs)
    let gene_metadata = read_row_attrs(&file, n_genes)?;
    
    // 读取 layers (/layers)
    let layers = read_layers(&file, n_cells, n_genes)?;
    
    // 读取图 (/col_graphs, /row_graphs)
    let (cell_pairwise, gene_pairwise) = read_graphs(&file, n_cells, n_genes)?;
    
    Ok(SingleCellData {
        expression,
        cell_metadata,
        gene_metadata,
        embeddings,
        layers,
        cell_pairwise,
        gene_pairwise,
        unstructured: None,
        spatial: None,
        metadata: DatasetMetadata::new(n_cells, n_genes, "loom".to_string()),
    })
}

/// 读取细胞属性 (/col_attrs)
/// 
/// 返回 (DataFrame, Option<HashMap<String, Embedding>>)
/// 嵌入坐标（如 _X, _Y）会被提取为 Embedding
fn read_col_attrs(file: &H5File, n_cells: usize) -> Result<(DataFrame, Option<HashMap<String, Embedding>>), CrossCellError> {
    if !file.link_exists("col_attrs") {
        return Ok((DataFrame::empty(n_cells), None));
    }
    
    let col_attrs = file.group("col_attrs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open /col_attrs: {}", e))
    })?;
    
    let member_names = col_attrs.member_names().map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to list /col_attrs members: {}", e))
    })?;
    
    let mut columns: Vec<String> = Vec::new();
    let mut data: Vec<ArrayRef> = Vec::new();
    let mut embedding_coords: HashMap<String, Vec<f64>> = HashMap::new();
    
    for name in member_names {
        // 检查是否是嵌入坐标（以 _ 开头，如 _X, _Y, _tSNE1, _tSNE2）
        if name.starts_with('_') && is_embedding_coord(&name) {
            if let Ok(ds) = col_attrs.dataset(&name) {
                if let Ok(values) = ds.read_1d::<f64>() {
                    embedding_coords.insert(name.clone(), values.to_vec());
                }
            }
            continue;
        }
        
        // 读取普通属性
        if let Ok(ds) = col_attrs.dataset(&name) {
            if let Some(array) = read_attribute_dataset(&ds, n_cells)? {
                columns.push(name);
                data.push(array);
            }
        }
    }
    
    // 构建 DataFrame
    let cell_metadata = if columns.is_empty() {
        DataFrame::empty(n_cells)
    } else {
        DataFrame::new(columns, data, n_cells).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create cell metadata DataFrame: {}", e))
        })?
    };
    
    // 构建嵌入
    let embeddings = build_embeddings_from_coords(embedding_coords, n_cells)?;
    
    Ok((cell_metadata, embeddings))
}


/// 检查属性名是否是嵌入坐标
fn is_embedding_coord(name: &str) -> bool {
    // 常见的嵌入坐标命名模式
    let patterns = ["_X", "_Y", "_Z", "_tSNE", "_UMAP", "_PCA", "_PC"];
    for pattern in patterns {
        if name.starts_with(pattern) || name.contains(pattern) {
            return true;
        }
    }
    false
}

/// 从坐标 HashMap 构建嵌入
fn build_embeddings_from_coords(
    coords: HashMap<String, Vec<f64>>,
    n_cells: usize,
) -> Result<Option<HashMap<String, Embedding>>, CrossCellError> {
    if coords.is_empty() {
        return Ok(None);
    }
    
    let mut embeddings: HashMap<String, Embedding> = HashMap::new();
    
    // 检查 UMAP 坐标 (_X, _Y)
    if coords.contains_key("_X") && coords.contains_key("_Y") {
        let x = coords.get("_X").unwrap();
        let y = coords.get("_Y").unwrap();
        
        if x.len() == n_cells && y.len() == n_cells {
            let mut data = Vec::with_capacity(n_cells * 2);
            for i in 0..n_cells {
                data.push(x[i]);
                data.push(y[i]);
            }
            
            let _n_dims = if coords.contains_key("_Z") {
                let z = coords.get("_Z").unwrap();
                if z.len() == n_cells {
                    // 3D 嵌入
                    let mut data_3d = Vec::with_capacity(n_cells * 3);
                    for i in 0..n_cells {
                        data_3d.push(x[i]);
                        data_3d.push(y[i]);
                        data_3d.push(z[i]);
                    }
                    embeddings.insert("X_umap".to_string(), Embedding {
                        name: "X_umap".to_string(),
                        data: data_3d,
                        n_rows: n_cells,
                        n_cols: 3,
                    });
                    3
                } else {
                    embeddings.insert("X_umap".to_string(), Embedding {
                        name: "X_umap".to_string(),
                        data,
                        n_rows: n_cells,
                        n_cols: 2,
                    });
                    2
                }
            } else {
                embeddings.insert("X_umap".to_string(), Embedding {
                    name: "X_umap".to_string(),
                    data,
                    n_rows: n_cells,
                    n_cols: 2,
                });
                2
            };
        }
    }
    
    // 检查 tSNE 坐标 (_tSNE1, _tSNE2)
    if coords.contains_key("_tSNE1") && coords.contains_key("_tSNE2") {
        let x = coords.get("_tSNE1").unwrap();
        let y = coords.get("_tSNE2").unwrap();
        
        if x.len() == n_cells && y.len() == n_cells {
            let mut data = Vec::with_capacity(n_cells * 2);
            for i in 0..n_cells {
                data.push(x[i]);
                data.push(y[i]);
            }
            
            embeddings.insert("X_tsne".to_string(), Embedding {
                name: "X_tsne".to_string(),
                data,
                n_rows: n_cells,
                n_cols: 2,
            });
        }
    }
    
    // 检查 PCA 坐标 (_PC1, _PC2, ...)
    let mut pca_coords: Vec<(usize, &Vec<f64>)> = Vec::new();
    for (name, values) in &coords {
        if name.starts_with("_PC") {
            if let Ok(idx) = name[3..].parse::<usize>() {
                pca_coords.push((idx, values));
            }
        }
    }
    
    if !pca_coords.is_empty() {
        pca_coords.sort_by_key(|(idx, _)| *idx);
        let n_pcs = pca_coords.len();
        
        // 验证所有 PC 有相同长度
        if pca_coords.iter().all(|(_, v)| v.len() == n_cells) {
            let mut data = Vec::with_capacity(n_cells * n_pcs);
            for i in 0..n_cells {
                for (_, values) in &pca_coords {
                    data.push(values[i]);
                }
            }
            
            embeddings.insert("X_pca".to_string(), Embedding {
                name: "X_pca".to_string(),
                data,
                n_rows: n_cells,
                n_cols: n_pcs,
            });
        }
    }
    
    if embeddings.is_empty() {
        Ok(None)
    } else {
        Ok(Some(embeddings))
    }
}

/// 读取属性数据集
fn read_attribute_dataset(ds: &hdf5::Dataset, expected_len: usize) -> Result<Option<ArrayRef>, CrossCellError> {
    let shape = ds.shape();
    
    // 只处理 1D 数组
    if shape.len() != 1 {
        return Ok(None);
    }
    
    let len = shape[0];
    if len != expected_len {
        return Err(CrossCellError::InvalidDimensions(format!(
            "Attribute has {} elements, expected {}",
            len, expected_len
        )));
    }
    
    // 尝试不同的数据类型
    // 1. Float64
    if let Ok(values) = ds.read_1d::<f64>() {
        let array = Float64Array::from(values.to_vec());
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 2. Int64
    if let Ok(values) = ds.read_1d::<i64>() {
        let array = Int64Array::from(values.to_vec());
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 3. Int32
    if let Ok(values) = ds.read_1d::<i32>() {
        let values_i64: Vec<i64> = values.iter().map(|&x| x as i64).collect();
        let array = Int64Array::from(values_i64);
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 4. Boolean
    if let Ok(values) = ds.read_1d::<bool>() {
        let array = BooleanArray::from(values.to_vec());
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 5. String (variable length)
    if let Ok(values) = ds.read_1d::<hdf5::types::VarLenUnicode>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 6. String (fixed length)
    if let Ok(values) = ds.read_1d::<hdf5::types::FixedUnicode<256>>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string().trim_end_matches('\0').to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Some(Arc::new(array) as ArrayRef));
    }
    
    // 无法识别的类型，跳过
    Ok(None)
}

/// 读取基因属性 (/row_attrs)
fn read_row_attrs(file: &H5File, n_genes: usize) -> Result<DataFrame, CrossCellError> {
    if !file.link_exists("row_attrs") {
        return Ok(DataFrame::empty(n_genes));
    }
    
    let row_attrs = file.group("row_attrs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open /row_attrs: {}", e))
    })?;
    
    let member_names = row_attrs.member_names().map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to list /row_attrs members: {}", e))
    })?;
    
    let mut columns: Vec<String> = Vec::new();
    let mut data: Vec<ArrayRef> = Vec::new();
    
    for name in member_names {
        if let Ok(ds) = row_attrs.dataset(&name) {
            if let Some(array) = read_attribute_dataset(&ds, n_genes)? {
                columns.push(name);
                data.push(array);
            }
        }
    }
    
    if columns.is_empty() {
        Ok(DataFrame::empty(n_genes))
    } else {
        DataFrame::new(columns, data, n_genes).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create gene metadata DataFrame: {}", e))
        })
    }
}


/// 读取 layers (/layers)
fn read_layers(file: &H5File, n_cells: usize, n_genes: usize) -> Result<Option<HashMap<String, ExpressionMatrix>>, CrossCellError> {
    if !file.link_exists("layers") {
        return Ok(None);
    }
    
    let layers_group = file.group("layers").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open /layers: {}", e))
    })?;
    
    let member_names = layers_group.member_names().map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to list /layers members: {}", e))
    })?;
    
    if member_names.is_empty() {
        return Ok(None);
    }
    
    let mut layers: HashMap<String, ExpressionMatrix> = HashMap::new();
    
    for name in member_names {
        if let Ok(ds) = layers_group.dataset(&name) {
            let shape = ds.shape();
            if shape.len() != 2 {
                continue;
            }
            
            let layer_n_genes = shape[0];
            let layer_n_cells = shape[1];
            
            if layer_n_genes != n_genes || layer_n_cells != n_cells {
                // 维度不匹配，跳过
                continue;
            }
            
            // 读取并转置 (genes × cells → cells × genes)
            if let Ok(matrix_data) = ds.read_raw::<f64>() {
                let mut transposed = vec![0.0; n_cells * n_genes];
                for i in 0..n_genes {
                    for j in 0..n_cells {
                        transposed[j * n_genes + i] = matrix_data[i * n_cells + j];
                    }
                }
                
                layers.insert(name, ExpressionMatrix::Dense(DenseMatrix {
                    n_rows: n_cells,
                    n_cols: n_genes,
                    data: transposed,
                }));
            }
        }
    }
    
    if layers.is_empty() {
        Ok(None)
    } else {
        Ok(Some(layers))
    }
}

/// 读取图数据 (/col_graphs, /row_graphs)
fn read_graphs(
    file: &H5File,
    n_cells: usize,
    n_genes: usize,
) -> Result<(Option<HashMap<String, PairwiseMatrix>>, Option<HashMap<String, PairwiseMatrix>>), CrossCellError> {
    let cell_pairwise = read_graph_group(file, "col_graphs", n_cells)?;
    let gene_pairwise = read_graph_group(file, "row_graphs", n_genes)?;
    
    Ok((cell_pairwise, gene_pairwise))
}

/// 读取图组
fn read_graph_group(
    file: &H5File,
    group_name: &str,
    n: usize,
) -> Result<Option<HashMap<String, PairwiseMatrix>>, CrossCellError> {
    if !file.link_exists(group_name) {
        return Ok(None);
    }
    
    let graphs_group = file.group(group_name).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open /{}: {}", group_name, e))
    })?;
    
    let member_names = graphs_group.member_names().map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to list /{} members: {}", group_name, e))
    })?;
    
    if member_names.is_empty() {
        return Ok(None);
    }
    
    let mut pairwise: HashMap<String, PairwiseMatrix> = HashMap::new();
    
    for name in member_names {
        // 每个图是一个 Group，包含 a, b, w 数据集
        if let Ok(graph_group) = graphs_group.group(&name) {
            if let Some(mut matrix) = read_graph_data(&graph_group, n)? {
                matrix.name = name.clone();
                pairwise.insert(name, matrix);
            }
        }
    }
    
    if pairwise.is_empty() {
        Ok(None)
    } else {
        Ok(Some(pairwise))
    }
}

/// 读取单个图数据
/// 
/// Loom 图结构：
/// - a: 边的起点索引
/// - b: 边的终点索引
/// - w: 边的权重
fn read_graph_data(group: &Group, n: usize) -> Result<Option<PairwiseMatrix>, CrossCellError> {
    // 检查必需的数据集
    if !group.link_exists("a") || !group.link_exists("b") {
        return Ok(None);
    }
    
    // 读取边的起点
    let a_ds = group.dataset("a").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open graph/a: {}", e))
    })?;
    let a: Vec<usize> = if let Ok(values) = a_ds.read_1d::<i64>() {
        values.iter().map(|&x| x as usize).collect()
    } else if let Ok(values) = a_ds.read_1d::<i32>() {
        values.iter().map(|&x| x as usize).collect()
    } else {
        return Ok(None);
    };
    
    // 读取边的终点
    let b_ds = group.dataset("b").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to open graph/b: {}", e))
    })?;
    let b: Vec<usize> = if let Ok(values) = b_ds.read_1d::<i64>() {
        values.iter().map(|&x| x as usize).collect()
    } else if let Ok(values) = b_ds.read_1d::<i32>() {
        values.iter().map(|&x| x as usize).collect()
    } else {
        return Ok(None);
    };
    
    // 读取边的权重（可选）
    let w: Vec<f64> = if group.link_exists("w") {
        if let Ok(w_ds) = group.dataset("w") {
            if let Ok(values) = w_ds.read_1d::<f64>() {
                values.to_vec()
            } else {
                vec![1.0; a.len()]
            }
        } else {
            vec![1.0; a.len()]
        }
    } else {
        vec![1.0; a.len()]
    };
    
    if a.len() != b.len() || a.len() != w.len() {
        return Err(CrossCellError::InvalidDimensions(
            "Graph edge arrays have different lengths".to_string()
        ));
    }
    
    // 构建 CSR 稀疏矩阵
    let nnz = a.len();
    
    // 计算每行的非零元素数量
    let mut row_counts = vec![0usize; n + 1];
    for &row in &a {
        if row < n {
            row_counts[row + 1] += 1;
        }
    }
    
    // 计算 indptr
    for i in 1..=n {
        row_counts[i] += row_counts[i - 1];
    }
    let indptr = row_counts;
    
    // 构建 indices 和 data
    let mut indices = vec![0usize; nnz];
    let mut data = vec![0.0; nnz];
    let mut current_pos = indptr.clone();
    
    for i in 0..nnz {
        let row = a[i];
        let col = b[i];
        let val = w[i];
        
        if row < n && col < n {
            let pos = current_pos[row];
            indices[pos] = col;
            data[pos] = val;
            current_pos[row] += 1;
        }
    }
    
    let csr = SparseMatrixCSR::new(data, indices, indptr, n, n).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create graph CSR matrix: {}", e))
    })?;
    
    Ok(Some(PairwiseMatrix {
        name: String::new(),  // Will be set by caller
        matrix: ExpressionMatrix::SparseCSR(csr),
    }))
}


/// 写入 Loom 文件
fn write_loom(data: &SingleCellData, path: &Path) -> Result<(), CrossCellError> {
    let file = H5File::create(path).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create loom file: {}", e))
    })?;
    
    let (n_cells, n_genes) = data.expression.shape();
    
    // 写入主矩阵 (需要转置: cells × genes → genes × cells)
    write_matrix(&file, "matrix", &data.expression)?;
    
    // 写入细胞属性 (/col_attrs)
    write_col_attrs(&file, &data.cell_metadata, &data.embeddings, n_cells)?;
    
    // 写入基因属性 (/row_attrs)
    write_row_attrs(&file, &data.gene_metadata, n_genes)?;
    
    // 写入 layers
    write_layers_group(&file, &data.layers, n_cells, n_genes)?;
    
    // 写入图
    write_graphs(&file, &data.cell_pairwise, &data.gene_pairwise)?;
    
    Ok(())
}

/// 写入主矩阵
fn write_matrix(
    file: &H5File,
    name: &str,
    matrix: &ExpressionMatrix,
) -> Result<(), CrossCellError> {
    let (n_rows, n_cols) = matrix.shape();
    
    // 获取矩阵数据
    let data = match matrix {
        ExpressionMatrix::Dense(dense) => dense.data.clone(),
        ExpressionMatrix::SparseCSR(sparse) => {
            // 转换为稠密矩阵
            let mut dense = vec![0.0; n_rows * n_cols];
            for row in 0..n_rows {
                let start = sparse.indptr[row];
                let end = sparse.indptr[row + 1];
                for idx in start..end {
                    let col = sparse.indices[idx];
                    dense[row * n_cols + col] = sparse.data[idx];
                }
            }
            dense
        }
        ExpressionMatrix::SparseCSC(sparse) => {
            // 转换为稠密矩阵
            let mut dense = vec![0.0; n_rows * n_cols];
            for col in 0..n_cols {
                let start = sparse.indptr[col];
                let end = sparse.indptr[col + 1];
                for idx in start..end {
                    let row = sparse.indices[idx];
                    dense[row * n_cols + col] = sparse.data[idx];
                }
            }
            dense
        }
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                return write_matrix(file, name, &cached);
            }
            return Err(CrossCellError::UnsupportedType { 
                type_name: "LazyMatrix without cached data".to_string() 
            });
        }
    };
    
    // 转置并写入 (cells × genes → genes × cells)
    let mut transposed = vec![0.0; n_rows * n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            transposed[j * n_rows + i] = data[i * n_cols + j];
        }
    }
    
    let ds = file
        .new_dataset::<f64>()
        .shape([n_cols, n_rows])  // genes × cells
        .create(name)
        .map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create dataset {}: {}", name, e))
        })?;
    
    ds.write_raw(&transposed).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to write dataset {}: {}", name, e))
    })?;
    
    Ok(())
}

/// 写入细胞属性 (/col_attrs)
fn write_col_attrs(
    file: &H5File,
    cell_metadata: &DataFrame,
    embeddings: &Option<HashMap<String, Embedding>>,
    _n_cells: usize,
) -> Result<(), CrossCellError> {
    let col_attrs = file.create_group("col_attrs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create /col_attrs: {}", e))
    })?;
    
    // 写入元数据列
    for (col_name, array) in cell_metadata.columns.iter().zip(cell_metadata.data.iter()) {
        write_attribute_to_group(&col_attrs, col_name, array)?;
    }
    
    // 写入嵌入坐标
    if let Some(embs) = embeddings {
        for (name, emb) in embs {
            write_embedding_as_attrs(&col_attrs, name, emb)?;
        }
    }
    
    Ok(())
}

/// 写入嵌入为属性
fn write_embedding_as_attrs(
    group: &hdf5::Group,
    name: &str,
    emb: &Embedding,
) -> Result<(), CrossCellError> {
    let n_cells = emb.n_rows;
    let n_dims = emb.n_cols;
    
    // 根据嵌入名称确定属性名前缀
    let prefix = if name.contains("umap") || name == "X_umap" {
        ""  // 使用 _X, _Y, _Z
    } else if name.contains("tsne") || name == "X_tsne" {
        "_tSNE"
    } else if name.contains("pca") || name == "X_pca" {
        "_PC"
    } else {
        // 其他嵌入使用原名
        name
    };
    
    if prefix.is_empty() {
        // UMAP 风格：_X, _Y, _Z
        let coord_names = ["_X", "_Y", "_Z"];
        for dim in 0..n_dims.min(3) {
            let mut values = Vec::with_capacity(n_cells);
            for i in 0..n_cells {
                values.push(emb.data[i * n_dims + dim]);
            }
            
            let ds = group
                .new_dataset::<f64>()
                .shape([n_cells])
                .create(coord_names[dim])
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", coord_names[dim], e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", coord_names[dim], e))
            })?;
        }
    } else if prefix == "_tSNE" {
        // tSNE 风格：_tSNE1, _tSNE2
        for dim in 0..n_dims.min(2) {
            let attr_name = format!("{}{}", prefix, dim + 1);
            let mut values = Vec::with_capacity(n_cells);
            for i in 0..n_cells {
                values.push(emb.data[i * n_dims + dim]);
            }
            
            let ds = group
                .new_dataset::<f64>()
                .shape([n_cells])
                .create(attr_name.as_str())
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", attr_name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", attr_name, e))
            })?;
        }
    } else if prefix == "_PC" {
        // PCA 风格：_PC1, _PC2, ...
        for dim in 0..n_dims {
            let attr_name = format!("{}{}", prefix, dim + 1);
            let mut values = Vec::with_capacity(n_cells);
            for i in 0..n_cells {
                values.push(emb.data[i * n_dims + dim]);
            }
            
            let ds = group
                .new_dataset::<f64>()
                .shape([n_cells])
                .create(attr_name.as_str())
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", attr_name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", attr_name, e))
            })?;
        }
    }
    
    Ok(())
}


/// 写入属性到组
fn write_attribute_to_group(
    group: &hdf5::Group,
    name: &str,
    array: &ArrayRef,
) -> Result<(), CrossCellError> {
    use arrow::array::{AsArray, Array};
    use arrow::datatypes::DataType;
    
    match array.data_type() {
        DataType::Float64 => {
            let float_array = array.as_primitive::<arrow::datatypes::Float64Type>();
            let values: Vec<f64> = float_array.values().to_vec();
            
            let ds = group
                .new_dataset::<f64>()
                .shape([values.len()])
                .create(name)
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
            })?;
        }
        DataType::Int64 => {
            let int_array = array.as_primitive::<arrow::datatypes::Int64Type>();
            let values: Vec<i64> = int_array.values().to_vec();
            
            let ds = group
                .new_dataset::<i64>()
                .shape([values.len()])
                .create(name)
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
            })?;
        }
        DataType::Int32 => {
            let int_array = array.as_primitive::<arrow::datatypes::Int32Type>();
            let values: Vec<i32> = int_array.values().to_vec();
            
            let ds = group
                .new_dataset::<i32>()
                .shape([values.len()])
                .create(name)
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
            })?;
        }
        DataType::Boolean => {
            let bool_array = array.as_boolean();
            let values: Vec<bool> = (0..bool_array.len())
                .map(|i| bool_array.value(i))
                .collect();
            
            let ds = group
                .new_dataset::<bool>()
                .shape([values.len()])
                .create(name)
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
            })?;
        }
        DataType::Utf8 => {
            let string_array = array.as_string::<i32>();
            let values: Vec<hdf5::types::VarLenUnicode> = (0..string_array.len())
                .map(|i| {
                    let s = string_array.value(i);
                    s.parse::<hdf5::types::VarLenUnicode>()
                        .unwrap_or_else(|_| "".parse().unwrap())
                })
                .collect();
            
            let ds = group
                .new_dataset::<hdf5::types::VarLenUnicode>()
                .shape([values.len()])
                .create(name)
                .map_err(|e| {
                    CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                })?;
            ds.write(&values).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
            })?;
        }
        DataType::Dictionary(_, _) => {
            // Categorical 列：展开为字符串
            use arrow::array::DictionaryArray;
            use arrow::datatypes::Int32Type;
            
            if let Some(dict_array) = array.as_any().downcast_ref::<DictionaryArray<Int32Type>>() {
                let keys = dict_array.keys();
                let values_arr = dict_array.values();
                
                if let Some(string_values) = values_arr.as_any().downcast_ref::<StringArray>() {
                    let strings: Vec<hdf5::types::VarLenUnicode> = (0..keys.len())
                        .map(|i| {
                            let key = keys.value(i) as usize;
                            let s = string_values.value(key);
                            s.parse::<hdf5::types::VarLenUnicode>()
                                .unwrap_or_else(|_| "".parse().unwrap())
                        })
                        .collect();
                    
                    let ds = group
                        .new_dataset::<hdf5::types::VarLenUnicode>()
                        .shape([strings.len()])
                        .create(name)
                        .map_err(|e| {
                            CrossCellError::InvalidFormat(format!("Failed to create {}: {}", name, e))
                        })?;
                    ds.write(&strings).map_err(|e| {
                        CrossCellError::InvalidFormat(format!("Failed to write {}: {}", name, e))
                    })?;
                }
            }
        }
        _ => {
            // 不支持的类型，跳过
        }
    }
    
    Ok(())
}

/// 写入基因属性 (/row_attrs)
fn write_row_attrs(
    file: &H5File,
    gene_metadata: &DataFrame,
    _n_genes: usize,
) -> Result<(), CrossCellError> {
    let row_attrs = file.create_group("row_attrs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create /row_attrs: {}", e))
    })?;
    
    // 写入元数据列
    for (col_name, array) in gene_metadata.columns.iter().zip(gene_metadata.data.iter()) {
        write_attribute_to_group(&row_attrs, col_name, array)?;
    }
    
    Ok(())
}

/// 写入 layers (/layers)
fn write_layers_group(
    file: &H5File,
    layers: &Option<HashMap<String, ExpressionMatrix>>,
    _n_cells: usize,
    _n_genes: usize,
) -> Result<(), CrossCellError> {
    let layers_group = file.create_group("layers").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create /layers: {}", e))
    })?;
    
    if let Some(layers_map) = layers {
        for (name, matrix) in layers_map {
            write_matrix_to_group(&layers_group, name, matrix)?;
        }
    }
    
    Ok(())
}

/// 写入矩阵到组
fn write_matrix_to_group(
    group: &hdf5::Group,
    name: &str,
    matrix: &ExpressionMatrix,
) -> Result<(), CrossCellError> {
    let (n_rows, n_cols) = matrix.shape();
    
    let data = match matrix {
        ExpressionMatrix::Dense(dense) => dense.data.clone(),
        ExpressionMatrix::SparseCSR(sparse) => {
            let mut dense = vec![0.0; n_rows * n_cols];
            for row in 0..n_rows {
                let start = sparse.indptr[row];
                let end = sparse.indptr[row + 1];
                for idx in start..end {
                    let col = sparse.indices[idx];
                    dense[row * n_cols + col] = sparse.data[idx];
                }
            }
            dense
        }
        ExpressionMatrix::SparseCSC(sparse) => {
            let mut dense = vec![0.0; n_rows * n_cols];
            for col in 0..n_cols {
                let start = sparse.indptr[col];
                let end = sparse.indptr[col + 1];
                for idx in start..end {
                    let row = sparse.indices[idx];
                    dense[row * n_cols + col] = sparse.data[idx];
                }
            }
            dense
        }
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                return write_matrix_to_group(group, name, &cached);
            }
            vec![0.0; n_rows * n_cols]
        }
    };
    
    // 转置 (cells × genes → genes × cells)
    let mut transposed = vec![0.0; n_rows * n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            transposed[j * n_rows + i] = data[i * n_cols + j];
        }
    }
    
    let ds = group
        .new_dataset::<f64>()
        .shape([n_cols, n_rows])  // genes × cells
        .create(name)
        .map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create layer {}: {}", name, e))
        })?;
    
    ds.write_raw(&transposed).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to write layer {}: {}", name, e))
    })?;
    
    Ok(())
}

/// 写入图数据 (/col_graphs, /row_graphs)
fn write_graphs(
    file: &H5File,
    cell_pairwise: &Option<HashMap<String, PairwiseMatrix>>,
    gene_pairwise: &Option<HashMap<String, PairwiseMatrix>>,
) -> Result<(), CrossCellError> {
    // 创建图组
    let col_graphs = file.create_group("col_graphs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create /col_graphs: {}", e))
    })?;
    
    let row_graphs = file.create_group("row_graphs").map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create /row_graphs: {}", e))
    })?;
    
    // 写入细胞图
    if let Some(cell_pw) = cell_pairwise {
        for (name, pw_matrix) in cell_pw {
            write_graph_to_group(&col_graphs, name, &pw_matrix.matrix)?;
        }
    }
    
    // 写入基因图
    if let Some(gene_pw) = gene_pairwise {
        for (name, pw_matrix) in gene_pw {
            write_graph_to_group(&row_graphs, name, &pw_matrix.matrix)?;
        }
    }
    
    Ok(())
}

/// 写入单个图到组
fn write_graph_to_group(
    parent: &hdf5::Group,
    name: &str,
    matrix: &ExpressionMatrix,
) -> Result<(), CrossCellError> {
    let graph_group = parent.create_group(name).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to create graph group {}: {}", name, e))
    })?;
    
    // 从矩阵提取边
    let (a, b, w) = match matrix {
        ExpressionMatrix::SparseCSR(csr) => {
            let mut a = Vec::new();
            let mut b = Vec::new();
            let mut w = Vec::new();
            
            for row in 0..csr.n_rows {
                let start = csr.indptr[row];
                let end = csr.indptr[row + 1];
                for idx in start..end {
                    a.push(row as i64);
                    b.push(csr.indices[idx] as i64);
                    w.push(csr.data[idx]);
                }
            }
            
            (a, b, w)
        }
        ExpressionMatrix::SparseCSC(csc) => {
            let mut a = Vec::new();
            let mut b = Vec::new();
            let mut w = Vec::new();
            
            for col in 0..csc.n_cols {
                let start = csc.indptr[col];
                let end = csc.indptr[col + 1];
                for idx in start..end {
                    a.push(csc.indices[idx] as i64);
                    b.push(col as i64);
                    w.push(csc.data[idx]);
                }
            }
            
            (a, b, w)
        }
        ExpressionMatrix::Dense(dense) => {
            let mut a = Vec::new();
            let mut b = Vec::new();
            let mut w = Vec::new();
            
            for row in 0..dense.n_rows {
                for col in 0..dense.n_cols {
                    let val = dense.data[row * dense.n_cols + col];
                    if val != 0.0 {
                        a.push(row as i64);
                        b.push(col as i64);
                        w.push(val);
                    }
                }
            }
            
            (a, b, w)
        }
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                return write_graph_to_group(parent, name, &cached);
            }
            (Vec::new(), Vec::new(), Vec::new())
        }
    };
    
    // 写入 a (起点)
    let a_ds = graph_group
        .new_dataset::<i64>()
        .shape([a.len()])
        .create("a")
        .map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create graph/a: {}", e))
        })?;
    a_ds.write(&a).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to write graph/a: {}", e))
    })?;
    
    // 写入 b (终点)
    let b_ds = graph_group
        .new_dataset::<i64>()
        .shape([b.len()])
        .create("b")
        .map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create graph/b: {}", e))
        })?;
    b_ds.write(&b).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to write graph/b: {}", e))
    })?;
    
    // 写入 w (权重)
    let w_ds = graph_group
        .new_dataset::<f64>()
        .shape([w.len()])
        .create("w")
        .map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to create graph/w: {}", e))
        })?;
    w_ds.write(&w).map_err(|e| {
        CrossCellError::InvalidFormat(format!("Failed to write graph/w: {}", e))
    })?;
    
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use arrow::array::StringArray;
    use tempfile::tempdir;
    
    #[test]
    fn test_loom_converter_info() {
        let converter = LoomConverter;
        
        assert_eq!(converter.name(), "loom");
        assert_eq!(converter.display_name(), "Loom");
        assert!(converter.can_read());
        assert!(converter.can_write());
        assert!(converter.extensions().contains(&".loom"));
    }
    
    #[test]
    fn test_loom_detect() {
        let converter = LoomConverter;
        
        assert!(converter.detect(Path::new("data.loom")));
        assert!(converter.detect(Path::new("path/to/data.LOOM")));
        assert!(!converter.detect(Path::new("data.h5ad")));
        assert!(!converter.detect(Path::new("data.rds")));
    }
    
    #[test]
    fn test_loom_roundtrip_basic() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.loom");
        
        // 创建测试数据
        let n_cells = 10;
        let n_genes = 5;
        
        let mut data = Vec::with_capacity(n_cells * n_genes);
        for i in 0..n_cells {
            for j in 0..n_genes {
                data.push((i * n_genes + j) as f64);
            }
        }
        
        let expression = ExpressionMatrix::Dense(DenseMatrix {
            n_rows: n_cells,
            n_cols: n_genes,
            data,
        });
        
        // 创建元数据
        let cell_ids: Vec<String> = (0..n_cells).map(|i| format!("cell_{}", i)).collect();
        let cell_metadata = DataFrame::new(
            vec!["CellID".to_string()],
            vec![Arc::new(StringArray::from(cell_ids)) as ArrayRef],
            n_cells,
        ).unwrap();
        
        let gene_names: Vec<String> = (0..n_genes).map(|i| format!("gene_{}", i)).collect();
        let gene_metadata = DataFrame::new(
            vec!["Gene".to_string()],
            vec![Arc::new(StringArray::from(gene_names)) as ArrayRef],
            n_genes,
        ).unwrap();
        
        let original = SingleCellData {
            expression,
            cell_metadata,
            gene_metadata,
            embeddings: None,
            layers: None,
            cell_pairwise: None,
            gene_pairwise: None,
            unstructured: None,
            spatial: None,
            metadata: DatasetMetadata::new(n_cells, n_genes, "loom".to_string()),
        };
        
        // 写入
        write_loom(&original, &path).unwrap();
        
        // 读取
        let loaded = read_loom(&path).unwrap();
        
        // 验证
        assert_eq!(loaded.metadata.n_cells, n_cells);
        assert_eq!(loaded.metadata.n_genes, n_genes);
        
        let (loaded_cells, loaded_genes) = loaded.expression.shape();
        assert_eq!(loaded_cells, n_cells);
        assert_eq!(loaded_genes, n_genes);
    }
    
    #[test]
    fn test_loom_roundtrip_with_layers() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test_layers.loom");
        
        let n_cells = 5;
        let n_genes = 3;
        
        // 主矩阵
        let main_data: Vec<f64> = (0..(n_cells * n_genes)).map(|i| i as f64).collect();
        let expression = ExpressionMatrix::Dense(DenseMatrix {
            n_rows: n_cells,
            n_cols: n_genes,
            data: main_data,
        });
        
        // Layers
        let spliced_data: Vec<f64> = (0..(n_cells * n_genes)).map(|i| (i * 2) as f64).collect();
        let unspliced_data: Vec<f64> = (0..(n_cells * n_genes)).map(|i| (i * 3) as f64).collect();
        
        let mut layers = HashMap::new();
        layers.insert("spliced".to_string(), ExpressionMatrix::Dense(DenseMatrix {
            n_rows: n_cells,
            n_cols: n_genes,
            data: spliced_data,
        }));
        layers.insert("unspliced".to_string(), ExpressionMatrix::Dense(DenseMatrix {
            n_rows: n_cells,
            n_cols: n_genes,
            data: unspliced_data,
        }));
        
        let original = SingleCellData {
            expression,
            cell_metadata: DataFrame::empty(n_cells),
            gene_metadata: DataFrame::empty(n_genes),
            embeddings: None,
            layers: Some(layers),
            cell_pairwise: None,
            gene_pairwise: None,
            unstructured: None,
            spatial: None,
            metadata: DatasetMetadata::new(n_cells, n_genes, "loom".to_string()),
        };
        
        // 写入
        write_loom(&original, &path).unwrap();
        
        // 读取
        let loaded = read_loom(&path).unwrap();
        
        // 验证 layers
        assert!(loaded.layers.is_some());
        let loaded_layers = loaded.layers.unwrap();
        assert!(loaded_layers.contains_key("spliced"));
        assert!(loaded_layers.contains_key("unspliced"));
    }
    
    #[test]
    fn test_loom_roundtrip_with_embeddings() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test_emb.loom");
        
        let n_cells = 10;
        let n_genes = 5;
        
        let main_data: Vec<f64> = (0..(n_cells * n_genes)).map(|i| i as f64).collect();
        let expression = ExpressionMatrix::Dense(DenseMatrix {
            n_rows: n_cells,
            n_cols: n_genes,
            data: main_data,
        });
        
        // 创建 UMAP 嵌入
        let mut umap_data = Vec::with_capacity(n_cells * 2);
        for i in 0..n_cells {
            umap_data.push(i as f64);  // X
            umap_data.push((i * 2) as f64);  // Y
        }
        
        let mut embeddings = HashMap::new();
        embeddings.insert("X_umap".to_string(), Embedding {
            name: "X_umap".to_string(),
            data: umap_data,
            n_rows: n_cells,
            n_cols: 2,
        });
        
        let original = SingleCellData {
            expression,
            cell_metadata: DataFrame::empty(n_cells),
            gene_metadata: DataFrame::empty(n_genes),
            embeddings: Some(embeddings),
            layers: None,
            cell_pairwise: None,
            gene_pairwise: None,
            unstructured: None,
            spatial: None,
            metadata: DatasetMetadata::new(n_cells, n_genes, "loom".to_string()),
        };
        
        // 写入
        write_loom(&original, &path).unwrap();
        
        // 读取
        let loaded = read_loom(&path).unwrap();
        
        // 验证嵌入
        assert!(loaded.embeddings.is_some());
        let loaded_embs = loaded.embeddings.unwrap();
        assert!(loaded_embs.contains_key("X_umap"));
        
        let umap = loaded_embs.get("X_umap").unwrap();
        assert_eq!(umap.n_rows, n_cells);
        assert_eq!(umap.n_cols, 2);
    }
    
    #[test]
    fn test_is_embedding_coord() {
        assert!(is_embedding_coord("_X"));
        assert!(is_embedding_coord("_Y"));
        assert!(is_embedding_coord("_Z"));
        assert!(is_embedding_coord("_tSNE1"));
        assert!(is_embedding_coord("_tSNE2"));
        assert!(is_embedding_coord("_PC1"));
        assert!(is_embedding_coord("_UMAP1"));
        
        assert!(!is_embedding_coord("CellID"));
        assert!(!is_embedding_coord("ClusterID"));
        assert!(!is_embedding_coord("Gene"));
    }
}
