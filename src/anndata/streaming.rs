//! 流式处理支持超大文件
//!
//! 实现流式读写，支持处理任意大小的文件，内存占用恒定。
//! 适用于超大数据集（10M+ 细胞）。
//!
//! # 设计原则
//! - 按块（chunk）读取/写入数据
//! - 内存占用 = chunk_size × 基因数 × 8 字节
//! - 支持 Iterator 模式
//!
//! # 使用示例
//! ```no_run
//! use crosscell::anndata::streaming::{StreamingH5adReader, DataChunk};
//! use std::path::Path;
//!
//! // 流式读取
//! let reader = StreamingH5adReader::new(Path::new("large.h5ad"), 10000).unwrap();
//! for chunk in reader {
//!     let chunk = chunk.unwrap();
//!     println!("Processing cells {}..{}", chunk.cell_indices.start, chunk.cell_indices.end);
//! }
//! ```

use super::{AnnDataError, Result};
use crate::ir::{
    DataFrame, DenseMatrix, Embedding, ExpressionMatrix, SparseMatrixCSR,
    SingleCellData, DatasetMetadata,
};
use arrow::array::{ArrayRef, Float64Array, Int32Array, Int64Array, StringArray, BooleanArray, DictionaryArray, Array};
use arrow::datatypes::Int32Type;
use hdf5::File;
use ndarray::s;
use std::collections::HashMap;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// 默认块大小（10000 个细胞）
pub const DEFAULT_CHUNK_SIZE: usize = 10000;

/// 数据块
///
/// 包含一个块的表达矩阵、细胞索引和细胞元数据。
#[derive(Debug, Clone)]
pub struct DataChunk {
    /// 表达矩阵（cells × genes）
    pub expression: ExpressionMatrix,
    /// 细胞索引范围
    pub cell_indices: Range<usize>,
    /// 细胞元数据（仅包含当前块的行）
    pub cell_metadata: DataFrame,
    /// 降维嵌入（仅包含当前块的行）
    pub embeddings: Option<HashMap<String, Embedding>>,
}

impl DataChunk {
    /// 获取块中的细胞数
    pub fn n_cells(&self) -> usize {
        self.cell_indices.end - self.cell_indices.start
    }

    /// 获取基因数
    pub fn n_genes(&self) -> usize {
        self.expression.shape().1
    }
}

/// 数据集元数据（用于流式处理）
#[derive(Debug, Clone)]
pub struct StreamingMetadata {
    /// 总细胞数
    pub n_cells: usize,
    /// 总基因数
    pub n_genes: usize,
    /// 是否为稀疏矩阵
    pub is_sparse: bool,
    /// 稀疏格式（csr_matrix 或 csc_matrix）
    pub sparse_format: Option<String>,
    /// 非零元素数量
    pub nnz: Option<usize>,
    /// 细胞元数据列名
    pub obs_columns: Vec<String>,
    /// 基因元数据列名
    pub var_columns: Vec<String>,
    /// 嵌入名称
    pub embedding_names: Vec<String>,
    /// 层名称
    pub layer_names: Vec<String>,
}

/// 流式 HDF5 读取器
///
/// 支持按块读取大型 .h5ad 文件，内存占用恒定。
pub struct StreamingH5adReader {
    /// HDF5 文件句柄
    file: File,
    /// 文件路径（用于错误消息）
    #[allow(dead_code)]
    path: PathBuf,
    /// 每块的细胞数
    chunk_size: usize,
    /// 当前块索引
    current_chunk: usize,
    /// 总块数
    total_chunks: usize,
    /// 数据集元数据
    metadata: StreamingMetadata,
    /// 是否读取元数据
    read_metadata: bool,
    /// 是否读取嵌入
    read_embeddings: bool,
    /// 基因元数据（一次性读取，不分块）
    gene_metadata: Option<DataFrame>,
    /// 完整嵌入数据（一次性读取，用于 h5ad_to_rds）
    full_embeddings: Option<HashMap<String, Embedding>>,
    /// 空间数据（一次性读取）
    spatial: Option<crate::ir::SpatialData>,
    /// 完整细胞元数据（旧格式 compound Dataset 时一次性读取）
    full_cell_metadata: Option<DataFrame>,
}

impl StreamingH5adReader {
    /// 创建新的流式读取器
    ///
    /// # 参数
    /// - `path`: .h5ad 文件路径
    /// - `chunk_size`: 每块的细胞数
    ///
    /// # 返回
    /// - `Ok(StreamingH5adReader)`: 成功创建读取器
    /// - `Err(AnnDataError)`: 创建失败
    pub fn new<P: AsRef<Path>>(path: P, chunk_size: usize) -> Result<Self> {
        let path_ref = path.as_ref();
        let file = File::open(path_ref)?;

        // 读取元数据
        let metadata = Self::read_metadata_internal(&file)?;

        // 计算总块数
        let total_chunks = (metadata.n_cells + chunk_size - 1) / chunk_size;

        // 一次性读取基因元数据（通常 < 50k 行，内存很小）
        let n_genes = metadata.n_genes;
        let gene_metadata = if file.link_exists("var") {
            Some(crate::anndata::reader::read_metadata_dataframe(&file, "var", n_genes)?)
        } else {
            None
        };

        // 一次性读取完整嵌入（n_cells × ~50 维，相对于表达矩阵很小）
        let n_cells = metadata.n_cells;
        let full_embeddings = if file.link_exists("obsm") {
            match crate::anndata::reader::read_embeddings(&file, n_cells) {
                Ok(emb) if !emb.is_empty() => Some(emb),
                _ => None,
            }
        } else {
            None
        };

        // 一次性读取空间数据
        let spatial = crate::anndata::reader::read_spatial_data(&file, n_cells)?;

        // 检测旧格式：obs 是 compound Dataset 而不是 Group
        // 旧格式无法分块读取，需要一次性读取完整细胞元数据
        let full_cell_metadata = if file.link_exists("obs") && file.group("obs").is_err() {
            // obs 存在但不是 Group → 旧格式 compound Dataset
            log::info!("Detected legacy compound Dataset for /obs, reading full cell metadata");
            Some(crate::anndata::reader::read_metadata_dataframe(&file, "obs", n_cells)?)
        } else {
            None
        };

        Ok(Self {
            file,
            path: path_ref.to_path_buf(),
            chunk_size,
            current_chunk: 0,
            total_chunks,
            metadata,
            read_metadata: true,
            read_embeddings: true,
            gene_metadata,
            full_embeddings,
            spatial,
            full_cell_metadata,
        })
    }

    /// 设置是否读取细胞元数据
    pub fn with_metadata(mut self, read_metadata: bool) -> Self {
        self.read_metadata = read_metadata;
        self
    }

    /// 设置是否读取嵌入
    pub fn with_embeddings(mut self, read_embeddings: bool) -> Self {
        self.read_embeddings = read_embeddings;
        self
    }

    /// 获取基因元数据引用
    pub fn gene_metadata(&self) -> Option<&DataFrame> {
        self.gene_metadata.as_ref()
    }

    /// 取出完整嵌入数据（消耗所有权）
    pub fn take_full_embeddings(&mut self) -> Option<HashMap<String, Embedding>> {
        self.full_embeddings.take()
    }

    /// 取出空间数据（消耗所有权）
    pub fn take_spatial(&mut self) -> Option<crate::ir::SpatialData> {
        self.spatial.take()
    }

    /// 获取数据集元数据
    pub fn metadata(&self) -> &StreamingMetadata {
        &self.metadata
    }

    /// 获取总块数
    pub fn total_chunks(&self) -> usize {
        self.total_chunks
    }

    /// 获取当前块索引
    pub fn current_chunk(&self) -> usize {
        self.current_chunk
    }

    /// 获取每块的细胞数
    pub fn chunk_size(&self) -> usize {
        self.chunk_size
    }

    /// 重置到第一个块
    pub fn reset(&mut self) {
        self.current_chunk = 0;
    }

    /// 跳转到指定块
    pub fn seek(&mut self, chunk_index: usize) -> Result<()> {
        if chunk_index >= self.total_chunks {
            return Err(AnnDataError::InvalidFormat(format!(
                "Chunk index {} out of range (total: {})",
                chunk_index, self.total_chunks
            )));
        }
        self.current_chunk = chunk_index;
        Ok(())
    }

    /// 读取下一个数据块
    pub fn next_chunk(&mut self) -> Option<Result<DataChunk>> {
        if self.current_chunk >= self.total_chunks {
            return None;
        }

        let result = self.read_chunk(self.current_chunk);
        self.current_chunk += 1;
        Some(result)
    }

    /// 读取指定块
    fn read_chunk(&self, chunk_index: usize) -> Result<DataChunk> {
        let start_row = chunk_index * self.chunk_size;
        let end_row = std::cmp::min((chunk_index + 1) * self.chunk_size, self.metadata.n_cells);
        let n_rows = end_row - start_row;

        // 读取表达矩阵块
        let expression = self.read_expression_chunk(start_row, end_row)?;

        // 读取细胞元数据块
        let cell_metadata = if self.read_metadata {
            self.read_metadata_chunk(start_row, end_row)?
        } else {
            DataFrame::empty(n_rows)
        };

        // 读取嵌入块
        let embeddings = if self.read_embeddings {
            self.read_embeddings_chunk(start_row, end_row)?
        } else {
            None
        };

        Ok(DataChunk {
            expression,
            cell_indices: start_row..end_row,
            cell_metadata,
            embeddings,
        })
    }

    /// 读取表达矩阵块
    fn read_expression_chunk(&self, start_row: usize, end_row: usize) -> Result<ExpressionMatrix> {
        if self.metadata.is_sparse {
            self.read_sparse_chunk(start_row, end_row)
        } else {
            self.read_dense_chunk(start_row, end_row)
        }
    }

    /// 读取稀疏矩阵块
    fn read_sparse_chunk(&self, start_row: usize, end_row: usize) -> Result<ExpressionMatrix> {
        let x_group = self.file.group("X")?;
        let n_rows = end_row - start_row;
        let n_cols = self.metadata.n_genes;

        // Slice read indptr（只读 start_row..=end_row）
        let indptr_dataset = x_group.dataset("indptr")?;
        let indptr_slice: Vec<i64> = if let Ok(slice_i32) = indptr_dataset.read_slice_1d::<i32, _>(s![start_row..=end_row]) {
            slice_i32.to_vec().into_iter().map(|x| x as i64).collect()
        } else {
            indptr_dataset.read_slice_1d::<i64, _>(s![start_row..=end_row])?.to_vec()
        };

        // 获取当前块的 data/indices 范围
        let data_start = indptr_slice[0] as usize;
        let data_end = indptr_slice[indptr_slice.len() - 1] as usize;

        // 调整 indptr 为相对偏移
        let indptr: Vec<usize> = indptr_slice
            .iter()
            .map(|&p| (p as usize) - data_start)
            .collect();

        // Slice read data 和 indices（只读需要的范围）
        let data_dataset = x_group.dataset("data")?;
        let indices_dataset = x_group.dataset("indices")?;

        let data: Vec<f64> = if data_start < data_end {
            data_dataset.read_slice_1d::<f64, _>(s![data_start..data_end])?.to_vec()
        } else {
            Vec::new()
        };

        let indices: Vec<usize> = if data_start < data_end {
            if let Ok(idx_i32) = indices_dataset.read_slice_1d::<i32, _>(s![data_start..data_end]) {
                idx_i32.to_vec().into_iter().map(|x| x as usize).collect()
            } else {
                indices_dataset.read_slice_1d::<i64, _>(s![data_start..data_end])?
                    .to_vec()
                    .into_iter()
                    .map(|x| x as usize)
                    .collect()
            }
        } else {
            Vec::new()
        };

        let csr = SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols)
            .map_err(|e| AnnDataError::InvalidFormat(e))?;

        Ok(ExpressionMatrix::SparseCSR(csr))
    }

    /// 读取稠密矩阵块
    fn read_dense_chunk(&self, start_row: usize, end_row: usize) -> Result<ExpressionMatrix> {
        let x_dataset = self.file.dataset("X")?;
        let n_rows = end_row - start_row;
        let n_cols = self.metadata.n_genes;

        // Slice read: 只读取需要的行
        let chunk_data = x_dataset.read_slice_2d::<f64, _>(s![start_row..end_row, ..])?;
        let data = chunk_data.into_raw_vec();

        let dense = DenseMatrix::new(data, n_rows, n_cols)
            .map_err(|e| AnnDataError::InvalidFormat(e))?;

        Ok(ExpressionMatrix::Dense(dense))
    }

    /// 读取元数据块
    fn read_metadata_chunk(&self, start_row: usize, end_row: usize) -> Result<DataFrame> {
        let n_rows = end_row - start_row;

        if !self.file.link_exists("obs") {
            return Ok(DataFrame::empty(n_rows));
        }

        // 旧格式：从预加载的完整细胞元数据中切片
        if let Some(ref full_meta) = self.full_cell_metadata {
            return Self::slice_dataframe(full_meta, start_row, end_row);
        }

        // 如果没有列，返回空 DataFrame
        if self.metadata.obs_columns.is_empty() {
            return Ok(DataFrame::empty(n_rows));
        }

        // 新版本格式：从 Group 中分块读取
        let obs_group = self.file.group("obs")?;
        let mut columns = Vec::new();
        let mut data: Vec<ArrayRef> = Vec::new();

        for col_name in &self.metadata.obs_columns {
            // 尝试作为 Dataset 打开
            if let Ok(dataset) = obs_group.dataset(col_name) {
                if let Ok(array) = self.read_column_chunk(&dataset, start_row, end_row) {
                    columns.push(col_name.clone());
                    data.push(array);
                }
            } else if let Ok(subgroup) = obs_group.group(col_name) {
                // 新版本 Categorical 列
                if subgroup.link_exists("codes") && subgroup.link_exists("categories") {
                    if let Ok(array) = self.read_categorical_chunk(&subgroup, start_row, end_row) {
                        columns.push(col_name.clone());
                        data.push(array);
                    }
                }
                // Nullable 列（nullable-integer, nullable-boolean）
                else if subgroup.link_exists("values") && subgroup.link_exists("mask") {
                    if let Ok(array) = self.read_nullable_chunk(&subgroup, start_row, end_row) {
                        columns.push(col_name.clone());
                        data.push(array);
                    }
                }
            }
        }

        DataFrame::new(columns, data, n_rows).map_err(|e| AnnDataError::InvalidFormat(e))
    }

    /// 读取列块
    fn read_column_chunk(
        &self,
        dataset: &hdf5::Dataset,
        start_row: usize,
        end_row: usize,
    ) -> Result<ArrayRef> {
        // Slice read 数值类型
        if let Ok(values) = dataset.read_slice_1d::<f64, _>(s![start_row..end_row]) {
            let array = Float64Array::from(values.to_vec());
            return Ok(Arc::new(array) as ArrayRef);
        }

        if let Ok(values) = dataset.read_slice_1d::<i64, _>(s![start_row..end_row]) {
            let array = Int64Array::from(values.to_vec());
            return Ok(Arc::new(array) as ArrayRef);
        }

        if let Ok(values) = dataset.read_slice_1d::<i32, _>(s![start_row..end_row]) {
            let array = Int32Array::from(values.to_vec());
            return Ok(Arc::new(array) as ArrayRef);
        }

        if let Ok(values) = dataset.read_slice_1d::<bool, _>(s![start_row..end_row]) {
            let array = BooleanArray::from(values.to_vec());
            return Ok(Arc::new(array) as ArrayRef);
        }

        // String 类型：VarLenUnicode 的 slice read 可能不被支持，保持全量读取再切片
        if let Ok(all_values) = dataset.read_1d::<hdf5::types::VarLenUnicode>() {
            let all_vec: Vec<_> = all_values.to_vec();
            let strings: Vec<String> = all_vec[start_row..end_row]
                .iter()
                .map(|s| s.to_string())
                .collect();
            let array = StringArray::from(strings);
            return Ok(Arc::new(array) as ArrayRef);
        }

        Err(AnnDataError::UnsupportedType(
            "Unsupported column data type".to_string(),
        ))
    }

    /// 读取 Categorical 列块
    fn read_categorical_chunk(
        &self,
        group: &hdf5::Group,
        start_row: usize,
        end_row: usize,
    ) -> Result<ArrayRef> {
        // Slice read codes（只读需要的范围）
        let codes_dataset = group.dataset("codes")?;
        let codes: Vec<i32> = if let Ok(c) = codes_dataset.read_slice_1d::<i8, _>(s![start_row..end_row]) {
            c.to_vec().into_iter().map(|x| x as i32).collect()
        } else if let Ok(c) = codes_dataset.read_slice_1d::<i16, _>(s![start_row..end_row]) {
            c.to_vec().into_iter().map(|x| x as i32).collect()
        } else if let Ok(c) = codes_dataset.read_slice_1d::<i32, _>(s![start_row..end_row]) {
            c.to_vec()
        } else {
            codes_dataset.read_slice_1d::<i64, _>(s![start_row..end_row])?
                .to_vec()
                .into_iter()
                .map(|x| x as i32)
                .collect()
        };

        // categories 保持全量读取（通常 < 1000 个类别，内存很小）
        let categories_dataset = group.dataset("categories")?;
        let categories: Vec<String> = if let Ok(cat_varlen) = categories_dataset.read_1d::<hdf5::types::VarLenUnicode>() {
            cat_varlen.iter().map(|s| s.to_string()).collect()
        } else if let Ok(cat_ascii) = categories_dataset.read_1d::<hdf5::types::VarLenAscii>() {
            cat_ascii.iter().map(|s| s.to_string()).collect()
        } else if let Ok(cat_fixed) = categories_dataset.read_1d::<hdf5::types::FixedUnicode<64>>() {
            cat_fixed.iter().map(|s| s.to_string()).collect()
        } else if let Ok(cat_fixed_ascii) = categories_dataset.read_1d::<hdf5::types::FixedAscii<64>>() {
            cat_fixed_ascii.iter().map(|s| s.to_string()).collect()
        } else if let Ok(cat_i64) = categories_dataset.read_1d::<i64>() {
            cat_i64.iter().map(|v| v.to_string()).collect()
        } else if let Ok(cat_f64) = categories_dataset.read_1d::<f64>() {
            cat_f64.iter().map(|v| v.to_string()).collect()
        } else if let Ok(cat_i32) = categories_dataset.read_1d::<i32>() {
            cat_i32.iter().map(|v| v.to_string()).collect()
        } else {
            return Err(AnnDataError::UnsupportedType(
                "Categorical labels must be string type".to_string(),
            ));
        };

        // 创建 Arrow Dictionary Array
        let keys = Int32Array::from(codes);
        let values = StringArray::from(categories);
        let dict_array = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values))
            .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to create dictionary array: {}", e)))?;

        Ok(Arc::new(dict_array) as ArrayRef)
    }

    /// 读取 Nullable 列块（nullable-integer, nullable-boolean）
    fn read_nullable_chunk(
        &self,
        subgroup: &hdf5::Group,
        start_row: usize,
        end_row: usize,
    ) -> Result<ArrayRef> {
        let values_dataset = subgroup.dataset("values")?;
        let mask_dataset = subgroup.dataset("mask")?;
        
        // Slice read mask（true = NA）
        let mask_slice: Vec<bool> = mask_dataset.read_slice_1d::<bool, _>(s![start_row..end_row])
            .map(|v| v.to_vec())
            .unwrap_or_else(|_| vec![false; end_row - start_row]);
        
        // validity: Arrow true = valid, anndata mask true = NA
        let validity: Vec<bool> = mask_slice.iter().map(|&m| !m).collect();
        
        // Slice read 数值类型
        if let Ok(vals) = values_dataset.read_slice_1d::<i64, _>(s![start_row..end_row]) {
            let vals = vals.to_vec();
            let array = Int64Array::from(
                vals.iter().zip(validity.iter())
                    .map(|(&v, &valid)| if valid { Some(v) } else { None })
                    .collect::<Vec<Option<i64>>>()
            );
            return Ok(Arc::new(array) as ArrayRef);
        }
        
        if let Ok(vals) = values_dataset.read_slice_1d::<i32, _>(s![start_row..end_row]) {
            let vals = vals.to_vec();
            let array = Int32Array::from(
                vals.iter().zip(validity.iter())
                    .map(|(&v, &valid)| if valid { Some(v) } else { None })
                    .collect::<Vec<Option<i32>>>()
            );
            return Ok(Arc::new(array) as ArrayRef);
        }
        
        if let Ok(vals) = values_dataset.read_slice_1d::<f64, _>(s![start_row..end_row]) {
            let vals = vals.to_vec();
            let array = Float64Array::from(
                vals.iter().zip(validity.iter())
                    .map(|(&v, &valid)| if valid { Some(v) } else { None })
                    .collect::<Vec<Option<f64>>>()
            );
            return Ok(Arc::new(array) as ArrayRef);
        }
        
        if let Ok(vals) = values_dataset.read_slice_1d::<bool, _>(s![start_row..end_row]) {
            let vals = vals.to_vec();
            let array = BooleanArray::from(
                vals.iter().zip(validity.iter())
                    .map(|(&v, &valid)| if valid { Some(v) } else { None })
                    .collect::<Vec<Option<bool>>>()
            );
            return Ok(Arc::new(array) as ArrayRef);
        }
        
        Err(AnnDataError::UnsupportedType(
            "Unsupported nullable column data type".to_string(),
        ))
    }

    /// 读取嵌入块
    fn read_embeddings_chunk(
        &self,
        start_row: usize,
        end_row: usize,
    ) -> Result<Option<HashMap<String, Embedding>>> {
        if !self.file.link_exists("obsm") {
            return Ok(None);
        }

        let obsm_group = self.file.group("obsm")?;
        let n_rows = end_row - start_row;
        let mut embeddings = HashMap::new();

        for name in &self.metadata.embedding_names {
            if let Ok(dataset) = obsm_group.dataset(name) {
                let shape = dataset.shape();
                if shape.len() != 2 {
                    continue;
                }
                let n_cols = shape[1] as usize;

                // Slice read: 只读取需要的行
                let chunk_data = if let Ok(arr) = dataset.read_slice_2d::<f64, _>(s![start_row..end_row, ..]) {
                    arr.into_raw_vec()
                } else if let Ok(arr) = dataset.read_slice_2d::<f32, _>(s![start_row..end_row, ..]) {
                    arr.into_raw_vec().into_iter().map(|x| x as f64).collect()
                } else {
                    continue;
                };

                let embedding = Embedding {
                    name: name.clone(),
                    data: chunk_data,
                    n_rows,
                    n_cols,
                };
                embeddings.insert(name.clone(), embedding);
            }
        }

        if embeddings.is_empty() {
            Ok(None)
        } else {
            Ok(Some(embeddings))
        }
    }

    /// 从 DataFrame 中切片指定行范围
    fn slice_dataframe(df: &DataFrame, start_row: usize, end_row: usize) -> Result<DataFrame> {
        let n_rows = end_row - start_row;
        let sliced_data: Vec<ArrayRef> = df.data.iter().map(|arr| {
            arr.slice(start_row, n_rows)
        }).collect();
        DataFrame::new(df.columns.clone(), sliced_data, n_rows)
            .map_err(|e| AnnDataError::InvalidFormat(e))
    }

    /// 内部方法：读取元数据
    fn read_metadata_internal(file: &File) -> Result<StreamingMetadata> {
        // 读取表达矩阵元数据
        let (n_cells, n_genes, is_sparse, sparse_format, nnz) = Self::read_expression_metadata(file)?;

        // 读取 obs 列名
        let obs_columns = if file.link_exists("obs") {
            Self::list_columns(file, "obs")?
        } else {
            Vec::new()
        };

        // 读取 var 列名
        let var_columns = if file.link_exists("var") {
            Self::list_columns(file, "var")?
        } else {
            Vec::new()
        };

        // 读取嵌入名称
        let embedding_names = if file.link_exists("obsm") {
            Self::list_group_members(file, "obsm")?
        } else {
            Vec::new()
        };

        // 读取层名称
        let layer_names = if file.link_exists("layers") {
            Self::list_group_members(file, "layers")?
        } else {
            Vec::new()
        };

        Ok(StreamingMetadata {
            n_cells,
            n_genes,
            is_sparse,
            sparse_format,
            nnz,
            obs_columns,
            var_columns,
            embedding_names,
            layer_names,
        })
    }

    /// 读取表达矩阵元数据
    fn read_expression_metadata(file: &File) -> Result<(usize, usize, bool, Option<String>, Option<usize>)> {
        if !file.link_exists("X") {
            return Err(AnnDataError::MissingField("/X".to_string()));
        }

        // 尝试作为 Dataset 打开（稠密矩阵）
        match file.dataset("X") {
            Ok(x_dataset) => {
                let shape = x_dataset.shape();
                if shape.len() != 2 {
                    return Err(AnnDataError::InvalidFormat(format!(
                        "Expected 2D matrix, got {}D",
                        shape.len()
                    )));
                }
                let (n_rows, n_cols) = (shape[0] as usize, shape[1] as usize);
                return Ok((n_rows, n_cols, false, None, None));
            }
            Err(e) => {
                log::debug!("file.dataset(\"X\") failed: {:?}", e);
            }
        }

        // 尝试作为 Group 打开（稀疏矩阵）
        match file.group("X") {
            Ok(x_group) => {
                // 尝试读取 shape 属性，支持 int64 和 usize 两种类型
                let shape: Vec<usize> = if let Ok(shape_attr) = x_group.attr("shape") {
                    // 先尝试 usize (u64)
                    if let Ok(arr) = shape_attr.read_1d::<usize>() {
                        arr.to_vec()
                    } else if let Ok(arr) = shape_attr.read_1d::<i64>() {
                        // numpy 默认存储为 int64，需要转换
                        log::debug!("/X shape attribute is int64, converting to usize");
                        arr.iter().map(|&v| v as usize).collect()
                    } else {
                        return Err(AnnDataError::InvalidFormat(
                            "Cannot read /X shape attribute as usize or i64".to_string()
                        ));
                    }
                } else if let Ok(shape_attr) = x_group.attr("h5sparse_shape") {
                    if let Ok(arr) = shape_attr.read_1d::<usize>() {
                        arr.to_vec()
                    } else if let Ok(arr) = shape_attr.read_1d::<i64>() {
                        arr.iter().map(|&v| v as usize).collect()
                    } else {
                        return Err(AnnDataError::InvalidFormat(
                            "Cannot read /X h5sparse_shape attribute".to_string()
                        ));
                    }
                } else {
                    return Err(AnnDataError::MissingField("/X shape attribute".to_string()));
                };

                if shape.len() != 2 {
                    return Err(AnnDataError::InvalidFormat(format!(
                        "Expected shape to have 2 dimensions, got {}",
                        shape.len()
                    )));
                }

                let (n_rows, n_cols) = (shape[0], shape[1]);

                let sparse_format = if let Ok(encoding_attr) = x_group.attr("encoding-type") {
                    if let Ok(encoding) = encoding_attr.read_scalar::<hdf5::types::VarLenUnicode>() {
                        Some(encoding.to_string())
                    } else {
                        Some("csr_matrix".to_string())
                    }
                } else {
                    Some("csr_matrix".to_string())
                };

                let nnz = if let Ok(data_dataset) = x_group.dataset("data") {
                    Some(data_dataset.shape()[0] as usize)
                } else {
                    None
                };

                return Ok((n_rows, n_cols, true, sparse_format, nnz));
            }
            Err(e) => {
                log::debug!("file.group(\"X\") failed: {:?}", e);
            }
        }

        // 如果两种方式都失败，输出诊断信息
        log::error!("Both file.dataset(\"X\") and file.group(\"X\") failed");
        log::error!("link_exists(\"X\") = {}", file.link_exists("X"));
        if let Ok(members) = file.member_names() {
            log::error!("Root group members: {:?}", members);
        }

        Err(AnnDataError::InvalidFormat(
            "/X is neither a Group (sparse) nor a Dataset (dense)".to_string(),
        ))
    }

    /// 列出 Group 中的列名
    fn list_columns(file: &File, group_name: &str) -> Result<Vec<String>> {
        if let Ok(group) = file.group(group_name) {
            let members = group.member_names()?;
            Ok(members
                .into_iter()
                .filter(|n| !n.starts_with("__"))
                .collect())
        } else {
            Ok(Vec::new())
        }
    }

    /// 列出 Group 中的成员名称
    fn list_group_members(file: &File, group_name: &str) -> Result<Vec<String>> {
        if let Ok(group) = file.group(group_name) {
            let members = group.member_names()?;
            Ok(members
                .into_iter()
                .filter(|n| !n.starts_with("__"))
                .collect())
        } else {
            Ok(Vec::new())
        }
    }
}

impl Iterator for StreamingH5adReader {
    type Item = Result<DataChunk>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_chunk()
    }
}


/// 流式 HDF5 写入器
///
/// 支持按块写入大型 .h5ad 文件。
/// 
/// 注意：由于 HDF5 crate 的限制，当前实现会在内存中累积数据，
/// 然后在 `finish()` 时一次性写入。对于真正的流式写入，
/// 需要使用支持可扩展数据集的 HDF5 库。
pub struct StreamingH5adWriter {
    /// 输出文件路径
    path: PathBuf,
    /// 数据集元数据
    metadata: StreamingMetadata,
    /// 累积的表达矩阵数据
    accumulated_data: Vec<f64>,
    /// 累积的列索引
    accumulated_indices: Vec<usize>,
    /// 累积的行指针
    accumulated_indptr: Vec<usize>,
    /// 已写入的行数
    written_rows: usize,
    /// 基因元数据
    gene_metadata: Option<DataFrame>,
    /// 累积的细胞元数据
    accumulated_cell_metadata: Vec<(String, Vec<ArrayRef>)>,
    /// 累积的嵌入
    accumulated_embeddings: HashMap<String, Vec<f64>>,
    /// 嵌入维度
    embedding_dims: HashMap<String, usize>,
}

impl StreamingH5adWriter {
    /// 创建新的流式写入器
    ///
    /// # 参数
    /// - `path`: 输出 .h5ad 文件路径
    /// - `metadata`: 数据集元数据
    /// - `gene_metadata`: 基因元数据（可选）
    ///
    /// # 返回
    /// - `Ok(StreamingH5adWriter)`: 成功创建写入器
    /// - `Err(AnnDataError)`: 创建失败
    pub fn new<P: AsRef<Path>>(
        path: P,
        metadata: StreamingMetadata,
        gene_metadata: Option<DataFrame>,
    ) -> Result<Self> {
        let path_ref = path.as_ref();

        Ok(Self {
            path: path_ref.to_path_buf(),
            metadata,
            accumulated_data: Vec::new(),
            accumulated_indices: Vec::new(),
            accumulated_indptr: vec![0],
            written_rows: 0,
            gene_metadata,
            accumulated_cell_metadata: Vec::new(),
            accumulated_embeddings: HashMap::new(),
            embedding_dims: HashMap::new(),
        })
    }

    /// 获取已写入的行数
    pub fn written_rows(&self) -> usize {
        self.written_rows
    }

    /// 写入一个数据块
    ///
    /// # 参数
    /// - `chunk`: 数据块
    ///
    /// # 返回
    /// - `Ok(())`: 成功写入
    /// - `Err(AnnDataError)`: 写入失败
    pub fn write_chunk(&mut self, chunk: &DataChunk) -> Result<()> {
        // 验证块
        let expected_start = self.written_rows;
        if chunk.cell_indices.start != expected_start {
            return Err(AnnDataError::InvalidFormat(format!(
                "Expected chunk starting at row {}, got {}",
                expected_start, chunk.cell_indices.start
            )));
        }

        // 累积表达矩阵数据
        self.accumulate_expression(&chunk.expression)?;

        // 累积细胞元数据
        self.accumulate_metadata(&chunk.cell_metadata)?;

        // 累积嵌入
        if let Some(ref embeddings) = chunk.embeddings {
            self.accumulate_embeddings(embeddings)?;
        }

        self.written_rows = chunk.cell_indices.end;
        Ok(())
    }

    /// 累积表达矩阵数据
    fn accumulate_expression(&mut self, expression: &ExpressionMatrix) -> Result<()> {
        match expression {
            ExpressionMatrix::SparseCSR(csr) => {
                // 累积 CSR 数据
                let offset = self.accumulated_data.len();
                self.accumulated_data.extend_from_slice(&csr.data);
                self.accumulated_indices.extend_from_slice(&csr.indices);
                
                // 调整 indptr 并追加（跳过第一个 0）
                for &ptr in csr.indptr.iter().skip(1) {
                    self.accumulated_indptr.push(ptr + offset);
                }
            }
            ExpressionMatrix::SparseCSC(csc) => {
                // 转换为 CSR 后累积
                let csr = crate::sparse::convert::csc_to_csr(csc);
                let offset = self.accumulated_data.len();
                self.accumulated_data.extend_from_slice(&csr.data);
                self.accumulated_indices.extend_from_slice(&csr.indices);
                
                for &ptr in csr.indptr.iter().skip(1) {
                    self.accumulated_indptr.push(ptr + offset);
                }
            }
            ExpressionMatrix::Dense(dense) => {
                // 转换为 CSR 后累积
                let csr = crate::sparse::convert::dense_to_csr(dense);
                let offset = self.accumulated_data.len();
                self.accumulated_data.extend_from_slice(&csr.data);
                self.accumulated_indices.extend_from_slice(&csr.indices);
                
                for &ptr in csr.indptr.iter().skip(1) {
                    self.accumulated_indptr.push(ptr + offset);
                }
            }
            ExpressionMatrix::Lazy(_) => {
                return Err(AnnDataError::InvalidFormat(
                    "Cannot write LazyMatrix chunk".to_string(),
                ));
            }
        }
        Ok(())
    }

    /// 累积元数据
    fn accumulate_metadata(&mut self, cell_metadata: &DataFrame) -> Result<()> {
        // 如果是第一个块，初始化列结构
        if self.accumulated_cell_metadata.is_empty() {
            for col_name in &cell_metadata.columns {
                self.accumulated_cell_metadata.push((col_name.clone(), Vec::new()));
            }
        }

        // 累积每列的数据
        for (i, array) in cell_metadata.data.iter().enumerate() {
            if i < self.accumulated_cell_metadata.len() {
                self.accumulated_cell_metadata[i].1.push(array.clone());
            }
        }

        Ok(())
    }

    /// 累积嵌入
    fn accumulate_embeddings(&mut self, embeddings: &HashMap<String, Embedding>) -> Result<()> {
        for (name, embedding) in embeddings {
            let entry = self.accumulated_embeddings.entry(name.clone()).or_insert_with(Vec::new);
            entry.extend_from_slice(&embedding.data);
            self.embedding_dims.insert(name.clone(), embedding.n_cols);
        }
        Ok(())
    }

    /// 完成写入
    ///
    /// 将累积的数据写入 HDF5 文件。
    pub fn finish(self) -> Result<()> {
        // 验证写入的行数
        if self.written_rows != self.metadata.n_cells {
            return Err(AnnDataError::InvalidFormat(format!(
                "Expected {} rows, wrote {}",
                self.metadata.n_cells, self.written_rows
            )));
        }

        // 创建 HDF5 文件
        let file = File::create(&self.path)?;

        // 写入表达矩阵
        self.write_expression_matrix(&file)?;

        // 写入基因元数据
        if let Some(ref gene_metadata) = self.gene_metadata {
            self.write_gene_metadata(&file, gene_metadata)?;
        } else {
            file.create_group("var")?;
        }

        // 写入细胞元数据
        self.write_cell_metadata(&file)?;

        // 写入嵌入
        self.write_embeddings(&file)?;

        Ok(())
    }

    /// 写入表达矩阵
    fn write_expression_matrix(&self, file: &File) -> Result<()> {
        let x_group = file.create_group("X")?;

        // 写入 shape 属性
        let shape = vec![self.metadata.n_cells as i64, self.metadata.n_genes as i64];
        x_group
            .new_attr::<i64>()
            .shape([2])
            .create("shape")?
            .write(&shape)?;

        // 写入 encoding-type 属性
        x_group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("encoding-type")?
            .write_scalar(
                &"csr_matrix"
                    .parse::<hdf5::types::VarLenUnicode>()
                    .map_err(|_| {
                        hdf5::Error::Internal("Failed to create VarLenUnicode".to_string())
                    })?,
            )?;

        // 写入 encoding-version 属性
        x_group
            .new_attr::<hdf5::types::VarLenUnicode>()
            .create("encoding-version")?
            .write_scalar(
                &"0.1.0"
                    .parse::<hdf5::types::VarLenUnicode>()
                    .map_err(|_| {
                        hdf5::Error::Internal("Failed to create VarLenUnicode".to_string())
                    })?,
            )?;

        // 写入 data
        let data_dataset = x_group
            .new_dataset::<f64>()
            .shape([self.accumulated_data.len()])
            .create("data")?;
        data_dataset.write(&self.accumulated_data)?;

        // 写入 indices
        let indices_i32: Vec<i32> = self.accumulated_indices.iter().map(|&x| x as i32).collect();
        let indices_dataset = x_group
            .new_dataset::<i32>()
            .shape([indices_i32.len()])
            .create("indices")?;
        indices_dataset.write(&indices_i32)?;

        // 写入 indptr
        let indptr_i32: Vec<i32> = self.accumulated_indptr.iter().map(|&x| x as i32).collect();
        let indptr_dataset = x_group
            .new_dataset::<i32>()
            .shape([indptr_i32.len()])
            .create("indptr")?;
        indptr_dataset.write(&indptr_i32)?;

        Ok(())
    }

    /// 写入基因元数据
    fn write_gene_metadata(&self, file: &File, gene_metadata: &DataFrame) -> Result<()> {
        use arrow::array::{Array, AsArray};
        use arrow::datatypes::DataType;

        let var_group = file.create_group("var")?;

        for (col_name, array) in gene_metadata.columns.iter().zip(gene_metadata.data.iter()) {
            match array.data_type() {
                DataType::Float64 => {
                    let float_array = array.as_primitive::<arrow::datatypes::Float64Type>();
                    let values: Vec<f64> = float_array.values().to_vec();
                    let dataset = var_group
                        .new_dataset::<f64>()
                        .shape([values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&values)?;
                }
                DataType::Int64 => {
                    let int_array = array.as_primitive::<arrow::datatypes::Int64Type>();
                    let values: Vec<i64> = int_array.values().to_vec();
                    let dataset = var_group
                        .new_dataset::<i64>()
                        .shape([values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&values)?;
                }
                DataType::Utf8 => {
                    let string_array = array.as_string::<i32>();
                    let values: Vec<hdf5::types::VarLenUnicode> = (0..string_array.len())
                        .map(|i| {
                            let s = string_array.value(i);
                            s.parse::<hdf5::types::VarLenUnicode>().map_err(|_| {
                                hdf5::Error::Internal("Failed to create VarLenUnicode".to_string())
                            })
                        })
                        .collect::<std::result::Result<Vec<_>, _>>()?;
                    let dataset = var_group
                        .new_dataset::<hdf5::types::VarLenUnicode>()
                        .shape([values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&values)?;
                }
                _ => {}
            }
        }

        Ok(())
    }

    /// 写入细胞元数据
    fn write_cell_metadata(&self, file: &File) -> Result<()> {
        let obs_group = file.create_group("obs")?;

        // 合并并写入每列
        for (col_name, arrays) in &self.accumulated_cell_metadata {
            if arrays.is_empty() {
                continue;
            }

            // 根据第一个数组的类型决定如何合并
            let first = &arrays[0];
            match first.data_type() {
                arrow::datatypes::DataType::Float64 => {
                    let mut all_values = Vec::new();
                    for arr in arrays {
                        let float_arr = arr.as_any().downcast_ref::<Float64Array>().unwrap();
                        all_values.extend(float_arr.values().iter().cloned());
                    }
                    let dataset = obs_group
                        .new_dataset::<f64>()
                        .shape([all_values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&all_values)?;
                }
                arrow::datatypes::DataType::Int64 => {
                    let mut all_values = Vec::new();
                    for arr in arrays {
                        let int_arr = arr.as_any().downcast_ref::<Int64Array>().unwrap();
                        all_values.extend(int_arr.values().iter().cloned());
                    }
                    let dataset = obs_group
                        .new_dataset::<i64>()
                        .shape([all_values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&all_values)?;
                }
                arrow::datatypes::DataType::Int32 => {
                    let mut all_values = Vec::new();
                    for arr in arrays {
                        let int_arr = arr.as_any().downcast_ref::<Int32Array>().unwrap();
                        all_values.extend(int_arr.values().iter().cloned());
                    }
                    let dataset = obs_group
                        .new_dataset::<i32>()
                        .shape([all_values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&all_values)?;
                }
                arrow::datatypes::DataType::Boolean => {
                    let mut all_values = Vec::new();
                    for arr in arrays {
                        let bool_arr = arr.as_any().downcast_ref::<BooleanArray>().unwrap();
                        for i in 0..bool_arr.len() {
                            all_values.push(bool_arr.value(i));
                        }
                    }
                    let dataset = obs_group
                        .new_dataset::<bool>()
                        .shape([all_values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&all_values)?;
                }
                arrow::datatypes::DataType::Utf8 => {
                    let mut all_values: Vec<hdf5::types::VarLenUnicode> = Vec::new();
                    for arr in arrays {
                        let str_arr = arr.as_any().downcast_ref::<StringArray>().unwrap();
                        for i in 0..str_arr.len() {
                            let s = str_arr.value(i);
                            let unicode = s.parse::<hdf5::types::VarLenUnicode>()
                                .map_err(|_| hdf5::Error::Internal("Failed to create VarLenUnicode".to_string()))?;
                            all_values.push(unicode);
                        }
                    }
                    let dataset = obs_group
                        .new_dataset::<hdf5::types::VarLenUnicode>()
                        .shape([all_values.len()])
                        .create(col_name.as_str())?;
                    dataset.write(&all_values)?;
                }
                _ => {}
            }
        }

        Ok(())
    }

    /// 写入嵌入
    fn write_embeddings(&self, file: &File) -> Result<()> {
        if self.accumulated_embeddings.is_empty() {
            return Ok(());
        }

        let obsm_group = file.create_group("obsm")?;

        for (name, data) in &self.accumulated_embeddings {
            let n_cols = self.embedding_dims.get(name).copied().unwrap_or(2);
            let n_rows = data.len() / n_cols;

            let dataset = obsm_group
                .new_dataset::<f64>()
                .shape([n_rows, n_cols])
                .create(name.as_str())?;
            dataset.write_raw(data)?;
        }

        Ok(())
    }
}

/// 流式格式转换
///
/// 支持 .h5ad → .h5ad 的流式转换（可用于重新分块或压缩）。
///
/// # 参数
/// - `input`: 输入文件路径
/// - `output`: 输出文件路径
/// - `chunk_size`: 每块的细胞数
/// - `progress_callback`: 进度回调函数（可选）
///
/// # 返回
/// - `Ok(())`: 成功转换
/// - `Err(AnnDataError)`: 转换失败
pub fn streaming_convert<P: AsRef<Path>, F>(
    input: P,
    output: P,
    chunk_size: usize,
    mut progress_callback: Option<F>,
) -> Result<()>
where
    F: FnMut(usize, usize), // (processed_cells, total_cells)
{
    // 创建读取器
    let reader = StreamingH5adReader::new(input, chunk_size)?;
    let metadata = reader.metadata().clone();
    let total_cells = metadata.n_cells;

    // 创建写入器
    let mut writer = StreamingH5adWriter::new(output, metadata, None)?;

    // 流式转换
    for chunk_result in reader {
        let chunk = chunk_result?;
        let processed = chunk.cell_indices.end;

        writer.write_chunk(&chunk)?;

        // 调用进度回调
        if let Some(ref mut callback) = progress_callback {
            callback(processed, total_cells);
        }
    }

    writer.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_metadata() {
        let metadata = StreamingMetadata {
            n_cells: 10000,
            n_genes: 5000,
            is_sparse: true,
            sparse_format: Some("csr_matrix".to_string()),
            nnz: Some(500000),
            obs_columns: vec!["cell_type".to_string(), "batch".to_string()],
            var_columns: vec!["gene_name".to_string()],
            embedding_names: vec!["X_pca".to_string(), "X_umap".to_string()],
            layer_names: vec!["counts".to_string()],
        };

        assert_eq!(metadata.n_cells, 10000);
        assert_eq!(metadata.n_genes, 5000);
        assert!(metadata.is_sparse);
    }

    #[test]
    fn test_data_chunk() {
        let data = vec![1.0, 2.0, 3.0];
        let indices = vec![0, 1, 2];
        let indptr = vec![0, 1, 2, 3];
        let csr = SparseMatrixCSR::new(data, indices, indptr, 3, 10).unwrap();

        let chunk = DataChunk {
            expression: ExpressionMatrix::SparseCSR(csr),
            cell_indices: 0..3,
            cell_metadata: DataFrame::empty(3),
            embeddings: None,
        };

        assert_eq!(chunk.n_cells(), 3);
        assert_eq!(chunk.n_genes(), 10);
    }
}


// ============================================================================
// Streaming RDS Support
// ============================================================================

/// 流式 RDS 读取器
///
/// RDS 格式不原生支持流式读取，因此这个实现采用以下策略：
/// 1. 先读取完整的 RDS 文件结构
/// 2. 将表达矩阵按块返回
///
/// 注意：由于 RDS 格式的限制，整个文件仍然需要先加载到内存中。
/// 这个实现主要用于与 StreamingH5adReader 保持一致的 API。
pub struct StreamingRdsReader {
    /// 已加载的 IR 数据
    data: crate::ir::SingleCellData,
    /// 每块的细胞数
    chunk_size: usize,
    /// 当前块索引
    current_chunk: usize,
    /// 总块数
    total_chunks: usize,
}

impl StreamingRdsReader {
    /// 创建新的流式 RDS 读取器
    ///
    /// # 参数
    /// - `path`: .rds 文件路径
    /// - `chunk_size`: 每块的细胞数
    ///
    /// # 返回
    /// - `Ok(StreamingRdsReader)`: 成功创建读取器
    /// - `Err(AnnDataError)`: 创建失败
    pub fn new<P: AsRef<Path>>(path: P, chunk_size: usize) -> Result<Self> {
        // 读取完整的 RDS 文件
        let path_str = path.as_ref().to_str()
            .ok_or_else(|| AnnDataError::InvalidFormat("Invalid path".to_string()))?;
        let data = crate::seurat::seurat_rds_to_ir(path_str)
            .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to read RDS: {}", e)))?;

        let n_cells = data.metadata.n_cells;
        let total_chunks = (n_cells + chunk_size - 1) / chunk_size;

        Ok(Self {
            data,
            chunk_size,
            current_chunk: 0,
            total_chunks,
        })
    }

    /// 获取数据集元数据
    pub fn metadata(&self) -> StreamingMetadata {
        StreamingMetadata {
            n_cells: self.data.metadata.n_cells,
            n_genes: self.data.metadata.n_genes,
            is_sparse: self.data.expression.is_sparse(),
            sparse_format: if self.data.expression.is_sparse() {
                Some("csc_matrix".to_string()) // R 默认使用 CSC
            } else {
                None
            },
            nnz: Some(self.data.expression.nnz()),
            obs_columns: self.data.cell_metadata.columns.clone(),
            var_columns: self.data.gene_metadata.columns.clone(),
            embedding_names: self.data.embeddings.as_ref()
                .map(|e| e.keys().cloned().collect())
                .unwrap_or_default(),
            layer_names: self.data.layers.as_ref()
                .map(|l| l.keys().cloned().collect())
                .unwrap_or_default(),
        }
    }

    /// 获取总块数
    pub fn total_chunks(&self) -> usize {
        self.total_chunks
    }

    /// 获取当前块索引
    pub fn current_chunk(&self) -> usize {
        self.current_chunk
    }

    /// 重置到第一个块
    pub fn reset(&mut self) {
        self.current_chunk = 0;
    }

    /// 读取下一个数据块
    pub fn next_chunk(&mut self) -> Option<Result<DataChunk>> {
        if self.current_chunk >= self.total_chunks {
            return None;
        }

        let result = self.read_chunk(self.current_chunk);
        self.current_chunk += 1;
        Some(result)
    }

    /// 读取指定块
    fn read_chunk(&self, chunk_index: usize) -> Result<DataChunk> {
        let start_row = chunk_index * self.chunk_size;
        let end_row = std::cmp::min((chunk_index + 1) * self.chunk_size, self.data.metadata.n_cells);
        let _n_rows = end_row - start_row;
        let n_cols = self.data.metadata.n_genes;

        // 提取表达矩阵块
        let expression = self.extract_expression_chunk(start_row, end_row, n_cols)?;

        // 提取细胞元数据块
        let cell_metadata = self.extract_metadata_chunk(start_row, end_row)?;

        // 提取嵌入块
        let embeddings = self.extract_embeddings_chunk(start_row, end_row)?;

        Ok(DataChunk {
            expression,
            cell_indices: start_row..end_row,
            cell_metadata,
            embeddings,
        })
    }

    /// 提取表达矩阵块
    fn extract_expression_chunk(&self, start_row: usize, end_row: usize, n_cols: usize) -> Result<ExpressionMatrix> {
        let n_rows = end_row - start_row;

        match &self.data.expression {
            ExpressionMatrix::SparseCSR(csr) => {
                let data_start = csr.indptr[start_row];
                let data_end = csr.indptr[end_row];

                let data: Vec<f64> = csr.data[data_start..data_end].to_vec();
                let indices: Vec<usize> = csr.indices[data_start..data_end].to_vec();
                let indptr: Vec<usize> = csr.indptr[start_row..=end_row]
                    .iter()
                    .map(|&p| p - data_start)
                    .collect();

                let chunk_csr = SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols)
                    .map_err(|e| AnnDataError::InvalidFormat(e))?;
                Ok(ExpressionMatrix::SparseCSR(chunk_csr))
            }
            ExpressionMatrix::SparseCSC(csc) => {
                // CSC 格式按列存储，提取行需要遍历所有列
                let mut new_data = Vec::new();
                let mut new_indices = Vec::new();
                let mut new_indptr = vec![0usize];

                for col in 0..n_cols {
                    let col_start = csc.indptr[col];
                    let col_end = csc.indptr[col + 1];

                    for idx in col_start..col_end {
                        let row = csc.indices[idx];
                        if row >= start_row && row < end_row {
                            new_data.push(csc.data[idx]);
                            new_indices.push(row - start_row);
                        }
                    }
                    new_indptr.push(new_data.len());
                }

                let chunk_csc = crate::ir::SparseMatrixCSC::new(new_data, new_indices, new_indptr, n_rows, n_cols)
                    .map_err(|e| AnnDataError::InvalidFormat(e))?;
                
                // 转换为 CSR 以保持一致性
                let chunk_csr = crate::sparse::convert::csc_to_csr(&chunk_csc);
                Ok(ExpressionMatrix::SparseCSR(chunk_csr))
            }
            ExpressionMatrix::Dense(dense) => {
                let mut data = Vec::with_capacity(n_rows * n_cols);
                for row in start_row..end_row {
                    let row_start = row * n_cols;
                    let row_end = row_start + n_cols;
                    data.extend_from_slice(&dense.data[row_start..row_end]);
                }

                let chunk_dense = DenseMatrix::new(data, n_rows, n_cols)
                    .map_err(|e| AnnDataError::InvalidFormat(e))?;
                Ok(ExpressionMatrix::Dense(chunk_dense))
            }
            ExpressionMatrix::Lazy(_) => {
                Err(AnnDataError::InvalidFormat("Cannot chunk LazyMatrix".to_string()))
            }
        }
    }

    /// 提取元数据块
    fn extract_metadata_chunk(&self, start_row: usize, end_row: usize) -> Result<DataFrame> {
        let n_rows = end_row - start_row;
        let mut columns = Vec::new();
        let mut data: Vec<ArrayRef> = Vec::new();

        for (col_name, array) in self.data.cell_metadata.columns.iter().zip(self.data.cell_metadata.data.iter()) {
            let sliced = self.slice_array(array, start_row, end_row)?;
            columns.push(col_name.clone());
            data.push(sliced);
        }

        DataFrame::new(columns, data, n_rows).map_err(|e| AnnDataError::InvalidFormat(e))
    }

    /// 切片 Arrow 数组
    fn slice_array(&self, array: &ArrayRef, start: usize, end: usize) -> Result<ArrayRef> {
        // Arrow 的 slice 方法返回一个视图，我们需要创建新的数组
        let sliced = array.slice(start, end - start);
        Ok(sliced)
    }

    /// 提取嵌入块
    fn extract_embeddings_chunk(&self, start_row: usize, end_row: usize) -> Result<Option<HashMap<String, Embedding>>> {
        let embeddings = match &self.data.embeddings {
            Some(emb) => emb,
            None => return Ok(None),
        };

        let n_rows = end_row - start_row;
        let mut result = HashMap::new();

        for (name, embedding) in embeddings {
            let n_cols = embedding.n_cols;
            let mut data = Vec::with_capacity(n_rows * n_cols);

            for row in start_row..end_row {
                let row_start = row * n_cols;
                let row_end = row_start + n_cols;
                if row_end <= embedding.data.len() {
                    data.extend_from_slice(&embedding.data[row_start..row_end]);
                }
            }

            result.insert(name.clone(), Embedding {
                name: name.clone(),
                data,
                n_rows,
                n_cols,
            });
        }

        if result.is_empty() {
            Ok(None)
        } else {
            Ok(Some(result))
        }
    }
}

impl Iterator for StreamingRdsReader {
    type Item = Result<DataChunk>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_chunk()
    }
}

/// 流式格式转换（RDS → H5AD）
///
/// 将 RDS 文件流式转换为 H5AD 文件。
///
/// # 参数
/// - `input`: 输入 RDS 文件路径
/// - `output`: 输出 H5AD 文件路径
/// - `chunk_size`: 每块的细胞数
/// - `progress_callback`: 进度回调函数（可选）
pub fn streaming_rds_to_h5ad<P: AsRef<Path>, F>(
    input: P,
    output: P,
    chunk_size: usize,
    mut progress_callback: Option<F>,
) -> Result<()>
where
    F: FnMut(usize, usize),
{
    let reader = StreamingRdsReader::new(input, chunk_size)?;
    let metadata = reader.metadata();
    let total_cells = metadata.n_cells;

    let mut writer = StreamingH5adWriter::new(output, metadata, None)?;

    for chunk_result in reader {
        let chunk = chunk_result?;
        let processed = chunk.cell_indices.end;

        writer.write_chunk(&chunk)?;

        if let Some(ref mut callback) = progress_callback {
            callback(processed, total_cells);
        }
    }

    writer.finish()?;
    Ok(())
}

/// 流式格式转换（H5AD → RDS）
///
/// 将 H5AD 文件流式转换为 RDS 文件。
///
/// 注意：由于 RDS 格式的限制，这个函数会在内存中累积所有数据，
/// 然后一次性写入 RDS 文件。
///
/// # 参数
/// - `input`: 输入 H5AD 文件路径
/// - `output`: 输出 RDS 文件路径
/// - `chunk_size`: 每块的细胞数（用于读取）
/// - `progress_callback`: 进度回调函数（可选）
pub fn streaming_h5ad_to_rds<P: AsRef<Path>, F>(
    input: P,
    output: P,
    chunk_size: usize,
    mut progress_callback: Option<F>,
) -> Result<()>
where
    F: FnMut(usize, usize),
{
    let mut reader = StreamingH5adReader::new(&input, chunk_size)?;
    let metadata = reader.metadata().clone();
    let total_cells = metadata.n_cells;
    let n_genes = metadata.n_genes;

    // 从 reader 获取已读取的非分块数据（不再二次读取！）
    let gene_metadata = reader.gene_metadata().cloned()
        .unwrap_or_else(|| DataFrame::empty(n_genes));
    let full_embeddings = reader.take_full_embeddings();
    let spatial = reader.take_spatial();

    // 累积表达矩阵和细胞元数据（分块读取）
    let mut all_data = Vec::new();
    let mut all_indices = Vec::new();
    let mut all_indptr = vec![0usize];
    let mut all_cell_metadata_columns: Option<Vec<String>> = None;
    let mut all_cell_metadata_arrays: Vec<Vec<ArrayRef>> = Vec::new();
    let mut total_accumulated_rows = 0usize;

    for chunk_result in &mut reader {
        let chunk = chunk_result?;
        let processed = chunk.cell_indices.end;

        // 累积表达矩阵
        match &chunk.expression {
            ExpressionMatrix::SparseCSR(csr) => {
                let offset = all_data.len();
                all_data.extend_from_slice(&csr.data);
                all_indices.extend_from_slice(&csr.indices);
                for &ptr in csr.indptr.iter().skip(1) {
                    all_indptr.push(ptr + offset);
                }
            }
            _ => {
                // 转换为 CSR
                let csr = match &chunk.expression {
                    ExpressionMatrix::SparseCSC(csc) => crate::sparse::convert::csc_to_csr(csc),
                    ExpressionMatrix::Dense(dense) => crate::sparse::convert::dense_to_csr(dense),
                    _ => continue,
                };
                let offset = all_data.len();
                all_data.extend_from_slice(&csr.data);
                all_indices.extend_from_slice(&csr.indices);
                for &ptr in csr.indptr.iter().skip(1) {
                    all_indptr.push(ptr + offset);
                }
            }
        }

        // 累积细胞元数据
        {
            let cell_meta = &chunk.cell_metadata;
            if all_cell_metadata_columns.is_none() {
                all_cell_metadata_columns = Some(cell_meta.columns.clone());
                all_cell_metadata_arrays = cell_meta.columns.iter().map(|_| Vec::new()).collect();
            }
            for (i, arr) in cell_meta.data.iter().enumerate() {
                if i < all_cell_metadata_arrays.len() {
                    all_cell_metadata_arrays[i].push(arr.clone());
                }
            }
            total_accumulated_rows += cell_meta.n_rows;
        }

        if let Some(ref mut callback) = progress_callback {
            callback(processed, total_cells);
        }
    }

    // 构建完整表达矩阵
    let expression = SparseMatrixCSR::new(all_data, all_indices, all_indptr, metadata.n_cells, metadata.n_genes)
        .map_err(|e| AnnDataError::InvalidFormat(e))?;

    // 合并细胞元数据
    let cell_metadata = if let Some(columns) = all_cell_metadata_columns {
        let mut merged_data: Vec<ArrayRef> = Vec::new();
        for col_arrays in &all_cell_metadata_arrays {
            if col_arrays.is_empty() {
                continue;
            }
            // 使用 arrow concat 合并数组
            let refs: Vec<&dyn arrow::array::Array> = col_arrays.iter().map(|a| a.as_ref()).collect();
            let merged = arrow::compute::concat(&refs)
                .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to concat arrays: {}", e)))?;
            merged_data.push(merged);
        }
        DataFrame::new(columns, merged_data, total_accumulated_rows)
            .map_err(|e| AnnDataError::InvalidFormat(e))?
    } else {
        DataFrame::empty(total_cells)
    };

    // 构建完整 IR（不再调用 read_h5ad）
    let dataset_metadata = DatasetMetadata::new(total_cells, n_genes, "anndata".to_string());
    let mut ir_data = SingleCellData::new(
        ExpressionMatrix::SparseCSR(expression),
        cell_metadata,
        gene_metadata,
        dataset_metadata,
    ).map_err(|e| AnnDataError::InvalidFormat(e))?;
    ir_data.embeddings = full_embeddings;
    ir_data.spatial = spatial;

    // 写入 RDS
    let output_str = output.as_ref().to_str()
        .ok_or_else(|| AnnDataError::InvalidFormat("Invalid output path".to_string()))?;
    crate::seurat::write_seurat_rds(&ir_data, output_str)
        .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to write RDS: {}", e)))?;

    Ok(())
}

#[cfg(test)]
mod rds_tests {
    use super::*;

    #[test]
    fn test_streaming_rds_metadata() {
        // 这个测试需要实际的 RDS 文件，跳过
    }
}
