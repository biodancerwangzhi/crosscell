//! 表达矩阵 IR 类型
//!
//! 定义了单细胞数据的表达矩阵结构，支持稠密和稀疏格式。
//! - DenseMatrix: 稠密矩阵（行优先存储）
//! - SparseMatrixCSR: 压缩稀疏行格式（Python scipy 默认）
//! - SparseMatrixCSC: 压缩稀疏列格式（R Matrix 默认）
//! - LazyMatrix: 延迟加载矩阵（支持超大数据集）
//! - ChunkedMatrix: 分块矩阵（支持超大数据集的分块存储和读取）

use std::fmt;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

/// 表达矩阵枚举，支持稠密和稀疏格式
#[derive(Debug, Clone, PartialEq)]
pub enum ExpressionMatrix {
    Dense(DenseMatrix),
    SparseCSR(SparseMatrixCSR),
    SparseCSC(SparseMatrixCSC),
    /// 延迟加载矩阵（用于超大数据集）
    Lazy(LazyMatrix),
}

/// 缓存策略
#[derive(Debug, Clone, PartialEq)]
pub enum CachePolicy {
    /// 不缓存，每次都从磁盘读取
    None,
    /// 完全缓存，第一次读取后保留在内存
    Full,
    /// LRU 缓存，限制内存使用（字节数）
    LRU(usize),
}

impl Default for CachePolicy {
    fn default() -> Self {
        CachePolicy::Full
    }
}

/// 存储后端类型
#[derive(Debug, Clone, PartialEq)]
pub enum BackendType {
    /// HDF5 格式（.h5ad）
    HDF5,
    /// RDS 格式（.rds）
    RDS,
    /// 内存后端（用于测试）
    Memory,
}

/// 延迟加载矩阵
/// 
/// 支持按需加载数据，显著降低内存占用。
/// 适用于超大数据集（100 万+ 细胞）。
#[derive(Debug)]
pub struct LazyMatrix {
    /// 存储后端类型
    pub backend_type: BackendType,
    /// 文件路径
    pub file_path: PathBuf,
    /// 数据集路径（在文件内的路径，如 "/X"）
    pub dataset_path: String,
    /// 矩阵维度 (n_rows, n_cols)
    pub shape: (usize, usize),
    /// 非零元素数量（如果已知）
    pub nnz: Option<usize>,
    /// 是否为稀疏矩阵
    pub is_sparse: bool,
    /// 稀疏格式（CSR 或 CSC）
    pub sparse_format: Option<SparseFormat>,
    /// 缓存策略
    pub cache_policy: CachePolicy,
    /// 缓存的数据（线程安全）
    cache: Arc<Mutex<Option<CachedMatrix>>>,
}

/// 稀疏格式类型
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SparseFormat {
    CSR,
    CSC,
}

/// 缓存的矩阵数据
#[derive(Debug, Clone)]
enum CachedMatrix {
    Dense(DenseMatrix),
    SparseCSR(SparseMatrixCSR),
    SparseCSC(SparseMatrixCSC),
}

impl Clone for LazyMatrix {
    fn clone(&self) -> Self {
        Self {
            backend_type: self.backend_type.clone(),
            file_path: self.file_path.clone(),
            dataset_path: self.dataset_path.clone(),
            shape: self.shape,
            nnz: self.nnz,
            is_sparse: self.is_sparse,
            sparse_format: self.sparse_format,
            cache_policy: self.cache_policy.clone(),
            cache: Arc::new(Mutex::new(None)), // 克隆时不复制缓存
        }
    }
}

impl PartialEq for LazyMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.backend_type == other.backend_type
            && self.file_path == other.file_path
            && self.dataset_path == other.dataset_path
            && self.shape == other.shape
            && self.is_sparse == other.is_sparse
            && self.sparse_format == other.sparse_format
    }
}

impl LazyMatrix {
    /// 创建新的延迟加载矩阵
    pub fn new(
        backend_type: BackendType,
        file_path: PathBuf,
        dataset_path: String,
        shape: (usize, usize),
        is_sparse: bool,
        sparse_format: Option<SparseFormat>,
    ) -> Self {
        Self {
            backend_type,
            file_path,
            dataset_path,
            shape,
            nnz: None,
            is_sparse,
            sparse_format,
            cache_policy: CachePolicy::default(),
            cache: Arc::new(Mutex::new(None)),
        }
    }

    /// 设置非零元素数量
    pub fn with_nnz(mut self, nnz: usize) -> Self {
        self.nnz = Some(nnz);
        self
    }

    /// 设置缓存策略
    pub fn with_cache_policy(mut self, policy: CachePolicy) -> Self {
        self.cache_policy = policy;
        self
    }

    /// 检查数据是否已缓存
    pub fn is_cached(&self) -> bool {
        self.cache.lock().unwrap().is_some()
    }

    /// 清除缓存
    pub fn clear_cache(&self) {
        let mut cache = self.cache.lock().unwrap();
        *cache = None;
    }

    /// 估算内存使用量（字节）
    pub fn estimate_memory_usage(&self) -> usize {
        let (n_rows, n_cols) = self.shape;
        if self.is_sparse {
            // 稀疏矩阵：data + indices + indptr
            // 假设 95% 稀疏度
            let estimated_nnz = self.nnz.unwrap_or((n_rows * n_cols) / 20);
            estimated_nnz * (8 + 8) + (n_rows + 1) * 8 // data(f64) + indices(usize) + indptr
        } else {
            // 稠密矩阵
            n_rows * n_cols * 8 // f64
        }
    }

    /// 获取矩阵形状
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// 缓存稠密矩阵数据
    pub fn cache_dense(&self, matrix: DenseMatrix) {
        let mut cache = self.cache.lock().unwrap();
        *cache = Some(CachedMatrix::Dense(matrix));
    }

    /// 缓存 CSR 稀疏矩阵数据
    pub fn cache_csr(&self, matrix: SparseMatrixCSR) {
        let mut cache = self.cache.lock().unwrap();
        *cache = Some(CachedMatrix::SparseCSR(matrix));
    }

    /// 缓存 CSC 稀疏矩阵数据
    pub fn cache_csc(&self, matrix: SparseMatrixCSC) {
        let mut cache = self.cache.lock().unwrap();
        *cache = Some(CachedMatrix::SparseCSC(matrix));
    }

    /// 获取缓存的矩阵（如果有）
    pub fn get_cached(&self) -> Option<ExpressionMatrix> {
        let cache = self.cache.lock().unwrap();
        cache.as_ref().map(|c| match c {
            CachedMatrix::Dense(m) => ExpressionMatrix::Dense(m.clone()),
            CachedMatrix::SparseCSR(m) => ExpressionMatrix::SparseCSR(m.clone()),
            CachedMatrix::SparseCSC(m) => ExpressionMatrix::SparseCSC(m.clone()),
        })
    }
}

/// 分块矩阵（用于超大数据集的分块存储和读取）
///
/// 借鉴 easySCF 的设计，将大矩阵分割为多个小块，每块默认 5000 个细胞。
/// 这种设计可以：
/// 1. 避免 R 的 int32 限制（最大 2^31-1 = 2,147,483,647）
/// 2. 支持超大数据集（10M+ 细胞）
/// 3. 降低内存峰值（只需加载单个块）
/// 4. 支持并行处理
#[derive(Debug, Clone)]
pub struct ChunkedMatrix {
    /// 每个块的最大行数（细胞数），默认 5000
    pub chunk_size: usize,
    /// 所有块的列表
    pub chunks: Vec<ExpressionMatrix>,
    /// 总行数（细胞数）
    pub total_rows: usize,
    /// 总列数（基因数）
    pub total_cols: usize,
    /// 是否为稀疏矩阵
    pub is_sparse: bool,
    /// 稀疏格式（如果是稀疏矩阵）
    pub sparse_format: Option<SparseFormat>,
}

impl ChunkedMatrix {
    /// 默认块大小（5000 个细胞，与 easySCF 一致）
    pub const DEFAULT_CHUNK_SIZE: usize = 5000;

    /// 创建新的分块矩阵
    pub fn new(
        chunk_size: usize,
        total_rows: usize,
        total_cols: usize,
        is_sparse: bool,
        sparse_format: Option<SparseFormat>,
    ) -> Self {
        Self {
            chunk_size,
            chunks: Vec::new(),
            total_rows,
            total_cols,
            is_sparse,
            sparse_format,
        }
    }

    /// 从完整矩阵创建分块矩阵
    ///
    /// # 参数
    /// - `matrix`: 完整的表达矩阵
    /// - `chunk_size`: 每个块的最大行数
    ///
    /// # 返回
    /// - `Ok(ChunkedMatrix)`: 分块后的矩阵
    /// - `Err(String)`: 分块失败
    pub fn from_matrix(matrix: ExpressionMatrix, chunk_size: usize) -> Result<Self, String> {
        let (total_rows, total_cols) = matrix.shape();
        let is_sparse = matrix.is_sparse();
        let sparse_format = match &matrix {
            ExpressionMatrix::SparseCSR(_) => Some(SparseFormat::CSR),
            ExpressionMatrix::SparseCSC(_) => Some(SparseFormat::CSC),
            _ => None,
        };

        let mut chunked = Self::new(chunk_size, total_rows, total_cols, is_sparse, sparse_format);

        // 计算块数
        let num_chunks = (total_rows + chunk_size - 1) / chunk_size;

        // 分割矩阵
        for i in 0..num_chunks {
            let start_row = i * chunk_size;
            let end_row = std::cmp::min((i + 1) * chunk_size, total_rows);
            let chunk = Self::extract_rows(&matrix, start_row, end_row)?;
            chunked.chunks.push(chunk);
        }

        Ok(chunked)
    }

    /// 从分块矩阵合并为完整矩阵
    ///
    /// # 返回
    /// - `Ok(ExpressionMatrix)`: 合并后的完整矩阵
    /// - `Err(String)`: 合并失败
    pub fn to_matrix(&self) -> Result<ExpressionMatrix, String> {
        if self.chunks.is_empty() {
            // 返回空矩阵
            if self.is_sparse {
                match self.sparse_format {
                    Some(SparseFormat::CSR) => {
                        return Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR::empty(
                            self.total_rows,
                            self.total_cols,
                        )));
                    }
                    Some(SparseFormat::CSC) => {
                        return Ok(ExpressionMatrix::SparseCSC(SparseMatrixCSC::empty(
                            self.total_rows,
                            self.total_cols,
                        )));
                    }
                    None => {
                        return Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR::empty(
                            self.total_rows,
                            self.total_cols,
                        )));
                    }
                }
            } else {
                return Ok(ExpressionMatrix::Dense(DenseMatrix::zeros(
                    self.total_rows,
                    self.total_cols,
                )));
            }
        }

        // 根据第一个块的类型决定合并方式
        match &self.chunks[0] {
            ExpressionMatrix::Dense(_) => self.merge_dense_chunks(),
            ExpressionMatrix::SparseCSR(_) => self.merge_csr_chunks(),
            ExpressionMatrix::SparseCSC(_) => self.merge_csc_chunks(),
            ExpressionMatrix::Lazy(_) => {
                Err("Cannot merge lazy chunks directly".to_string())
            }
        }
    }

    /// 获取块数
    pub fn num_chunks(&self) -> usize {
        self.chunks.len()
    }

    /// 获取指定块
    pub fn get_chunk(&self, index: usize) -> Option<&ExpressionMatrix> {
        self.chunks.get(index)
    }

    /// 添加一个块
    pub fn add_chunk(&mut self, chunk: ExpressionMatrix) -> Result<(), String> {
        let (chunk_rows, chunk_cols) = chunk.shape();
        
        // 验证列数一致
        if chunk_cols != self.total_cols {
            return Err(format!(
                "Chunk column count {} does not match total columns {}",
                chunk_cols, self.total_cols
            ));
        }

        // 验证行数不超过 chunk_size（最后一个块可以小于 chunk_size）
        if chunk_rows > self.chunk_size {
            return Err(format!(
                "Chunk row count {} exceeds chunk size {}",
                chunk_rows, self.chunk_size
            ));
        }

        self.chunks.push(chunk);
        Ok(())
    }

    /// 获取矩阵形状
    pub fn shape(&self) -> (usize, usize) {
        (self.total_rows, self.total_cols)
    }

    /// 估算内存使用量（字节）
    pub fn estimate_memory_usage(&self) -> usize {
        self.chunks.iter().map(|c| c.estimate_memory_usage()).sum()
    }

    /// 估算单个块的内存使用量（字节）
    pub fn estimate_chunk_memory(&self) -> usize {
        let chunk_rows = std::cmp::min(self.chunk_size, self.total_rows);
        if self.is_sparse {
            // 假设 95% 稀疏度
            let estimated_nnz = (chunk_rows * self.total_cols) / 20;
            estimated_nnz * (8 + 8) + (chunk_rows + 1) * 8
        } else {
            chunk_rows * self.total_cols * 8
        }
    }

    /// 从矩阵中提取指定行范围
    fn extract_rows(
        matrix: &ExpressionMatrix,
        start_row: usize,
        end_row: usize,
    ) -> Result<ExpressionMatrix, String> {
        match matrix {
            ExpressionMatrix::Dense(m) => {
                let chunk_rows = end_row - start_row;
                let start_idx = start_row * m.n_cols;
                let end_idx = end_row * m.n_cols;
                let data = m.data[start_idx..end_idx].to_vec();
                Ok(ExpressionMatrix::Dense(DenseMatrix::new(
                    data,
                    chunk_rows,
                    m.n_cols,
                )?))
            }
            ExpressionMatrix::SparseCSR(m) => {
                let chunk_rows = end_row - start_row;
                let start_ptr = m.indptr[start_row];
                let end_ptr = m.indptr[end_row];

                let data = m.data[start_ptr..end_ptr].to_vec();
                let indices = m.indices[start_ptr..end_ptr].to_vec();
                let indptr: Vec<usize> = m.indptr[start_row..=end_row]
                    .iter()
                    .map(|&p| p - start_ptr)
                    .collect();

                Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR::new(
                    data,
                    indices,
                    indptr,
                    chunk_rows,
                    m.n_cols,
                )?))
            }
            ExpressionMatrix::SparseCSC(m) => {
                // CSC 格式按列存储，提取行需要遍历所有列
                let chunk_rows = end_row - start_row;
                let mut new_data = Vec::new();
                let mut new_indices = Vec::new();
                let mut new_indptr = vec![0usize];

                for col in 0..m.n_cols {
                    let col_start = m.indptr[col];
                    let col_end = m.indptr[col + 1];

                    for idx in col_start..col_end {
                        let row = m.indices[idx];
                        if row >= start_row && row < end_row {
                            new_data.push(m.data[idx]);
                            new_indices.push(row - start_row);
                        }
                    }
                    new_indptr.push(new_data.len());
                }

                Ok(ExpressionMatrix::SparseCSC(SparseMatrixCSC::new(
                    new_data,
                    new_indices,
                    new_indptr,
                    chunk_rows,
                    m.n_cols,
                )?))
            }
            ExpressionMatrix::Lazy(_) => {
                Err("Cannot extract rows from lazy matrix".to_string())
            }
        }
    }

    /// 合并稠密块
    fn merge_dense_chunks(&self) -> Result<ExpressionMatrix, String> {
        let mut data = Vec::with_capacity(self.total_rows * self.total_cols);

        for chunk in &self.chunks {
            if let ExpressionMatrix::Dense(m) = chunk {
                data.extend_from_slice(&m.data);
            } else {
                return Err("Mixed chunk types".to_string());
            }
        }

        Ok(ExpressionMatrix::Dense(DenseMatrix::new(
            data,
            self.total_rows,
            self.total_cols,
        )?))
    }

    /// 合并 CSR 块
    fn merge_csr_chunks(&self) -> Result<ExpressionMatrix, String> {
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for chunk in &self.chunks {
            if let ExpressionMatrix::SparseCSR(m) = chunk {
                let offset = data.len();
                data.extend_from_slice(&m.data);
                indices.extend_from_slice(&m.indices);

                // 更新 indptr（跳过第一个 0）
                for &ptr in m.indptr.iter().skip(1) {
                    indptr.push(ptr + offset);
                }
            } else {
                return Err("Mixed chunk types".to_string());
            }
        }

        Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR::new(
            data,
            indices,
            indptr,
            self.total_rows,
            self.total_cols,
        )?))
    }

    /// 合并 CSC 块
    fn merge_csc_chunks(&self) -> Result<ExpressionMatrix, String> {
        // CSC 格式合并比较复杂，需要按列合并
        let mut col_data: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.total_cols];
        let mut current_row_offset = 0;

        for chunk in &self.chunks {
            if let ExpressionMatrix::SparseCSC(m) = chunk {
                for col in 0..m.n_cols {
                    let col_start = m.indptr[col];
                    let col_end = m.indptr[col + 1];

                    for idx in col_start..col_end {
                        let row = m.indices[idx] + current_row_offset;
                        let value = m.data[idx];
                        col_data[col].push((row, value));
                    }
                }
                current_row_offset += m.n_rows;
            } else {
                return Err("Mixed chunk types".to_string());
            }
        }

        // 构建合并后的 CSC 矩阵
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for col_entries in col_data {
            for (row, value) in col_entries {
                indices.push(row);
                data.push(value);
            }
            indptr.push(data.len());
        }

        Ok(ExpressionMatrix::SparseCSC(SparseMatrixCSC::new(
            data,
            indices,
            indptr,
            self.total_rows,
            self.total_cols,
        )?))
    }
}

impl PartialEq for ChunkedMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.chunk_size == other.chunk_size
            && self.total_rows == other.total_rows
            && self.total_cols == other.total_cols
            && self.is_sparse == other.is_sparse
            && self.sparse_format == other.sparse_format
            && self.chunks.len() == other.chunks.len()
    }
}

impl fmt::Display for ChunkedMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let format_str = if self.is_sparse {
            match self.sparse_format {
                Some(SparseFormat::CSR) => "CSR",
                Some(SparseFormat::CSC) => "CSC",
                None => "Sparse",
            }
        } else {
            "Dense"
        };
        write!(
            f,
            "ChunkedMatrix({} × {}, {} chunks of {}, {})",
            self.total_rows,
            self.total_cols,
            self.chunks.len(),
            self.chunk_size,
            format_str
        )
    }
}

impl ExpressionMatrix {
    /// 获取矩阵形状 (n_rows, n_cols)
    pub fn shape(&self) -> (usize, usize) {
        match self {
            ExpressionMatrix::Dense(m) => (m.n_rows, m.n_cols),
            ExpressionMatrix::SparseCSR(m) => (m.n_rows, m.n_cols),
            ExpressionMatrix::SparseCSC(m) => (m.n_rows, m.n_cols),
            ExpressionMatrix::Lazy(m) => m.shape,
        }
    }

    /// 验证矩阵结构的一致性
    pub fn validate(&self) -> Result<(), String> {
        match self {
            ExpressionMatrix::Dense(m) => m.validate(),
            ExpressionMatrix::SparseCSR(m) => m.validate(),
            ExpressionMatrix::SparseCSC(m) => m.validate(),
            ExpressionMatrix::Lazy(m) => {
                // Lazy 矩阵只验证元数据
                if m.shape.0 == 0 && m.shape.1 == 0 {
                    Ok(())
                } else if m.shape.0 == 0 || m.shape.1 == 0 {
                    // 允许空矩阵，但不允许只有一个维度为 0
                    Err(format!(
                        "LazyMatrix: invalid shape {:?}, both dimensions must be 0 or both non-zero",
                        m.shape
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }

    /// 获取非零元素数量
    pub fn nnz(&self) -> usize {
        match self {
            ExpressionMatrix::Dense(m) => m.data.iter().filter(|&&x| x != 0.0).count(),
            ExpressionMatrix::SparseCSR(m) => m.data.len(),
            ExpressionMatrix::SparseCSC(m) => m.data.len(),
            ExpressionMatrix::Lazy(m) => m.nnz.unwrap_or(0),
        }
    }

    /// 检查是否为延迟加载矩阵
    pub fn is_lazy(&self) -> bool {
        matches!(self, ExpressionMatrix::Lazy(_))
    }

    /// 检查是否为稀疏矩阵
    pub fn is_sparse(&self) -> bool {
        match self {
            ExpressionMatrix::Dense(_) => false,
            ExpressionMatrix::SparseCSR(_) => true,
            ExpressionMatrix::SparseCSC(_) => true,
            ExpressionMatrix::Lazy(m) => m.is_sparse,
        }
    }

    /// 估算内存使用量（字节）
    pub fn estimate_memory_usage(&self) -> usize {
        match self {
            ExpressionMatrix::Dense(m) => m.data.len() * 8,
            ExpressionMatrix::SparseCSR(m) => {
                m.data.len() * 8 + m.indices.len() * 8 + m.indptr.len() * 8
            }
            ExpressionMatrix::SparseCSC(m) => {
                m.data.len() * 8 + m.indices.len() * 8 + m.indptr.len() * 8
            }
            ExpressionMatrix::Lazy(m) => m.estimate_memory_usage(),
        }
    }
}

impl fmt::Display for ExpressionMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExpressionMatrix::Dense(m) => write!(f, "DenseMatrix({} × {})", m.n_rows, m.n_cols),
            ExpressionMatrix::SparseCSR(m) => {
                write!(
                    f,
                    "SparseCSR({} × {}, {} nnz)",
                    m.n_rows,
                    m.n_cols,
                    m.data.len()
                )
            }
            ExpressionMatrix::SparseCSC(m) => {
                write!(
                    f,
                    "SparseCSC({} × {}, {} nnz)",
                    m.n_rows,
                    m.n_cols,
                    m.data.len()
                )
            }
            ExpressionMatrix::Lazy(m) => {
                let format_str = if m.is_sparse {
                    match m.sparse_format {
                        Some(SparseFormat::CSR) => "CSR",
                        Some(SparseFormat::CSC) => "CSC",
                        None => "Sparse",
                    }
                } else {
                    "Dense"
                };
                let cached_str = if m.is_cached() { ", cached" } else { "" };
                let nnz_str = m.nnz.map(|n| format!(", {} nnz", n)).unwrap_or_default();
                write!(
                    f,
                    "LazyMatrix({} × {}, {}{}{})",
                    m.shape.0, m.shape.1, format_str, nnz_str, cached_str
                )
            }
        }
    }
}

/// 稠密矩阵（行优先存储）
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    /// 矩阵数据（行优先，cells × genes）
    pub data: Vec<f64>,
    /// 行数（细胞数）
    pub n_rows: usize,
    /// 列数（基因数）
    pub n_cols: usize,
}

impl DenseMatrix {
    /// 创建新的稠密矩阵
    pub fn new(data: Vec<f64>, n_rows: usize, n_cols: usize) -> Result<Self, String> {
        let matrix = Self {
            data,
            n_rows,
            n_cols,
        };
        matrix.validate()?;
        Ok(matrix)
    }

    /// 创建零矩阵
    pub fn zeros(n_rows: usize, n_cols: usize) -> Self {
        Self {
            data: vec![0.0; n_rows * n_cols],
            n_rows,
            n_cols,
        }
    }

    /// 验证矩阵维度一致性
    pub fn validate(&self) -> Result<(), String> {
        let expected_len = self.n_rows * self.n_cols;
        if self.data.len() != expected_len {
            return Err(format!(
                "DenseMatrix: data length {} does not match dimensions {} × {} = {}",
                self.data.len(),
                self.n_rows,
                self.n_cols,
                expected_len
            ));
        }
        Ok(())
    }

    /// 获取指定位置的元素
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.n_rows || col >= self.n_cols {
            return None;
        }
        Some(self.data[row * self.n_cols + col])
    }
}

/// CSR 格式稀疏矩阵（Compressed Sparse Row，Python scipy 默认）
///
/// 存储格式：
/// - data: 非零元素值
/// - indices: 每个非零元素的列索引
/// - indptr: 每行的起始位置（长度 = n_rows + 1）
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrixCSR {
    /// 非零元素值
    pub data: Vec<f64>,
    /// 列索引（对应每个非零元素）
    pub indices: Vec<usize>,
    /// 行指针（长度 = n_rows + 1）
    pub indptr: Vec<usize>,
    /// 行数
    pub n_rows: usize,
    /// 列数
    pub n_cols: usize,
}

impl SparseMatrixCSR {
    /// 创建新的 CSR 稀疏矩阵
    pub fn new(
        data: Vec<f64>,
        indices: Vec<usize>,
        indptr: Vec<usize>,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Self, String> {
        let matrix = Self {
            data,
            indices,
            indptr,
            n_rows,
            n_cols,
        };
        matrix.validate()?;
        Ok(matrix)
    }

    /// 创建空的 CSR 矩阵
    pub fn empty(n_rows: usize, n_cols: usize) -> Self {
        Self {
            data: Vec::new(),
            indices: Vec::new(),
            indptr: vec![0; n_rows + 1],
            n_rows,
            n_cols,
        }
    }

    /// 验证 CSR 矩阵结构的一致性
    pub fn validate(&self) -> Result<(), String> {
        // 检查 data 和 indices 长度一致
        if self.data.len() != self.indices.len() {
            return Err(format!(
                "SparseCSR: data length {} != indices length {}",
                self.data.len(),
                self.indices.len()
            ));
        }

        // 检查 indptr 长度
        if self.indptr.len() != self.n_rows + 1 {
            return Err(format!(
                "SparseCSR: indptr length {} != n_rows + 1 = {}",
                self.indptr.len(),
                self.n_rows + 1
            ));
        }

        // 检查 indptr 单调递增
        for i in 0..self.indptr.len() - 1 {
            if self.indptr[i] > self.indptr[i + 1] {
                return Err(format!(
                    "SparseCSR: indptr not monotonic at index {}: {} > {}",
                    i,
                    self.indptr[i],
                    self.indptr[i + 1]
                ));
            }
        }

        // 检查 indptr 最后一个元素等于 nnz
        if let Some(&last) = self.indptr.last() {
            if last != self.data.len() {
                return Err(format!(
                    "SparseCSR: indptr[-1] = {} != data.len() = {}",
                    last,
                    self.data.len()
                ));
            }
        }

        // 检查列索引在范围内
        for (i, &col_idx) in self.indices.iter().enumerate() {
            if col_idx >= self.n_cols {
                return Err(format!(
                    "SparseCSR: column index {} at position {} >= n_cols = {}",
                    col_idx, i, self.n_cols
                ));
            }
        }

        Ok(())
    }
}

/// CSC 格式稀疏矩阵（Compressed Sparse Column，R Matrix 默认）
///
/// 存储格式：
/// - data: 非零元素值
/// - indices: 每个非零元素的行索引
/// - indptr: 每列的起始位置（长度 = n_cols + 1）
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrixCSC {
    /// 非零元素值
    pub data: Vec<f64>,
    /// 行索引（对应每个非零元素）
    pub indices: Vec<usize>,
    /// 列指针（长度 = n_cols + 1）
    pub indptr: Vec<usize>,
    /// 行数
    pub n_rows: usize,
    /// 列数
    pub n_cols: usize,
}

impl SparseMatrixCSC {
    /// 创建新的 CSC 稀疏矩阵
    pub fn new(
        data: Vec<f64>,
        indices: Vec<usize>,
        indptr: Vec<usize>,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Self, String> {
        let matrix = Self {
            data,
            indices,
            indptr,
            n_rows,
            n_cols,
        };
        matrix.validate()?;
        Ok(matrix)
    }

    /// 创建空的 CSC 矩阵
    pub fn empty(n_rows: usize, n_cols: usize) -> Self {
        Self {
            data: Vec::new(),
            indices: Vec::new(),
            indptr: vec![0; n_cols + 1],
            n_rows,
            n_cols,
        }
    }

    /// 验证 CSC 矩阵结构的一致性
    pub fn validate(&self) -> Result<(), String> {
        // 检查 data 和 indices 长度一致
        if self.data.len() != self.indices.len() {
            return Err(format!(
                "SparseCSC: data length {} != indices length {}",
                self.data.len(),
                self.indices.len()
            ));
        }

        // 检查 indptr 长度
        if self.indptr.len() != self.n_cols + 1 {
            return Err(format!(
                "SparseCSC: indptr length {} != n_cols + 1 = {}",
                self.indptr.len(),
                self.n_cols + 1
            ));
        }

        // 检查 indptr 单调递增
        for i in 0..self.indptr.len() - 1 {
            if self.indptr[i] > self.indptr[i + 1] {
                return Err(format!(
                    "SparseCSC: indptr not monotonic at index {}: {} > {}",
                    i,
                    self.indptr[i],
                    self.indptr[i + 1]
                ));
            }
        }

        // 检查 indptr 最后一个元素等于 nnz
        if let Some(&last) = self.indptr.last() {
            if last != self.data.len() {
                return Err(format!(
                    "SparseCSC: indptr[-1] = {} != data.len() = {}",
                    last,
                    self.data.len()
                ));
            }
        }

        // 检查行索引在范围内
        for (i, &row_idx) in self.indices.iter().enumerate() {
            if row_idx >= self.n_rows {
                return Err(format!(
                    "SparseCSC: row index {} at position {} >= n_rows = {}",
                    row_idx, i, self.n_rows
                ));
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dense_matrix_creation() {
        // 2×3 矩阵
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let matrix = DenseMatrix::new(data, 2, 3).unwrap();

        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 3);
        assert_eq!(matrix.get(0, 0), Some(1.0));
        assert_eq!(matrix.get(1, 2), Some(6.0));
    }

    #[test]
    fn test_dense_matrix_validation() {
        // 维度不匹配
        let data = vec![1.0, 2.0, 3.0];
        let result = DenseMatrix::new(data, 2, 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_sparse_csr_creation() {
        // 简单的 2×3 CSR 矩阵
        // [[1.0, 0.0, 2.0],
        //  [0.0, 3.0, 0.0]]
        let data = vec![1.0, 2.0, 3.0];
        let indices = vec![0, 2, 1];
        let indptr = vec![0, 2, 3];

        let matrix = SparseMatrixCSR::new(data, indices, indptr, 2, 3).unwrap();
        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 3);
        assert_eq!(matrix.data.len(), 3);
    }

    #[test]
    fn test_sparse_csr_validation() {
        // data 和 indices 长度不匹配
        let data = vec![1.0, 2.0];
        let indices = vec![0, 1, 2];
        let indptr = vec![0, 2, 3];
        let result = SparseMatrixCSR::new(data, indices, indptr, 2, 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_sparse_csc_creation() {
        // 简单的 2×3 CSC 矩阵
        // [[1.0, 0.0, 2.0],
        //  [0.0, 3.0, 0.0]]
        let data = vec![1.0, 3.0, 2.0];
        let indices = vec![0, 1, 0];
        let indptr = vec![0, 1, 2, 3];

        let matrix = SparseMatrixCSC::new(data, indices, indptr, 2, 3).unwrap();
        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 3);
        assert_eq!(matrix.data.len(), 3);
    }

    #[test]
    fn test_expression_matrix_shape() {
        let dense = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        assert_eq!(dense.shape(), (10, 20));

        let csr = ExpressionMatrix::SparseCSR(SparseMatrixCSR::empty(100, 50));
        assert_eq!(csr.shape(), (100, 50));
    }

    #[test]
    fn test_empty_matrix() {
        // 0×0 矩阵
        let dense = DenseMatrix::zeros(0, 0);
        assert_eq!(dense.validate(), Ok(()));

        let csr = SparseMatrixCSR::empty(0, 0);
        assert_eq!(csr.validate(), Ok(()));
    }

    #[test]
    fn test_single_row_matrix() {
        // 1×N 矩阵
        let dense = DenseMatrix::zeros(1, 1000);
        assert_eq!(dense.validate(), Ok(()));

        let csr = SparseMatrixCSR::empty(1, 1000);
        assert_eq!(csr.validate(), Ok(()));
    }

    #[test]
    fn test_single_column_matrix() {
        // N×1 矩阵
        let dense = DenseMatrix::zeros(1000, 1);
        assert_eq!(dense.validate(), Ok(()));

        let csc = SparseMatrixCSC::empty(1000, 1);
        assert_eq!(csc.validate(), Ok(()));
    }

    #[test]
    fn test_lazy_matrix_creation() {
        // 创建延迟加载矩阵
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (1000, 500),
            true,
            Some(SparseFormat::CSR),
        );

        assert_eq!(lazy.shape(), (1000, 500));
        assert!(lazy.is_sparse);
        assert_eq!(lazy.sparse_format, Some(SparseFormat::CSR));
        assert!(!lazy.is_cached());
    }

    #[test]
    fn test_lazy_matrix_with_nnz() {
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (1000, 500),
            true,
            Some(SparseFormat::CSR),
        )
        .with_nnz(25000);

        assert_eq!(lazy.nnz, Some(25000));
    }

    #[test]
    fn test_lazy_matrix_cache() {
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (10, 5),
            false,
            None,
        );

        // 初始状态：未缓存
        assert!(!lazy.is_cached());
        assert!(lazy.get_cached().is_none());

        // 缓存数据
        let dense = DenseMatrix::zeros(10, 5);
        lazy.cache_dense(dense);

        // 验证缓存
        assert!(lazy.is_cached());
        let cached = lazy.get_cached();
        assert!(cached.is_some());
        assert!(matches!(cached.unwrap(), ExpressionMatrix::Dense(_)));

        // 清除缓存
        lazy.clear_cache();
        assert!(!lazy.is_cached());
    }

    #[test]
    fn test_lazy_matrix_memory_estimation() {
        // 稀疏矩阵
        let lazy_sparse = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (10000, 5000),
            true,
            Some(SparseFormat::CSR),
        )
        .with_nnz(250000);

        let mem = lazy_sparse.estimate_memory_usage();
        // 250000 * (8 + 8) + (10000 + 1) * 8 = 4,000,000 + 80,008 = 4,080,008
        assert!(mem > 4_000_000);

        // 稠密矩阵
        let lazy_dense = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (1000, 500),
            false,
            None,
        );

        let mem_dense = lazy_dense.estimate_memory_usage();
        // 1000 * 500 * 8 = 4,000,000
        assert_eq!(mem_dense, 4_000_000);
    }

    #[test]
    fn test_expression_matrix_lazy() {
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (100, 50),
            true,
            Some(SparseFormat::CSC),
        )
        .with_nnz(500);

        let expr = ExpressionMatrix::Lazy(lazy);

        assert_eq!(expr.shape(), (100, 50));
        assert!(expr.is_lazy());
        assert!(expr.is_sparse());
        assert_eq!(expr.nnz(), 500);
        assert!(expr.validate().is_ok());
    }

    #[test]
    fn test_expression_matrix_display() {
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (1000, 500),
            true,
            Some(SparseFormat::CSR),
        )
        .with_nnz(25000);

        let expr = ExpressionMatrix::Lazy(lazy);
        let display = format!("{}", expr);
        assert!(display.contains("LazyMatrix"));
        assert!(display.contains("1000"));
        assert!(display.contains("500"));
        assert!(display.contains("CSR"));
        assert!(display.contains("25000 nnz"));
    }

    #[test]
    fn test_cache_policy() {
        let lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from("test.h5ad"),
            "/X".to_string(),
            (100, 50),
            true,
            None,
        )
        .with_cache_policy(CachePolicy::LRU(1024 * 1024)); // 1 MB

        assert_eq!(lazy.cache_policy, CachePolicy::LRU(1024 * 1024));
    }

    // ==================== ChunkedMatrix Tests ====================

    #[test]
    fn test_chunked_matrix_creation() {
        let chunked = ChunkedMatrix::new(5000, 10000, 500, true, Some(SparseFormat::CSR));

        assert_eq!(chunked.chunk_size, 5000);
        assert_eq!(chunked.total_rows, 10000);
        assert_eq!(chunked.total_cols, 500);
        assert!(chunked.is_sparse);
        assert_eq!(chunked.sparse_format, Some(SparseFormat::CSR));
        assert_eq!(chunked.num_chunks(), 0);
    }

    #[test]
    fn test_chunked_matrix_from_dense() {
        // 创建 100×50 的稠密矩阵
        let data: Vec<f64> = (0..5000).map(|i| i as f64).collect();
        let dense = DenseMatrix::new(data, 100, 50).unwrap();
        let matrix = ExpressionMatrix::Dense(dense);

        // 分块，每块 30 行
        let chunked = ChunkedMatrix::from_matrix(matrix, 30).unwrap();

        assert_eq!(chunked.total_rows, 100);
        assert_eq!(chunked.total_cols, 50);
        assert_eq!(chunked.chunk_size, 30);
        assert_eq!(chunked.num_chunks(), 4); // 30 + 30 + 30 + 10 = 100
        assert!(!chunked.is_sparse);
    }

    #[test]
    fn test_chunked_matrix_from_csr() {
        // 创建简单的 CSR 矩阵
        // 10×5 矩阵，每行有 2 个非零元素
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0];

        for row in 0..10 {
            data.push((row * 2) as f64);
            data.push((row * 2 + 1) as f64);
            indices.push(row % 5);
            indices.push((row + 1) % 5);
            indptr.push(data.len());
        }

        let csr = SparseMatrixCSR::new(data, indices, indptr, 10, 5).unwrap();
        let matrix = ExpressionMatrix::SparseCSR(csr);

        // 分块，每块 3 行
        let chunked = ChunkedMatrix::from_matrix(matrix, 3).unwrap();

        assert_eq!(chunked.total_rows, 10);
        assert_eq!(chunked.total_cols, 5);
        assert_eq!(chunked.num_chunks(), 4); // 3 + 3 + 3 + 1 = 10
        assert!(chunked.is_sparse);
        assert_eq!(chunked.sparse_format, Some(SparseFormat::CSR));
    }

    #[test]
    fn test_chunked_matrix_roundtrip_dense() {
        // 创建 100×50 的稠密矩阵
        let data: Vec<f64> = (0..5000).map(|i| i as f64).collect();
        let dense = DenseMatrix::new(data.clone(), 100, 50).unwrap();
        let original = ExpressionMatrix::Dense(dense);

        // 分块
        let chunked = ChunkedMatrix::from_matrix(original.clone(), 30).unwrap();

        // 合并
        let merged = chunked.to_matrix().unwrap();

        // 验证
        if let ExpressionMatrix::Dense(m) = merged {
            assert_eq!(m.n_rows, 100);
            assert_eq!(m.n_cols, 50);
            assert_eq!(m.data, data);
        } else {
            panic!("Expected Dense matrix");
        }
    }

    #[test]
    fn test_chunked_matrix_roundtrip_csr() {
        // 创建简单的 CSR 矩阵
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0];

        for row in 0..10 {
            data.push((row * 2) as f64);
            data.push((row * 2 + 1) as f64);
            indices.push(row % 5);
            indices.push((row + 1) % 5);
            indptr.push(data.len());
        }

        let csr = SparseMatrixCSR::new(data.clone(), indices.clone(), indptr.clone(), 10, 5).unwrap();
        let original = ExpressionMatrix::SparseCSR(csr);

        // 分块
        let chunked = ChunkedMatrix::from_matrix(original.clone(), 3).unwrap();

        // 合并
        let merged = chunked.to_matrix().unwrap();

        // 验证
        if let ExpressionMatrix::SparseCSR(m) = merged {
            assert_eq!(m.n_rows, 10);
            assert_eq!(m.n_cols, 5);
            assert_eq!(m.data, data);
            assert_eq!(m.indices, indices);
            assert_eq!(m.indptr, indptr);
        } else {
            panic!("Expected SparseCSR matrix");
        }
    }

    #[test]
    fn test_chunked_matrix_add_chunk() {
        let mut chunked = ChunkedMatrix::new(5000, 10000, 500, false, None);

        // 添加有效的块
        let chunk = ExpressionMatrix::Dense(DenseMatrix::zeros(5000, 500));
        assert!(chunked.add_chunk(chunk).is_ok());
        assert_eq!(chunked.num_chunks(), 1);

        // 添加列数不匹配的块
        let bad_chunk = ExpressionMatrix::Dense(DenseMatrix::zeros(5000, 100));
        assert!(chunked.add_chunk(bad_chunk).is_err());

        // 添加行数超过 chunk_size 的块
        let big_chunk = ExpressionMatrix::Dense(DenseMatrix::zeros(6000, 500));
        assert!(chunked.add_chunk(big_chunk).is_err());
    }

    #[test]
    fn test_chunked_matrix_memory_estimation() {
        // 创建 100×50 的稠密矩阵
        let data: Vec<f64> = (0..5000).map(|i| i as f64).collect();
        let dense = DenseMatrix::new(data, 100, 50).unwrap();
        let matrix = ExpressionMatrix::Dense(dense);

        let chunked = ChunkedMatrix::from_matrix(matrix, 30).unwrap();

        // 估算内存
        let mem = chunked.estimate_memory_usage();
        // 100 * 50 * 8 = 40000 bytes
        assert_eq!(mem, 40000);

        // 估算单个块内存
        let chunk_mem = chunked.estimate_chunk_memory();
        // 30 * 50 * 8 = 12000 bytes
        assert_eq!(chunk_mem, 12000);
    }

    #[test]
    fn test_chunked_matrix_display() {
        let chunked = ChunkedMatrix::new(5000, 10000, 500, true, Some(SparseFormat::CSR));
        let display = format!("{}", chunked);
        
        assert!(display.contains("ChunkedMatrix"));
        assert!(display.contains("10000"));
        assert!(display.contains("500"));
        assert!(display.contains("5000"));
        assert!(display.contains("CSR"));
    }

    #[test]
    fn test_chunked_matrix_empty() {
        let chunked = ChunkedMatrix::new(5000, 0, 500, true, Some(SparseFormat::CSR));
        
        let merged = chunked.to_matrix().unwrap();
        assert_eq!(merged.shape(), (0, 500));
    }

    #[test]
    fn test_chunked_matrix_default_chunk_size() {
        assert_eq!(ChunkedMatrix::DEFAULT_CHUNK_SIZE, 5000);
    }
}
