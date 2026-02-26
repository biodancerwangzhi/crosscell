//! 磁盘支持的矩阵实现
//!
//! 提供超大数据集的磁盘存储支持，内存占用稳定在单个块大小。

use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::error::{CrossCellError, Result};
use crate::ir::expression::{ExpressionMatrix, SparseFormat, SparseMatrixCSC, SparseMatrixCSR};
use crate::storage::chunk_store::ChunkStore;

/// 磁盘支持的矩阵配置
#[derive(Debug, Clone)]
pub struct DiskBackedConfig {
    /// 临时目录
    pub temp_dir: PathBuf,
    /// 每块的行数（细胞数）
    pub chunk_size: usize,
    /// 是否启用压缩
    pub compression: bool,
    /// 是否在完成后清理临时文件
    pub cleanup: bool,
}

impl Default for DiskBackedConfig {
    fn default() -> Self {
        Self {
            temp_dir: std::env::temp_dir(),
            chunk_size: 5000,
            compression: true,
            cleanup: true,
        }
    }
}

impl DiskBackedConfig {
    /// 创建新配置
    pub fn new() -> Self {
        Self::default()
    }

    /// 设置临时目录
    pub fn with_temp_dir<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.temp_dir = path.as_ref().to_path_buf();
        self
    }

    /// 设置块大小
    pub fn with_chunk_size(mut self, size: usize) -> Self {
        self.chunk_size = size;
        self
    }

    /// 设置是否压缩
    pub fn with_compression(mut self, enabled: bool) -> Self {
        self.compression = enabled;
        self
    }

    /// 设置是否清理
    pub fn with_cleanup(mut self, enabled: bool) -> Self {
        self.cleanup = enabled;
        self
    }
}

/// 磁盘支持的矩阵
///
/// 将大矩阵分块存储到磁盘，支持流式处理超大数据集。
///
/// ## 使用场景
/// - 10M+ 细胞的超大数据集
/// - 内存受限的环境
/// - 需要稳定内存占用的批处理
///
/// ## 内存使用
/// - 峰值内存 ≈ 2 × chunk_size × n_cols × 8 bytes
/// - 例如：chunk_size=5000, n_cols=30000 → ~2.4 GB
#[derive(Debug)]
pub struct DiskBackedMatrix {
    /// 配置
    config: DiskBackedConfig,
    /// 分块存储
    store: ChunkStore,
    /// 矩阵总行数
    total_rows: usize,
    /// 矩阵总列数
    total_cols: usize,
    /// 稀疏格式
    sparse_format: SparseFormat,
}

impl DiskBackedMatrix {
    /// 创建新的磁盘支持矩阵
    pub fn create(
        total_rows: usize,
        total_cols: usize,
        sparse_format: SparseFormat,
        config: DiskBackedConfig,
    ) -> Result<Self> {
        // 创建临时文件
        let file_name = format!(
            "crosscell_{}_{}.chunks",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let path = config.temp_dir.join(file_name);

        let store = ChunkStore::create(
            &path,
            total_rows,
            total_cols,
            config.chunk_size,
            sparse_format,
            config.compression,
        )?;

        Ok(Self {
            config,
            store,
            total_rows,
            total_cols,
            sparse_format,
        })
    }

    /// 从现有存储打开
    pub fn open<P: AsRef<Path>>(path: P, config: DiskBackedConfig) -> Result<Self> {
        let store = ChunkStore::open(path)?;

        Ok(Self {
            total_rows: store.total_rows,
            total_cols: store.total_cols,
            sparse_format: store.sparse_format,
            config,
            store,
        })
    }

    /// 从内存矩阵创建磁盘支持矩阵（流式写入）
    pub fn from_matrix(matrix: &ExpressionMatrix, config: DiskBackedConfig) -> Result<Self> {
        let (total_rows, total_cols) = matrix.shape();
        let sparse_format = match matrix {
            ExpressionMatrix::SparseCSR(_) => SparseFormat::CSR,
            ExpressionMatrix::SparseCSC(_) => SparseFormat::CSC,
            _ => {
                return Err(CrossCellError::InvalidFormat(
                    "DiskBackedMatrix only supports sparse matrices".to_string(),
                ))
            }
        };

        let mut disk_matrix = Self::create(total_rows, total_cols, sparse_format, config)?;

        // 分块写入
        let chunk_size = disk_matrix.config.chunk_size;
        let num_chunks = (total_rows + chunk_size - 1) / chunk_size;

        log::info!(
            "Writing {} chunks to disk (chunk_size={})",
            num_chunks,
            chunk_size
        );

        for chunk_idx in 0..num_chunks {
            let start_row = chunk_idx * chunk_size;
            let end_row = std::cmp::min((chunk_idx + 1) * chunk_size, total_rows);

            match matrix {
                ExpressionMatrix::SparseCSR(m) => {
                    let chunk = extract_csr_rows(m, start_row, end_row)?;
                    disk_matrix.store.write_csr_chunk(chunk_idx, &chunk)?;
                }
                ExpressionMatrix::SparseCSC(m) => {
                    let chunk = extract_csc_rows(m, start_row, end_row)?;
                    disk_matrix.store.write_csc_chunk(chunk_idx, &chunk)?;
                }
                _ => unreachable!(),
            }

            if (chunk_idx + 1) % 10 == 0 || chunk_idx == num_chunks - 1 {
                log::info!(
                    "Progress: {}/{} chunks written ({:.1}%)",
                    chunk_idx + 1,
                    num_chunks,
                    (chunk_idx + 1) as f64 / num_chunks as f64 * 100.0
                );
            }
        }

        Ok(disk_matrix)
    }

    /// 获取矩阵形状
    pub fn shape(&self) -> (usize, usize) {
        (self.total_rows, self.total_cols)
    }

    /// 获取块数量
    pub fn num_chunks(&self) -> usize {
        self.store.num_chunks()
    }

    /// 获取块大小
    pub fn chunk_size(&self) -> usize {
        self.config.chunk_size
    }

    /// 获取稀疏格式
    pub fn sparse_format(&self) -> SparseFormat {
        self.sparse_format
    }

    /// 获取磁盘使用量（字节）
    pub fn disk_usage(&self) -> u64 {
        self.store.disk_usage()
    }

    /// 估算单个块的内存使用量（字节）
    pub fn estimate_chunk_memory(&self) -> usize {
        // 假设 95% 稀疏度
        let chunk_rows = std::cmp::min(self.config.chunk_size, self.total_rows);
        let estimated_nnz = (chunk_rows * self.total_cols) / 20;
        estimated_nnz * (8 + 8) + (chunk_rows + 1) * 8
    }

    /// 读取指定块（CSR 格式）
    pub fn read_csr_chunk(&self, chunk_idx: usize) -> Result<SparseMatrixCSR> {
        self.store.read_csr_chunk(chunk_idx)
    }

    /// 读取指定块（CSC 格式）
    pub fn read_csc_chunk(&self, chunk_idx: usize) -> Result<SparseMatrixCSC> {
        self.store.read_csc_chunk(chunk_idx)
    }

    /// 写入指定块（CSR 格式）
    pub fn write_csr_chunk(&mut self, chunk_idx: usize, chunk: &SparseMatrixCSR) -> Result<()> {
        self.store.write_csr_chunk(chunk_idx, chunk)
    }

    /// 写入指定块（CSC 格式）
    pub fn write_csc_chunk(&mut self, chunk_idx: usize, chunk: &SparseMatrixCSC) -> Result<()> {
        self.store.write_csc_chunk(chunk_idx, chunk)
    }

    /// 转换为内存矩阵（合并所有块）
    ///
    /// 警告：这会将整个矩阵加载到内存，可能导致 OOM
    pub fn to_matrix(&self) -> Result<ExpressionMatrix> {
        log::warn!(
            "Loading entire disk-backed matrix to memory ({} rows × {} cols)",
            self.total_rows,
            self.total_cols
        );

        match self.sparse_format {
            SparseFormat::CSR => self.merge_csr_chunks(),
            SparseFormat::CSC => self.merge_csc_chunks(),
        }
    }

    /// 流式处理每个块
    ///
    /// 这是处理超大数据集的推荐方式，内存占用稳定。
    pub fn process_chunks<F, T>(&self, mut processor: F) -> Result<Vec<T>>
    where
        F: FnMut(usize, ExpressionMatrix) -> Result<T>,
    {
        let num_chunks = self.num_chunks();
        let mut results = Vec::with_capacity(num_chunks);

        for chunk_idx in 0..num_chunks {
            let chunk = match self.sparse_format {
                SparseFormat::CSR => ExpressionMatrix::SparseCSR(self.read_csr_chunk(chunk_idx)?),
                SparseFormat::CSC => ExpressionMatrix::SparseCSC(self.read_csc_chunk(chunk_idx)?),
            };

            let result = processor(chunk_idx, chunk)?;
            results.push(result);
        }

        Ok(results)
    }

    /// 并行处理每个块
    pub fn process_chunks_parallel<F, T>(&self, processor: F) -> Result<Vec<T>>
    where
        F: Fn(usize, ExpressionMatrix) -> Result<T> + Sync,
        T: Send,
    {
        let num_chunks = self.num_chunks();
        let chunk_indices: Vec<usize> = (0..num_chunks).collect();

        chunk_indices
            .into_par_iter()
            .map(|chunk_idx| {
                let chunk = match self.sparse_format {
                    SparseFormat::CSR => {
                        ExpressionMatrix::SparseCSR(self.store.read_csr_chunk(chunk_idx)?)
                    }
                    SparseFormat::CSC => {
                        ExpressionMatrix::SparseCSC(self.store.read_csc_chunk(chunk_idx)?)
                    }
                };
                processor(chunk_idx, chunk)
            })
            .collect()
    }

    /// 获取存储路径
    pub fn path(&self) -> &Path {
        &self.store.path
    }

    // ========== 私有方法 ==========

    fn merge_csr_chunks(&self) -> Result<ExpressionMatrix> {
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for chunk_idx in 0..self.num_chunks() {
            let chunk = self.read_csr_chunk(chunk_idx)?;
            let offset = data.len();

            data.extend_from_slice(&chunk.data);
            indices.extend_from_slice(&chunk.indices);

            for &ptr in chunk.indptr.iter().skip(1) {
                indptr.push(ptr + offset);
            }
        }

        let matrix = SparseMatrixCSR::new(data, indices, indptr, self.total_rows, self.total_cols)
            .map_err(|e| CrossCellError::InvalidFormat(e))?;

        Ok(ExpressionMatrix::SparseCSR(matrix))
    }

    fn merge_csc_chunks(&self) -> Result<ExpressionMatrix> {
        // CSC 合并比较复杂，需要按列合并
        let mut col_data: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.total_cols];
        let mut current_row_offset = 0;

        for chunk_idx in 0..self.num_chunks() {
            let chunk = self.read_csc_chunk(chunk_idx)?;

            for col in 0..chunk.n_cols {
                let col_start = chunk.indptr[col];
                let col_end = chunk.indptr[col + 1];

                for idx in col_start..col_end {
                    let row = chunk.indices[idx] + current_row_offset;
                    let value = chunk.data[idx];
                    col_data[col].push((row, value));
                }
            }
            current_row_offset += chunk.n_rows;
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

        let matrix = SparseMatrixCSC::new(data, indices, indptr, self.total_rows, self.total_cols)
            .map_err(|e| CrossCellError::InvalidFormat(e))?;

        Ok(ExpressionMatrix::SparseCSC(matrix))
    }
}

impl Drop for DiskBackedMatrix {
    fn drop(&mut self) {
        if self.config.cleanup {
            if let Err(e) = std::fs::remove_file(&self.store.path) {
                log::warn!("Failed to cleanup temp file {:?}: {}", self.store.path, e);
            }
        }
    }
}

// ========== 辅助函数 ==========

/// 从 CSR 矩阵提取指定行范围
fn extract_csr_rows(
    matrix: &SparseMatrixCSR,
    start_row: usize,
    end_row: usize,
) -> Result<SparseMatrixCSR> {
    let chunk_rows = end_row - start_row;
    let start_ptr = matrix.indptr[start_row];
    let end_ptr = matrix.indptr[end_row];

    let data = matrix.data[start_ptr..end_ptr].to_vec();
    let indices = matrix.indices[start_ptr..end_ptr].to_vec();
    let indptr: Vec<usize> = matrix.indptr[start_row..=end_row]
        .iter()
        .map(|&p| p - start_ptr)
        .collect();

    SparseMatrixCSR::new(data, indices, indptr, chunk_rows, matrix.n_cols)
        .map_err(|e| CrossCellError::InvalidFormat(e))
}

/// 从 CSC 矩阵提取指定行范围
fn extract_csc_rows(
    matrix: &SparseMatrixCSC,
    start_row: usize,
    end_row: usize,
) -> Result<SparseMatrixCSC> {
    let chunk_rows = end_row - start_row;
    let mut new_data = Vec::new();
    let mut new_indices = Vec::new();
    let mut new_indptr = vec![0usize];

    for col in 0..matrix.n_cols {
        let col_start = matrix.indptr[col];
        let col_end = matrix.indptr[col + 1];

        for idx in col_start..col_end {
            let row = matrix.indices[idx];
            if row >= start_row && row < end_row {
                new_data.push(matrix.data[idx]);
                new_indices.push(row - start_row);
            }
        }
        new_indptr.push(new_data.len());
    }

    SparseMatrixCSC::new(new_data, new_indices, new_indptr, chunk_rows, matrix.n_cols)
        .map_err(|e| CrossCellError::InvalidFormat(e))
}
