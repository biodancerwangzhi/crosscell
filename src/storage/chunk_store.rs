//! 分块存储实现
//!
//! 将稀疏矩阵分块存储到磁盘，支持按需加载。

use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::error::{CrossCellError, Result};
use crate::ir::expression::{SparseFormat, SparseMatrixCSC, SparseMatrixCSR};

/// 块索引条目
#[derive(Debug, Clone)]
pub struct ChunkIndex {
    /// 块在文件中的偏移量
    pub offset: u64,
    /// 块的压缩大小（字节）
    pub compressed_size: u64,
    /// 块的原始大小（字节）
    pub original_size: u64,
    /// 块的行数
    pub n_rows: usize,
    /// 非零元素数量
    pub nnz: usize,
}

/// 分块存储
///
/// 将稀疏矩阵分块存储到单个文件中，支持随机访问。
#[derive(Debug)]
pub struct ChunkStore {
    /// 存储文件路径
    pub path: PathBuf,
    /// 矩阵总行数
    pub total_rows: usize,
    /// 矩阵总列数
    pub total_cols: usize,
    /// 每块的行数
    pub chunk_size: usize,
    /// 稀疏格式
    pub sparse_format: SparseFormat,
    /// 块索引
    pub indices: Vec<ChunkIndex>,
    /// 是否启用压缩
    pub compression: bool,
}

/// 文件头魔数
const MAGIC: &[u8; 8] = b"CROSSCEL";
/// 文件版本
const VERSION: u32 = 1;

impl ChunkStore {
    /// 创建新的分块存储
    pub fn create<P: AsRef<Path>>(
        path: P,
        total_rows: usize,
        total_cols: usize,
        chunk_size: usize,
        sparse_format: SparseFormat,
        compression: bool,
    ) -> Result<Self> {
        let path = path.as_ref().to_path_buf();

        // 创建文件并写入头部
        let mut file = BufWriter::new(File::create(&path)?);

        // 写入魔数
        file.write_all(MAGIC)?;

        // 写入版本
        file.write_all(&VERSION.to_le_bytes())?;

        // 写入元数据
        file.write_all(&(total_rows as u64).to_le_bytes())?;
        file.write_all(&(total_cols as u64).to_le_bytes())?;
        file.write_all(&(chunk_size as u64).to_le_bytes())?;

        // 写入稀疏格式 (0 = CSR, 1 = CSC)
        let format_byte: u8 = match sparse_format {
            SparseFormat::CSR => 0,
            SparseFormat::CSC => 1,
        };
        file.write_all(&[format_byte])?;

        // 写入压缩标志
        file.write_all(&[if compression { 1 } else { 0 }])?;

        // 预留索引区域大小（后续写入）
        let num_chunks = (total_rows + chunk_size - 1) / chunk_size;
        file.write_all(&(num_chunks as u64).to_le_bytes())?;

        // 预留索引空间（每个索引 40 字节）
        let index_size = num_chunks * 40;
        let zeros = vec![0u8; index_size];
        file.write_all(&zeros)?;

        file.flush()?;

        Ok(Self {
            path,
            total_rows,
            total_cols,
            chunk_size,
            sparse_format,
            indices: Vec::with_capacity(num_chunks),
            compression,
        })
    }

    /// 打开现有的分块存储
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let mut file = BufReader::new(File::open(&path)?);

        // 读取并验证魔数
        let mut magic = [0u8; 8];
        file.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(CrossCellError::InvalidFormat(
                "Invalid chunk store file: wrong magic number".to_string(),
            ));
        }

        // 读取版本
        let mut version_bytes = [0u8; 4];
        file.read_exact(&mut version_bytes)?;
        let version = u32::from_le_bytes(version_bytes);
        if version != VERSION {
            return Err(CrossCellError::InvalidFormat(format!(
                "Unsupported chunk store version: {}",
                version
            )));
        }

        // 读取元数据
        let mut buf8 = [0u8; 8];

        file.read_exact(&mut buf8)?;
        let total_rows = u64::from_le_bytes(buf8) as usize;

        file.read_exact(&mut buf8)?;
        let total_cols = u64::from_le_bytes(buf8) as usize;

        file.read_exact(&mut buf8)?;
        let chunk_size = u64::from_le_bytes(buf8) as usize;

        // 读取稀疏格式
        let mut format_byte = [0u8; 1];
        file.read_exact(&mut format_byte)?;
        let sparse_format = match format_byte[0] {
            0 => SparseFormat::CSR,
            1 => SparseFormat::CSC,
            _ => {
                return Err(CrossCellError::InvalidFormat(
                    "Invalid sparse format in chunk store".to_string(),
                ))
            }
        };

        // 读取压缩标志
        let mut compression_byte = [0u8; 1];
        file.read_exact(&mut compression_byte)?;
        let compression = compression_byte[0] != 0;

        // 读取块数量
        file.read_exact(&mut buf8)?;
        let num_chunks = u64::from_le_bytes(buf8) as usize;

        // 读取索引
        let mut indices = Vec::with_capacity(num_chunks);
        for _ in 0..num_chunks {
            file.read_exact(&mut buf8)?;
            let offset = u64::from_le_bytes(buf8);

            file.read_exact(&mut buf8)?;
            let compressed_size = u64::from_le_bytes(buf8);

            file.read_exact(&mut buf8)?;
            let original_size = u64::from_le_bytes(buf8);

            file.read_exact(&mut buf8)?;
            let n_rows = u64::from_le_bytes(buf8) as usize;

            file.read_exact(&mut buf8)?;
            let nnz = u64::from_le_bytes(buf8) as usize;

            // 只添加有效的索引（offset > 0 表示已写入）
            if offset > 0 {
                indices.push(ChunkIndex {
                    offset,
                    compressed_size,
                    original_size,
                    n_rows,
                    nnz,
                });
            }
        }

        Ok(Self {
            path,
            total_rows,
            total_cols,
            chunk_size,
            sparse_format,
            indices,
            compression,
        })
    }

    /// 写入一个 CSR 块
    pub fn write_csr_chunk(&mut self, chunk_idx: usize, matrix: &SparseMatrixCSR) -> Result<()> {
        if self.sparse_format != SparseFormat::CSR {
            return Err(CrossCellError::InvalidFormat(
                "Cannot write CSR chunk to CSC store".to_string(),
            ));
        }

        // 序列化块数据
        let serialized = self.serialize_csr(matrix)?;

        // 压缩（如果启用）
        let (data, compressed_size, original_size) = if self.compression {
            let compressed = lz4_flex::compress_prepend_size(&serialized);
            let compressed_size = compressed.len() as u64;
            let original_size = serialized.len() as u64;
            (compressed, compressed_size, original_size)
        } else {
            let size = serialized.len() as u64;
            (serialized, size, size)
        };

        // 打开文件追加写入
        let mut file = OpenOptions::new().read(true).write(true).open(&self.path)?;

        // 移动到文件末尾
        let offset = file.seek(SeekFrom::End(0))?;

        // 写入数据
        file.write_all(&data)?;

        // 创建索引条目
        let index = ChunkIndex {
            offset,
            compressed_size,
            original_size,
            n_rows: matrix.n_rows,
            nnz: matrix.data.len(),
        };

        // 更新索引区域
        self.write_index(&mut file, chunk_idx, &index)?;

        // 更新内存中的索引
        while self.indices.len() <= chunk_idx {
            self.indices.push(ChunkIndex {
                offset: 0,
                compressed_size: 0,
                original_size: 0,
                n_rows: 0,
                nnz: 0,
            });
        }
        self.indices[chunk_idx] = index;

        Ok(())
    }

    /// 写入一个 CSC 块
    pub fn write_csc_chunk(&mut self, chunk_idx: usize, matrix: &SparseMatrixCSC) -> Result<()> {
        if self.sparse_format != SparseFormat::CSC {
            return Err(CrossCellError::InvalidFormat(
                "Cannot write CSC chunk to CSR store".to_string(),
            ));
        }

        // 序列化块数据
        let serialized = self.serialize_csc(matrix)?;

        // 压缩（如果启用）
        let (data, compressed_size, original_size) = if self.compression {
            let compressed = lz4_flex::compress_prepend_size(&serialized);
            let compressed_size = compressed.len() as u64;
            let original_size = serialized.len() as u64;
            (compressed, compressed_size, original_size)
        } else {
            let size = serialized.len() as u64;
            (serialized, size, size)
        };

        // 打开文件追加写入
        let mut file = OpenOptions::new().read(true).write(true).open(&self.path)?;

        // 移动到文件末尾
        let offset = file.seek(SeekFrom::End(0))?;

        // 写入数据
        file.write_all(&data)?;

        // 创建索引条目
        let index = ChunkIndex {
            offset,
            compressed_size,
            original_size,
            n_rows: matrix.n_rows,
            nnz: matrix.data.len(),
        };

        // 更新索引区域
        self.write_index(&mut file, chunk_idx, &index)?;

        // 更新内存中的索引
        while self.indices.len() <= chunk_idx {
            self.indices.push(ChunkIndex {
                offset: 0,
                compressed_size: 0,
                original_size: 0,
                n_rows: 0,
                nnz: 0,
            });
        }
        self.indices[chunk_idx] = index;

        Ok(())
    }

    /// 读取一个 CSR 块
    pub fn read_csr_chunk(&self, chunk_idx: usize) -> Result<SparseMatrixCSR> {
        if self.sparse_format != SparseFormat::CSR {
            return Err(CrossCellError::InvalidFormat(
                "Cannot read CSR chunk from CSC store".to_string(),
            ));
        }

        if chunk_idx >= self.indices.len() {
            return Err(CrossCellError::InvalidFormat(format!(
                "Chunk index {} out of range",
                chunk_idx
            )));
        }

        let index = &self.indices[chunk_idx];

        // 读取数据
        let mut file = BufReader::new(File::open(&self.path)?);
        file.seek(SeekFrom::Start(index.offset))?;

        let mut data = vec![0u8; index.compressed_size as usize];
        file.read_exact(&mut data)?;

        // 解压缩（如果需要）
        let decompressed = if self.compression {
            lz4_flex::decompress_size_prepended(&data).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Decompression failed: {}", e))
            })?
        } else {
            data
        };

        // 反序列化
        self.deserialize_csr(&decompressed, index.n_rows, index.nnz)
    }

    /// 读取一个 CSC 块
    pub fn read_csc_chunk(&self, chunk_idx: usize) -> Result<SparseMatrixCSC> {
        if self.sparse_format != SparseFormat::CSC {
            return Err(CrossCellError::InvalidFormat(
                "Cannot read CSC chunk from CSR store".to_string(),
            ));
        }

        if chunk_idx >= self.indices.len() {
            return Err(CrossCellError::InvalidFormat(format!(
                "Chunk index {} out of range",
                chunk_idx
            )));
        }

        let index = &self.indices[chunk_idx];

        // 读取数据
        let mut file = BufReader::new(File::open(&self.path)?);
        file.seek(SeekFrom::Start(index.offset))?;

        let mut data = vec![0u8; index.compressed_size as usize];
        file.read_exact(&mut data)?;

        // 解压缩（如果需要）
        let decompressed = if self.compression {
            lz4_flex::decompress_size_prepended(&data).map_err(|e| {
                CrossCellError::InvalidFormat(format!("Decompression failed: {}", e))
            })?
        } else {
            data
        };

        // 反序列化
        self.deserialize_csc(&decompressed, index.n_rows, index.nnz)
    }

    /// 获取块数量
    pub fn num_chunks(&self) -> usize {
        (self.total_rows + self.chunk_size - 1) / self.chunk_size
    }

    /// 获取已写入的块数量
    pub fn written_chunks(&self) -> usize {
        self.indices.iter().filter(|i| i.offset > 0).count()
    }

    /// 估算磁盘使用量（字节）
    pub fn disk_usage(&self) -> u64 {
        self.indices.iter().map(|i| i.compressed_size).sum()
    }

    // ========== 私有方法 ==========

    /// 写入索引到文件
    fn write_index(&self, file: &mut File, chunk_idx: usize, index: &ChunkIndex) -> Result<()> {
        // 索引区域起始位置：魔数(8) + 版本(4) + 元数据(8*3) + 格式(1) + 压缩(1) + 块数(8) = 46
        let index_start = 46u64;
        let index_offset = index_start + (chunk_idx as u64) * 40;

        file.seek(SeekFrom::Start(index_offset))?;

        file.write_all(&index.offset.to_le_bytes())?;
        file.write_all(&index.compressed_size.to_le_bytes())?;
        file.write_all(&index.original_size.to_le_bytes())?;
        file.write_all(&(index.n_rows as u64).to_le_bytes())?;
        file.write_all(&(index.nnz as u64).to_le_bytes())?;

        Ok(())
    }

    /// 序列化 CSR 矩阵
    fn serialize_csr(&self, matrix: &SparseMatrixCSR) -> Result<Vec<u8>> {
        let mut buf = Vec::new();

        // 写入维度
        buf.extend_from_slice(&(matrix.n_rows as u64).to_le_bytes());
        buf.extend_from_slice(&(matrix.n_cols as u64).to_le_bytes());
        buf.extend_from_slice(&(matrix.data.len() as u64).to_le_bytes());

        // 写入 data
        for &val in &matrix.data {
            buf.extend_from_slice(&val.to_le_bytes());
        }

        // 写入 indices
        for &idx in &matrix.indices {
            buf.extend_from_slice(&(idx as u64).to_le_bytes());
        }

        // 写入 indptr
        for &ptr in &matrix.indptr {
            buf.extend_from_slice(&(ptr as u64).to_le_bytes());
        }

        Ok(buf)
    }

    /// 序列化 CSC 矩阵
    fn serialize_csc(&self, matrix: &SparseMatrixCSC) -> Result<Vec<u8>> {
        let mut buf = Vec::new();

        // 写入维度
        buf.extend_from_slice(&(matrix.n_rows as u64).to_le_bytes());
        buf.extend_from_slice(&(matrix.n_cols as u64).to_le_bytes());
        buf.extend_from_slice(&(matrix.data.len() as u64).to_le_bytes());

        // 写入 data
        for &val in &matrix.data {
            buf.extend_from_slice(&val.to_le_bytes());
        }

        // 写入 indices
        for &idx in &matrix.indices {
            buf.extend_from_slice(&(idx as u64).to_le_bytes());
        }

        // 写入 indptr
        for &ptr in &matrix.indptr {
            buf.extend_from_slice(&(ptr as u64).to_le_bytes());
        }

        Ok(buf)
    }

    /// 反序列化 CSR 矩阵
    fn deserialize_csr(&self, buf: &[u8], n_rows: usize, nnz: usize) -> Result<SparseMatrixCSR> {
        let mut offset = 0;

        // 读取维度
        let n_rows_read = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;
        let n_cols = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;
        let nnz_read = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;

        if n_rows_read != n_rows || nnz_read != nnz {
            return Err(CrossCellError::InvalidFormat(
                "Chunk metadata mismatch".to_string(),
            ));
        }

        // 读取 data
        let mut data = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let val = f64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap());
            data.push(val);
            offset += 8;
        }

        // 读取 indices
        let mut indices = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let idx = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
            indices.push(idx);
            offset += 8;
        }

        // 读取 indptr
        let mut indptr = Vec::with_capacity(n_rows + 1);
        for _ in 0..=n_rows {
            let ptr = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
            indptr.push(ptr);
            offset += 8;
        }

        SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols)
            .map_err(|e| CrossCellError::InvalidFormat(e))
    }

    /// 反序列化 CSC 矩阵
    fn deserialize_csc(&self, buf: &[u8], n_rows: usize, nnz: usize) -> Result<SparseMatrixCSC> {
        let mut offset = 0;

        // 读取维度
        let n_rows_read = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;
        let n_cols = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;
        let nnz_read = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
        offset += 8;

        if n_rows_read != n_rows || nnz_read != nnz {
            return Err(CrossCellError::InvalidFormat(
                "Chunk metadata mismatch".to_string(),
            ));
        }

        // 读取 data
        let mut data = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let val = f64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap());
            data.push(val);
            offset += 8;
        }

        // 读取 indices
        let mut indices = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let idx = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
            indices.push(idx);
            offset += 8;
        }

        // 读取 indptr
        let mut indptr = Vec::with_capacity(n_cols + 1);
        for _ in 0..=n_cols {
            let ptr = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap()) as usize;
            indptr.push(ptr);
            offset += 8;
        }

        SparseMatrixCSC::new(data, indices, indptr, n_rows, n_cols)
            .map_err(|e| CrossCellError::InvalidFormat(e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_chunk_store_csr_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.chunks");

        // 创建存储
        let mut store = ChunkStore::create(&path, 100, 50, 25, SparseFormat::CSR, true).unwrap();

        // 创建测试矩阵
        let matrix =
            SparseMatrixCSR::new(vec![1.0, 2.0, 3.0], vec![0, 2, 1], vec![0, 2, 3], 2, 50).unwrap();

        // 写入
        store.write_csr_chunk(0, &matrix).unwrap();

        // 重新打开并读取
        let store2 = ChunkStore::open(&path).unwrap();
        let matrix2 = store2.read_csr_chunk(0).unwrap();

        assert_eq!(matrix.data, matrix2.data);
        assert_eq!(matrix.indices, matrix2.indices);
        assert_eq!(matrix.indptr, matrix2.indptr);
    }

    #[test]
    fn test_chunk_store_csc_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.chunks");

        // 创建存储
        let mut store = ChunkStore::create(&path, 100, 50, 25, SparseFormat::CSC, true).unwrap();

        // 创建测试矩阵
        let matrix =
            SparseMatrixCSC::new(vec![1.0, 2.0, 3.0], vec![0, 1, 0], vec![0, 1, 2, 3], 100, 3)
                .unwrap();

        // 写入
        store.write_csc_chunk(0, &matrix).unwrap();

        // 重新打开并读取
        let store2 = ChunkStore::open(&path).unwrap();
        let matrix2 = store2.read_csc_chunk(0).unwrap();

        assert_eq!(matrix.data, matrix2.data);
        assert_eq!(matrix.indices, matrix2.indices);
        assert_eq!(matrix.indptr, matrix2.indptr);
    }
}
