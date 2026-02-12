// Backend abstraction layer for CrossCell
// Inspired by anndata-rs backend design
// Provides unified interface for different storage formats (HDF5, RDS, Memory)

use std::path::{Path, PathBuf};
use crate::error::Result;

/// Scalar data types supported by backends
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarType {
    I8,
    I16,
    I32,
    I64,
    U8,
    U16,
    U32,
    U64,
    F32,
    F64,
    Bool,
    String,
}

/// Selection information for slicing datasets
#[derive(Debug, Clone)]
pub struct SelectInfo {
    pub start: usize,
    pub end: usize,
    pub step: usize,
}

impl SelectInfo {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end, step: 1 }
    }
    
    pub fn with_step(start: usize, end: usize, step: usize) -> Self {
        Self { start, end, step }
    }
}

/// Dynamic array type for backend operations
#[derive(Debug, Clone)]
pub enum DynArray {
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
    Bool(Vec<bool>),
    String(Vec<String>),
}

/// Generic value type for attributes
#[derive(Debug, Clone)]
pub enum Value {
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    Bool(bool),
    String(String),
    Array(DynArray),
}

/// Main storage backend trait
/// 
/// This trait defines the interface for different storage backends.
/// Each backend must implement this trait to provide storage operations.
pub trait StorageBackend: 'static + Sized {
    /// Backend name (e.g., "hdf5", "rds", "memory")
    const NAME: &'static str;
    
    /// Store type (file handle)
    type Store: StoreOp<Self> + GroupOp<Self>;
    
    /// Group type (container for datasets)
    type Group: GroupOp<Self> + AttributeOp<Self>;
    
    /// Dataset type (actual data storage)
    type Dataset: DatasetOp<Self> + AttributeOp<Self>;
    
    /// Create a new store at the given path
    fn new<P: AsRef<Path>>(path: P) -> Result<Self::Store>;
    
    /// Open an existing store in read-only mode
    fn open<P: AsRef<Path>>(path: P) -> Result<Self::Store>;
    
    /// Open an existing store in read-write mode
    fn open_rw<P: AsRef<Path>>(path: P) -> Result<Self::Store>;
}

/// Store operations (file-level operations)
pub trait StoreOp<B: StorageBackend> {
    /// Get the file path
    fn filename(&self) -> PathBuf;
    
    /// Close the store and flush all changes
    fn close(self) -> Result<()>;
}

/// Group operations (container operations)
pub trait GroupOp<B: StorageBackend> {
    /// List all items in this group
    fn list(&self) -> Result<Vec<String>>;
    
    /// Create a new subgroup
    fn new_group(&self, name: &str) -> Result<B::Group>;
    
    /// Open an existing subgroup
    fn open_group(&self, name: &str) -> Result<B::Group>;
    
    /// Create a new dataset
    fn new_dataset(&self, name: &str, shape: &[usize], dtype: ScalarType) -> Result<B::Dataset>;
    
    /// Open an existing dataset
    fn open_dataset(&self, name: &str) -> Result<B::Dataset>;
    
    /// Check if an item exists
    fn exists(&self, name: &str) -> Result<bool>;
    
    /// Delete an item
    fn delete(&self, name: &str) -> Result<()>;
}

/// Dataset operations (data read/write)
pub trait DatasetOp<B: StorageBackend> {
    /// Get the data type
    fn dtype(&self) -> Result<ScalarType>;
    
    /// Get the shape (dimensions)
    fn shape(&self) -> Vec<usize>;
    
    /// Read a slice of data
    fn read_slice(&self, selection: &[SelectInfo]) -> Result<DynArray>;
    
    /// Write a slice of data
    fn write_slice(&self, data: &DynArray, selection: &[SelectInfo]) -> Result<()>;
    
    /// Read all data
    fn read_all(&self) -> Result<DynArray> {
        let shape = self.shape();
        let selection: Vec<SelectInfo> = shape.iter()
            .map(|&dim| SelectInfo::new(0, dim))
            .collect();
        self.read_slice(&selection)
    }
    
    /// Write all data
    fn write_all(&self, data: &DynArray) -> Result<()> {
        let shape = self.shape();
        let selection: Vec<SelectInfo> = shape.iter()
            .map(|&dim| SelectInfo::new(0, dim))
            .collect();
        self.write_slice(data, &selection)
    }
}

/// 分块读取器
/// 
/// 用于按块读取大型矩阵，降低内存占用
pub struct ChunkedReader {
    /// 矩阵总维度
    pub total_shape: (usize, usize),
    /// 每块的大小
    pub chunk_size: (usize, usize),
    /// 当前块的索引 (row_chunk, col_chunk)
    current_chunk: (usize, usize),
    /// 总块数 (row_chunks, col_chunks)
    total_chunks: (usize, usize),
}

impl ChunkedReader {
    /// 创建新的分块读取器
    pub fn new(total_shape: (usize, usize), chunk_size: (usize, usize)) -> Self {
        let row_chunks = (total_shape.0 + chunk_size.0 - 1) / chunk_size.0.max(1);
        let col_chunks = (total_shape.1 + chunk_size.1 - 1) / chunk_size.1.max(1);
        
        Self {
            total_shape,
            chunk_size,
            current_chunk: (0, 0),
            total_chunks: (row_chunks.max(1), col_chunks.max(1)),
        }
    }

    /// 获取下一个块的选择范围
    pub fn next_selection(&mut self) -> Option<ChunkSelection> {
        if self.current_chunk.0 >= self.total_chunks.0 {
            return None;
        }

        let row_start = self.current_chunk.0 * self.chunk_size.0;
        let row_end = ((self.current_chunk.0 + 1) * self.chunk_size.0).min(self.total_shape.0);
        let col_start = self.current_chunk.1 * self.chunk_size.1;
        let col_end = ((self.current_chunk.1 + 1) * self.chunk_size.1).min(self.total_shape.1);

        let selection = ChunkSelection {
            row_range: (row_start, row_end),
            col_range: (col_start, col_end),
            chunk_index: self.current_chunk,
        };

        // 移动到下一个块
        self.current_chunk.1 += 1;
        if self.current_chunk.1 >= self.total_chunks.1 {
            self.current_chunk.1 = 0;
            self.current_chunk.0 += 1;
        }

        Some(selection)
    }

    /// 重置到第一个块
    pub fn reset(&mut self) {
        self.current_chunk = (0, 0);
    }

    /// 获取总块数
    pub fn total_chunk_count(&self) -> usize {
        self.total_chunks.0 * self.total_chunks.1
    }

    /// 检查是否还有更多块
    pub fn has_more(&self) -> bool {
        self.current_chunk.0 < self.total_chunks.0
    }
}

/// 块选择信息
#[derive(Debug, Clone)]
pub struct ChunkSelection {
    /// 行范围 (start, end)
    pub row_range: (usize, usize),
    /// 列范围 (start, end)
    pub col_range: (usize, usize),
    /// 块索引 (row_chunk, col_chunk)
    pub chunk_index: (usize, usize),
}

impl ChunkSelection {
    /// 获取块的形状
    pub fn shape(&self) -> (usize, usize) {
        (
            self.row_range.1 - self.row_range.0,
            self.col_range.1 - self.col_range.0,
        )
    }

    /// 转换为 SelectInfo 数组
    pub fn to_select_info(&self) -> Vec<SelectInfo> {
        vec![
            SelectInfo::new(self.row_range.0, self.row_range.1),
            SelectInfo::new(self.col_range.0, self.col_range.1),
        ]
    }
}

/// Attribute operations (metadata)
pub trait AttributeOp<B: StorageBackend> {
    /// Get an attribute value
    fn get_attr(&self, name: &str) -> Result<Value>;
    
    /// Set an attribute value
    fn set_attr(&mut self, name: &str, value: &Value) -> Result<()>;
    
    /// List all attribute names
    fn list_attrs(&self) -> Result<Vec<String>>;
    
    /// Check if an attribute exists
    fn has_attr(&self, name: &str) -> Result<bool>;
}

// Re-export submodules
pub mod hdf5;
pub mod rds;
pub mod memory;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunked_reader_basic() {
        // 10x10 矩阵，每块 3x3
        let mut reader = ChunkedReader::new((10, 10), (3, 3));
        
        // 应该有 4x4 = 16 个块
        assert_eq!(reader.total_chunk_count(), 16);
        assert!(reader.has_more());
        
        // 第一个块
        let chunk = reader.next_selection().unwrap();
        assert_eq!(chunk.row_range, (0, 3));
        assert_eq!(chunk.col_range, (0, 3));
        assert_eq!(chunk.chunk_index, (0, 0));
        assert_eq!(chunk.shape(), (3, 3));
    }

    #[test]
    fn test_chunked_reader_all_chunks() {
        let mut reader = ChunkedReader::new((10, 10), (3, 3));
        
        let mut count = 0;
        while let Some(chunk) = reader.next_selection() {
            count += 1;
            // 验证范围有效
            assert!(chunk.row_range.0 < chunk.row_range.1);
            assert!(chunk.col_range.0 < chunk.col_range.1);
            assert!(chunk.row_range.1 <= 10);
            assert!(chunk.col_range.1 <= 10);
        }
        
        assert_eq!(count, 16);
        assert!(!reader.has_more());
    }

    #[test]
    fn test_chunked_reader_last_chunk() {
        // 10x10 矩阵，每块 3x3
        // 最后一个块应该是 (9, 10) x (9, 10)，即 1x1
        let mut reader = ChunkedReader::new((10, 10), (3, 3));
        
        let mut last_chunk = None;
        while let Some(chunk) = reader.next_selection() {
            last_chunk = Some(chunk);
        }
        
        let last = last_chunk.unwrap();
        assert_eq!(last.row_range, (9, 10));
        assert_eq!(last.col_range, (9, 10));
        assert_eq!(last.shape(), (1, 1));
    }

    #[test]
    fn test_chunked_reader_exact_fit() {
        // 9x9 矩阵，每块 3x3，正好整除
        let mut reader = ChunkedReader::new((9, 9), (3, 3));
        
        assert_eq!(reader.total_chunk_count(), 9);
        
        let mut count = 0;
        while let Some(chunk) = reader.next_selection() {
            count += 1;
            // 每个块都应该是 3x3
            assert_eq!(chunk.shape(), (3, 3));
        }
        
        assert_eq!(count, 9);
    }

    #[test]
    fn test_chunked_reader_single_chunk() {
        // 小矩阵，一个块就够了
        let mut reader = ChunkedReader::new((5, 5), (10, 10));
        
        assert_eq!(reader.total_chunk_count(), 1);
        
        let chunk = reader.next_selection().unwrap();
        assert_eq!(chunk.row_range, (0, 5));
        assert_eq!(chunk.col_range, (0, 5));
        assert_eq!(chunk.shape(), (5, 5));
        
        assert!(reader.next_selection().is_none());
    }

    #[test]
    fn test_chunked_reader_reset() {
        let mut reader = ChunkedReader::new((10, 10), (5, 5));
        
        // 读取所有块
        while reader.next_selection().is_some() {}
        assert!(!reader.has_more());
        
        // 重置
        reader.reset();
        assert!(reader.has_more());
        
        // 第一个块应该又是 (0, 0)
        let chunk = reader.next_selection().unwrap();
        assert_eq!(chunk.chunk_index, (0, 0));
    }

    #[test]
    fn test_chunked_reader_row_only() {
        // 只按行分块
        let mut reader = ChunkedReader::new((100, 50), (10, 50));
        
        assert_eq!(reader.total_chunk_count(), 10);
        
        let chunk = reader.next_selection().unwrap();
        assert_eq!(chunk.row_range, (0, 10));
        assert_eq!(chunk.col_range, (0, 50));
    }

    #[test]
    fn test_chunk_selection_to_select_info() {
        let selection = ChunkSelection {
            row_range: (10, 20),
            col_range: (5, 15),
            chunk_index: (1, 0),
        };
        
        let select_info = selection.to_select_info();
        assert_eq!(select_info.len(), 2);
        assert_eq!(select_info[0].start, 10);
        assert_eq!(select_info[0].end, 20);
        assert_eq!(select_info[1].start, 5);
        assert_eq!(select_info[1].end, 15);
    }

    #[test]
    fn test_chunked_reader_empty_matrix() {
        // 空矩阵
        let mut reader = ChunkedReader::new((0, 0), (10, 10));
        
        // 应该至少有一个块（即使是空的）
        assert_eq!(reader.total_chunk_count(), 1);
    }

    #[test]
    fn test_chunked_reader_large_matrix() {
        // 大矩阵：100000 x 30000，每块 10000 x 5000
        let mut reader = ChunkedReader::new((100000, 30000), (10000, 5000));
        
        // 10 x 6 = 60 个块
        assert_eq!(reader.total_chunk_count(), 60);
        
        let mut count = 0;
        while reader.next_selection().is_some() {
            count += 1;
        }
        assert_eq!(count, 60);
    }
}
