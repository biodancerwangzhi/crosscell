//! 磁盘支持的矩阵存储模块
//!
//! 提供超大数据集（10M+ 细胞）的磁盘存储支持，借鉴 BPCells 设计。
//!
//! ## 核心特性
//! - 内存映射文件 (mmap) 实现零拷贝读取
//! - 分块存储降低内存峰值
//! - LZ4 压缩减少磁盘占用
//! - 流式处理支持超大数据集

pub mod chunk_store;
pub mod disk_backed;
pub mod mmap;

pub use chunk_store::ChunkStore;
pub use disk_backed::{DiskBackedConfig, DiskBackedMatrix};
pub use mmap::MmapReader;
