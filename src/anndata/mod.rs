//! AnnData (.h5ad) 文件格式支持
//!
//! 实现 HDF5 格式的 AnnData 文件读写功能。
//!
//! AnnData 文件结构：
//! - /X: 主表达矩阵（稀疏 CSR 或稠密）
//! - /obs: 细胞元数据
//! - /var: 基因元数据
//! - /obsm: 降维嵌入（PCA, UMAP 等）
//! - /layers: 多层表达矩阵
//! - /obsp: 细胞-细胞成对矩阵
//! - /varp: 基因-基因成对矩阵
//! - /uns: 非结构化数据

pub mod reader;
pub mod streaming;
pub mod writer;

pub use reader::{inspect_h5ad, read_h5ad, read_h5ad_partial, H5adInfo, PartialLoadOptions};
pub use streaming::{
    streaming_convert, streaming_h5ad_to_rds, streaming_rds_to_h5ad, DataChunk,
    StreamingH5adReader, StreamingH5adWriter, StreamingMetadata, StreamingRdsReader,
    DEFAULT_CHUNK_SIZE,
};
pub use writer::write_h5ad;

use thiserror::Error;

/// AnnData 错误类型
#[derive(Error, Debug)]
pub enum AnnDataError {
    #[error("HDF5 error: {0}")]
    Hdf5Error(#[from] hdf5::Error),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Invalid format: {0}")]
    InvalidFormat(String),

    #[error("Missing required field: {0}")]
    MissingField(String),

    #[error("Dimension mismatch: {0}")]
    DimensionMismatch(String),

    #[error("Unsupported data type: {0}")]
    UnsupportedType(String),
}

pub type Result<T> = std::result::Result<T, AnnDataError>;
