//! CrossCell - 单细胞转录组数据格式转换器
//!
//! 提供 H5AD 和 RDS 格式之间的高性能双向转换（纯 Rust 实现）。

pub mod anndata;
pub mod backend;
pub mod diagnostics;
pub mod error;
pub mod formats;
pub mod ir;
pub mod rds;
pub mod sce;
pub mod seurat;
pub mod sparse;
pub mod storage;
pub mod transform;
pub mod types;
pub mod validation;

// 重新导出核心类型
pub use anndata::{
    inspect_h5ad, read_h5ad, read_h5ad_partial, write_h5ad, H5adInfo, PartialLoadOptions,
};
pub use diagnostics::{detect_issues, DataCleaner, DiagnosticReport};
pub use error::{CrossCellError, Result};
pub use formats::{create_default_registry, FormatConverter, FormatInfo, FormatRegistry};
pub use ir::SingleCellData;
pub use rds::{RObject, RdsError, RdsFile};
pub use seurat::seurat_to_ir;
pub use sparse::{csc_to_csr, csr_to_csc};
pub use storage::{ChunkStore, DiskBackedConfig, DiskBackedMatrix};
pub use validation::ValidationReport;

// 提供便捷的 RDS 读写函数
use std::path::Path;

/// 读取 RDS 文件
pub fn read_rds(path: &Path) -> Result<RObject> {
    let file = rds::parse_rds(path).map_err(|e| {
        CrossCellError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            e.to_string(),
        ))
    })?;
    Ok(file.object)
}

/// 写入 RDS 文件
pub fn write_rds(path: &Path, obj: &RObject) -> Result<()> {
    let file = RdsFile::with_object(obj.clone());
    rds::write_rds(&file, path).map_err(|e| {
        CrossCellError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            e.to_string(),
        ))
    })
}
