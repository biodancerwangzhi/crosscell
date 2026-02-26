//! Seurat (.rds) 格式转换器

use std::path::Path;

use crate::error::CrossCellError;
use crate::ir::SingleCellData;
use crate::seurat::{seurat_rds_to_ir, write_seurat_rds};

use super::FormatConverter;

/// Seurat 格式转换器
///
/// 支持读写 .rds 文件（R 的序列化格式，用于 Seurat 对象）。
pub struct SeuratConverter;

impl FormatConverter for SeuratConverter {
    fn name(&self) -> &str {
        "seurat"
    }

    fn display_name(&self) -> &str {
        "Seurat"
    }

    fn extensions(&self) -> &[&str] {
        &[".rds"]
    }

    fn can_read(&self) -> bool {
        true
    }

    fn can_write(&self) -> bool {
        true
    }

    fn description(&self) -> &str {
        "Seurat RDS format for R"
    }

    fn read(&self, path: &Path) -> Result<SingleCellData, CrossCellError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path encoding".to_string()))?;

        seurat_rds_to_ir(path_str).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to read Seurat object: {}", e))
        })
    }

    fn write(&self, data: &SingleCellData, path: &Path) -> Result<(), CrossCellError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path encoding".to_string()))?;

        write_seurat_rds(data, path_str).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to write Seurat object: {}", e))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seurat_converter_info() {
        let converter = SeuratConverter;

        assert_eq!(converter.name(), "seurat");
        assert_eq!(converter.display_name(), "Seurat");
        assert!(converter.can_read());
        assert!(converter.can_write());
        assert!(converter.extensions().contains(&".rds"));
    }

    #[test]
    fn test_seurat_detect() {
        let converter = SeuratConverter;

        assert!(converter.detect(Path::new("data.rds")));
        assert!(converter.detect(Path::new("path/to/data.RDS")));
        assert!(!converter.detect(Path::new("data.h5ad")));
        assert!(!converter.detect(Path::new("data.loom")));
    }
}
