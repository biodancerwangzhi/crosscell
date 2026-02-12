//! AnnData (.h5ad) 格式转换器

use std::path::Path;

use crate::anndata::{read_h5ad, write_h5ad};
use crate::error::CrossCellError;
use crate::ir::SingleCellData;

use super::FormatConverter;

/// AnnData 格式转换器
///
/// 支持读写 .h5ad 文件（AnnData 的 HDF5 格式）。
pub struct AnndataConverter;

impl FormatConverter for AnndataConverter {
    fn name(&self) -> &str {
        "anndata"
    }
    
    fn display_name(&self) -> &str {
        "AnnData"
    }
    
    fn extensions(&self) -> &[&str] {
        &[".h5ad"]
    }
    
    fn can_read(&self) -> bool {
        true
    }
    
    fn can_write(&self) -> bool {
        true
    }
    
    fn description(&self) -> &str {
        "AnnData HDF5 format for Python/Scanpy"
    }
    
    fn read(&self, path: &Path) -> Result<SingleCellData, CrossCellError> {
        read_h5ad(path).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to read H5AD file: {}", e))
        })
    }
    
    fn write(&self, data: &SingleCellData, path: &Path) -> Result<(), CrossCellError> {
        write_h5ad(data, path).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to write H5AD file: {}", e))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_anndata_converter_info() {
        let converter = AnndataConverter;
        
        assert_eq!(converter.name(), "anndata");
        assert_eq!(converter.display_name(), "AnnData");
        assert!(converter.can_read());
        assert!(converter.can_write());
        assert!(converter.extensions().contains(&".h5ad"));
    }
    
    #[test]
    fn test_anndata_detect() {
        let converter = AnndataConverter;
        
        assert!(converter.detect(Path::new("data.h5ad")));
        assert!(converter.detect(Path::new("path/to/data.H5AD")));
        assert!(!converter.detect(Path::new("data.rds")));
        assert!(!converter.detect(Path::new("data.loom")));
    }
}
