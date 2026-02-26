//! SingleCellExperiment (.rds) 格式转换器

use std::path::Path;

use crate::error::CrossCellError;
use crate::ir::SingleCellData;
use crate::sce::{sce_rds_to_ir, write_sce_rds};

use super::FormatConverter;

/// SingleCellExperiment 格式转换器
///
/// 支持读写 SingleCellExperiment 对象（Bioconductor 生态系统的标准格式）。
///
/// 注意：SCE 和 Seurat 都使用 .rds 扩展名，需要通过内容检测来区分。
/// 可以使用 .sce.rds 扩展名来明确指定 SCE 格式。
pub struct SceConverter;

impl FormatConverter for SceConverter {
    fn name(&self) -> &str {
        "sce"
    }

    fn display_name(&self) -> &str {
        "SingleCellExperiment"
    }

    fn extensions(&self) -> &[&str] {
        // 使用 .sce.rds 作为明确的扩展名
        // 普通 .rds 文件需要通过内容检测
        &[".sce.rds"]
    }

    fn can_read(&self) -> bool {
        true
    }

    fn can_write(&self) -> bool {
        true
    }

    fn description(&self) -> &str {
        "SingleCellExperiment RDS format for Bioconductor"
    }

    fn read(&self, path: &Path) -> Result<SingleCellData, CrossCellError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path encoding".to_string()))?;

        sce_rds_to_ir(path_str)
            .map_err(|e| CrossCellError::InvalidFormat(format!("Failed to read SCE object: {}", e)))
    }

    fn write(&self, data: &SingleCellData, path: &Path) -> Result<(), CrossCellError> {
        let path_str = path
            .to_str()
            .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path encoding".to_string()))?;

        write_sce_rds(data, path_str).map_err(|e| {
            CrossCellError::InvalidFormat(format!("Failed to write SCE object: {}", e))
        })
    }

    fn detect(&self, path: &Path) -> bool {
        let path_str = path.to_string_lossy().to_lowercase();

        // 检查 .sce.rds 扩展名
        if path_str.ends_with(".sce.rds") {
            return true;
        }

        // 对于普通 .rds 文件，尝试读取并检测是否为 SCE
        if path_str.ends_with(".rds") && path.exists() {
            // 尝试解析 RDS 文件并检查类
            if let Ok(file) = crate::rds::parse_rds(path) {
                return crate::sce::is_sce_object(&file.object);
            }
        }

        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sce_converter_info() {
        let converter = SceConverter;

        assert_eq!(converter.name(), "sce");
        assert_eq!(converter.display_name(), "SingleCellExperiment");
        assert!(converter.can_read());
        assert!(converter.can_write());
        assert!(converter.extensions().contains(&".sce.rds"));
    }

    #[test]
    fn test_sce_detect_extension() {
        let converter = SceConverter;

        // .sce.rds 应该被检测为 SCE
        assert!(converter.detect(Path::new("data.sce.rds")));
        assert!(converter.detect(Path::new("path/to/data.SCE.RDS")));

        // 其他扩展名不应该被检测
        assert!(!converter.detect(Path::new("data.h5ad")));
        assert!(!converter.detect(Path::new("data.loom")));
    }

    #[test]
    fn test_sce_read_real_file() {
        let converter = SceConverter;
        let path = Path::new("tests/data/sce_minimal.rds");

        if path.exists() {
            let result = converter.read(path);
            assert!(result.is_ok(), "Should read SCE file successfully");

            let data = result.unwrap();
            assert!(data.metadata.n_cells > 0);
            assert!(data.metadata.n_genes > 0);
        }
    }
}
