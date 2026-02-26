//! 往返验证
//!
//! 提供高级往返验证函数，用于验证数据转换的一致性。

use super::{compare::compare_single_cell_data, ValidationReport};
use crate::ir::SingleCellData;

/// 验证两个 SingleCellData 对象的一致性
///
/// # 参数
/// - `original`: 原始数据
/// - `converted`: 转换后的数据
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `ValidationReport`: 验证报告
///
/// # 示例
/// ```ignore
/// use crosscell::validation::roundtrip::validate_roundtrip;
///
/// let report = validate_roundtrip(&original_data, &converted_data, 1e-7);
/// if report.passed() {
///     println!("Roundtrip validation passed!");
/// } else {
///     println!("Validation failed: {}", report.summary());
/// }
/// ```
pub fn validate_roundtrip(
    original: &SingleCellData,
    converted: &SingleCellData,
    tolerance: f64,
) -> ValidationReport {
    compare_single_cell_data(original, converted, tolerance)
}

/// 验证 .h5ad 往返一致性
///
/// 读取两个 .h5ad 文件并比较它们的内容
///
/// # 参数
/// - `original_path`: 原始 .h5ad 文件路径
/// - `converted_path`: 转换后的 .h5ad 文件路径
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `Result<ValidationReport, String>`: 验证报告或错误信息
///
/// # 示例
/// ```ignore
/// use crosscell::validation::roundtrip::validate_h5ad_roundtrip;
///
/// match validate_h5ad_roundtrip("original.h5ad", "converted.h5ad", 1e-7) {
///     Ok(report) => {
///         if report.passed() {
///             println!("✓ H5AD roundtrip validation passed");
///         } else {
///             println!("✗ Validation failed:\n{}", report.detailed_report());
///         }
///     }
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn validate_h5ad_roundtrip(
    original_path: &str,
    converted_path: &str,
    tolerance: f64,
) -> Result<ValidationReport, String> {
    // 读取原始文件
    let original = crate::anndata::reader::read_h5ad(original_path)
        .map_err(|e| format!("Failed to read original file: {}", e))?;

    // 读取转换后的文件
    let converted = crate::anndata::reader::read_h5ad(converted_path)
        .map_err(|e| format!("Failed to read converted file: {}", e))?;

    // 比较
    Ok(validate_roundtrip(&original, &converted, tolerance))
}

/// 验证 .rds 往返一致性
///
/// 读取两个 .rds 文件并比较它们的内容
///
/// # 参数
/// - `original_path`: 原始 .rds 文件路径
/// - `converted_path`: 转换后的 .rds 文件路径
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `Result<ValidationReport, String>`: 验证报告或错误信息
///
/// # 示例
/// ```ignore
/// use crosscell::validation::roundtrip::validate_rds_roundtrip;
///
/// match validate_rds_roundtrip("original.rds", "converted.rds", 1e-7) {
///     Ok(report) => {
///         if report.passed() {
///             println!("✓ RDS roundtrip validation passed");
///         } else {
///             println!("✗ Validation failed:\n{}", report.detailed_report());
///         }
///     }
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn validate_rds_roundtrip(
    original_path: &str,
    converted_path: &str,
    tolerance: f64,
) -> Result<ValidationReport, String> {
    // 读取原始文件
    let original = crate::seurat::seurat_to_ir::seurat_rds_to_ir(original_path)
        .map_err(|e| format!("Failed to read original file: {}", e))?;

    // 读取转换后的文件
    let converted = crate::seurat::seurat_to_ir::seurat_rds_to_ir(converted_path)
        .map_err(|e| format!("Failed to read converted file: {}", e))?;

    // 比较
    Ok(validate_roundtrip(&original, &converted, tolerance))
}

/// 验证跨格式转换一致性
///
/// 比较 .h5ad 和 .rds 文件的内容
///
/// # 参数
/// - `h5ad_path`: .h5ad 文件路径
/// - `rds_path`: .rds 文件路径
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `Result<ValidationReport, String>`: 验证报告或错误信息
///
/// # 示例
/// ```ignore
/// use crosscell::validation::roundtrip::validate_cross_format;
///
/// match validate_cross_format("data.h5ad", "data.rds", 1e-7) {
///     Ok(report) => {
///         if report.passed() {
///             println!("✓ Cross-format validation passed");
///         } else {
///             println!("✗ Validation failed:\n{}", report.detailed_report());
///         }
///     }
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn validate_cross_format(
    h5ad_path: &str,
    rds_path: &str,
    tolerance: f64,
) -> Result<ValidationReport, String> {
    // 读取 .h5ad 文件
    let h5ad_data = crate::anndata::reader::read_h5ad(h5ad_path)
        .map_err(|e| format!("Failed to read H5AD file: {}", e))?;

    // 读取 .rds 文件
    let rds_data = crate::seurat::seurat_to_ir::seurat_rds_to_ir(rds_path)
        .map_err(|e| format!("Failed to read RDS file: {}", e))?;

    // 比较
    Ok(validate_roundtrip(&h5ad_data, &rds_data, tolerance))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::*;

    #[test]
    fn test_validate_roundtrip_identical() {
        // 创建两个相同的数据集
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta = DataFrame::empty(10);
        let gene_meta = DataFrame::empty(20);
        let metadata = DatasetMetadata::new(10, 20, "test".to_string());

        let data1 = SingleCellData::new(
            expr.clone(),
            cell_meta.clone(),
            gene_meta.clone(),
            metadata.clone(),
        )
        .unwrap();

        let data2 = SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap();

        let report = validate_roundtrip(&data1, &data2, 1e-7);
        assert!(report.passed(), "Identical data should pass validation");
    }

    #[test]
    fn test_validate_roundtrip_dimension_mismatch() {
        // 创建维度不匹配的数据集
        let expr1 = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta1 = DataFrame::empty(10);
        let gene_meta1 = DataFrame::empty(20);
        let metadata1 = DatasetMetadata::new(10, 20, "test".to_string());

        let data1 = SingleCellData::new(expr1, cell_meta1, gene_meta1, metadata1).unwrap();

        let expr2 = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 25));
        let cell_meta2 = DataFrame::empty(10);
        let gene_meta2 = DataFrame::empty(25);
        let metadata2 = DatasetMetadata::new(10, 25, "test".to_string());

        let data2 = SingleCellData::new(expr2, cell_meta2, gene_meta2, metadata2).unwrap();

        let report = validate_roundtrip(&data1, &data2, 1e-7);
        assert!(
            !report.passed(),
            "Dimension mismatch should fail validation"
        );
    }
}
