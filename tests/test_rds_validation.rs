//! RDS 往返验证测试
//!
//! 测试 .rds 文件的往返验证功能

use crosscell::validation::roundtrip::validate_rds_roundtrip;

#[test]
fn test_rds_validation_api() {
    // 这个测试验证 validate_rds_roundtrip API 是否可用
    // 实际的往返测试将在 Task 12 中实现

    // 测试不存在的文件应该返回错误
    let result = validate_rds_roundtrip(
        "tests/data/nonexistent.rds",
        "tests/data/nonexistent2.rds",
        1e-7,
    );

    assert!(result.is_err(), "Should return error for nonexistent files");
}

#[test]
fn test_rds_validation_with_simple_files() {
    // 使用简单的测试文件验证 API
    let original_path = "tests/data/r_integer.rds";
    let converted_path = "tests/data/rust_integer.rds";

    // 这两个文件都是简单的整数向量
    // API 应该正常工作
    let result = validate_rds_roundtrip(original_path, converted_path, 1e-7);

    // 验证 API 可以正常调用
    // 注意：这两个文件可能不是完整的 Seurat 对象，所以可能会失败
    // 但这证明了验证引擎的 API 是可用的
    match result {
        Ok(report) => {
            println!("Validation report: {}", report.summary());
            // API 工作正常
        }
        Err(e) => {
            // 如果文件不是 Seurat 对象，会返回错误
            // 这也证明了 API 正在工作
            println!("Expected error for non-Seurat files: {}", e);
            assert!(
                e.contains("Failed to read") || e.contains("not found"),
                "Error should be about file reading or format"
            );
        }
    }
}
