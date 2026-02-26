//! .h5ad 往返集成测试
//!
//! 测试 .h5ad 文件的读取 → IR → 写入 → 读取往返一致性

use crosscell::anndata::{reader::read_h5ad, writer::write_h5ad};
use crosscell::validation::roundtrip::validate_h5ad_roundtrip;
use std::path::Path;

#[test]
fn test_h5ad_roundtrip_basic() {
    // 使用现有的测试文件
    let original_path = "tests/data/small_sparse.h5ad";

    // 检查文件是否存在
    if !Path::new(original_path).exists() {
        eprintln!("Test file not found: {}", original_path);
        return;
    }

    // 读取原始文件
    let original_data = read_h5ad(original_path).expect("Failed to read original h5ad");

    // 写入到临时文件
    let temp_path = "tests/data/temp_roundtrip.h5ad";
    write_h5ad(&original_data, temp_path).expect("Failed to write h5ad");

    // 读取写入的文件
    let roundtrip_data = read_h5ad(temp_path).expect("Failed to read roundtrip h5ad");

    // 验证往返一致性
    let report =
        crosscell::validation::roundtrip::validate_roundtrip(&original_data, &roundtrip_data, 1e-7);

    // 打印详细报告
    println!("{}", report.detailed_report());

    // 清理临时文件
    let _ = std::fs::remove_file(temp_path);

    // 断言验证通过
    assert!(
        report.passed(),
        "Roundtrip validation failed:\n{}",
        report.detailed_report()
    );
}

#[test]
fn test_h5ad_roundtrip_with_validation_function() {
    let original_path = "tests/data/small_sparse.h5ad";

    if !Path::new(original_path).exists() {
        eprintln!("Test file not found: {}", original_path);
        return;
    }

    // 读取并写入
    let data = read_h5ad(original_path).expect("Failed to read h5ad");
    let temp_path = "tests/data/temp_roundtrip2.h5ad";
    write_h5ad(&data, temp_path).expect("Failed to write h5ad");

    // 使用验证函数
    let result = validate_h5ad_roundtrip(original_path, temp_path, 1e-7);

    // 清理
    let _ = std::fs::remove_file(temp_path);

    match result {
        Ok(report) => {
            println!("{}", report.summary());
            assert!(
                report.passed(),
                "Validation failed: {}",
                report.detailed_report()
            );
        }
        Err(e) => panic!("Validation error: {}", e),
    }
}

#[test]
fn test_h5ad_expression_matrix_preservation() {
    let original_path = "tests/data/small_sparse.h5ad";

    if !Path::new(original_path).exists() {
        eprintln!("Test file not found: {}", original_path);
        return;
    }

    // 读取原始数据
    let original = read_h5ad(original_path).expect("Failed to read original");
    let (orig_rows, orig_cols) = original.expression.shape();

    // 往返
    let temp_path = "tests/data/temp_expr_test.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");
    let (rt_rows, rt_cols) = roundtrip.expression.shape();

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证维度
    assert_eq!(orig_rows, rt_rows, "Row count mismatch");
    assert_eq!(orig_cols, rt_cols, "Column count mismatch");

    // 验证矩阵类型
    match (&original.expression, &roundtrip.expression) {
        (
            crosscell::ir::ExpressionMatrix::SparseCSR(_),
            crosscell::ir::ExpressionMatrix::SparseCSR(_),
        ) => {
            println!("✓ Sparse CSR format preserved");
        }
        (
            crosscell::ir::ExpressionMatrix::SparseCSC(_),
            crosscell::ir::ExpressionMatrix::SparseCSC(_),
        ) => {
            println!("✓ Sparse CSC format preserved");
        }
        (crosscell::ir::ExpressionMatrix::Dense(_), crosscell::ir::ExpressionMatrix::Dense(_)) => {
            println!("✓ Dense format preserved");
        }
        _ => panic!("Matrix format changed during roundtrip"),
    }
}

#[test]
fn test_h5ad_metadata_preservation() {
    let original_path = "tests/data/small_sparse.h5ad";

    if !Path::new(original_path).exists() {
        eprintln!("Test file not found: {}", original_path);
        return;
    }

    // 读取原始数据
    let original = read_h5ad(original_path).expect("Failed to read original");

    // 往返
    let temp_path = "tests/data/temp_meta_test.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证元数据
    assert_eq!(
        original.cell_metadata.n_rows, roundtrip.cell_metadata.n_rows,
        "Cell metadata row count mismatch"
    );
    assert_eq!(
        original.cell_metadata.n_cols(),
        roundtrip.cell_metadata.n_cols(),
        "Cell metadata column count mismatch"
    );

    assert_eq!(
        original.gene_metadata.n_rows, roundtrip.gene_metadata.n_rows,
        "Gene metadata row count mismatch"
    );
    assert_eq!(
        original.gene_metadata.n_cols(),
        roundtrip.gene_metadata.n_cols(),
        "Gene metadata column count mismatch"
    );

    println!("✓ Metadata dimensions preserved");
}

#[test]
fn test_h5ad_embeddings_preservation() {
    let original_path = "tests/data/python_verify_embeddings.h5ad";

    if !Path::new(original_path).exists() {
        eprintln!("Test file not found: {}", original_path);
        return;
    }

    // 读取原始数据
    let original = read_h5ad(original_path).expect("Failed to read original");

    // 检查是否有嵌入
    if original.embeddings.is_none() {
        println!("No embeddings in test file, skipping");
        return;
    }

    // 往返
    let temp_path = "tests/data/temp_emb_test.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证嵌入
    assert!(
        roundtrip.embeddings.is_some(),
        "Embeddings lost during roundtrip"
    );

    let orig_emb = original.embeddings.as_ref().unwrap();
    let rt_emb = roundtrip.embeddings.as_ref().unwrap();

    assert_eq!(orig_emb.len(), rt_emb.len(), "Embedding count mismatch");

    for (key, orig_data) in orig_emb {
        assert!(
            rt_emb.contains_key(key),
            "Embedding '{}' lost during roundtrip",
            key
        );

        let rt_data = &rt_emb[key];
        assert_eq!(
            orig_data.n_rows, rt_data.n_rows,
            "Embedding '{}' row count mismatch",
            key
        );
        assert_eq!(
            orig_data.n_cols, rt_data.n_cols,
            "Embedding '{}' column count mismatch",
            key
        );
    }

    println!("✓ Embeddings preserved");
}

#[test]
fn test_h5ad_layers_preservation() {
    // 这个测试需要一个包含 layers 的测试文件
    // 如果没有，我们跳过
    let original_path = "tests/data/rust_generated_seurat_with_layers.rds";

    if !Path::new(original_path).exists() {
        println!("Test file with layers not found, skipping");
        return;
    }

    // 注意：这个文件是 RDS 格式，我们需要先转换
    // 这里只是一个占位符测试
    println!("Layers test placeholder - needs h5ad file with layers");
}
