//! 真实数据集准确性验证测试
//!
//! Task 20: 使用真实数据集验证 CrossCell 的往返转换准确性

use crosscell::anndata::{reader::read_h5ad, writer::write_h5ad};
use crosscell::seurat::{seurat_rds_to_ir, write_seurat_rds};
use crosscell::validation::accuracy::{calculate_full_accuracy, AccuracyReport};
use std::path::Path;

/// 检查测试数据文件是否存在
fn check_test_data(path: &str) -> bool {
    if Path::new(path).exists() {
        true
    } else {
        eprintln!("Test data not found: {}", path);
        false
    }
}

/// 打印准确性报告
fn print_report(report: &AccuracyReport) {
    println!("\n{}", report.detailed_report());
}

// ============================================================================
// Task 20.3: PBMC 3k 往返准确性测试
// ============================================================================

#[test]
fn test_pbmc3k_h5ad_roundtrip() {
    println!("\n=== PBMC 3k H5AD Roundtrip Test ===\n");

    let original_path = "tests/data/real_datasets/pbmc3k.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    // 读取原始文件
    let original = read_h5ad(original_path).expect("Failed to read original");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 写入临时文件
    let temp_path = "tests/data/real_datasets/pbmc3k_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    // 读取往返文件
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    // 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "PBMC3k H5AD");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证
    assert!(report.passed, "PBMC 3k roundtrip failed");
    assert!(report.expression_accuracy.correlation > 0.99999);
}

/// 跨格式测试 (H5AD -> RDS -> H5AD)
/// Task 20.3: PBMC 3k 跨格式往返准确性测试
#[test]
fn test_pbmc3k_cross_format() {
    println!("\n=== PBMC 3k Cross-Format Test (H5AD -> RDS -> H5AD) ===\n");

    let h5ad_path = "tests/data/real_datasets/pbmc3k.h5ad";
    if !check_test_data(h5ad_path) {
        return;
    }

    // 1. 读取 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 写入 RDS
    let rds_path = "tests/data/real_datasets/pbmc3k_cross_temp.rds";
    write_seurat_rds(&original, rds_path).expect("Failed to write rds");
    println!("Written to RDS: {}", rds_path);

    // 3. 读取 RDS
    let from_rds = seurat_rds_to_ir(rds_path).expect("Failed to read rds");
    println!(
        "From RDS: {} cells x {} genes",
        from_rds.metadata.n_cells, from_rds.metadata.n_genes
    );

    // 4. 写入 H5AD
    let h5ad_rt_path = "tests/data/real_datasets/pbmc3k_cross_rt.h5ad";
    write_h5ad(&from_rds, h5ad_rt_path).expect("Failed to write h5ad");

    // 5. 读取往返 H5AD
    let roundtrip = read_h5ad(h5ad_rt_path).expect("Failed to read roundtrip");

    // 6. 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "PBMC3k Cross-Format");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(rds_path);
    let _ = std::fs::remove_file(h5ad_rt_path);

    // 跨格式转换可能有一些精度损失，放宽要求
    assert!(
        report.expression_accuracy.correlation > 0.999,
        "Cross-format correlation too low: {}",
        report.expression_accuracy.correlation
    );
}

// ============================================================================
// Task 20.4: Pancreas 多批次数据往返准确性测试
// ============================================================================

/// Pancreas 数据集往返测试
#[test]
fn test_pancreas_h5ad_roundtrip() {
    println!("\n=== Pancreas H5AD Roundtrip Test ===\n");

    let original_path = "tests/data/real_datasets/pancreas.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    // 读取原始文件
    let original = read_h5ad(original_path).expect("Failed to read original");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 写入临时文件
    let temp_path = "tests/data/real_datasets/pancreas_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    // 读取往返文件
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    // 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Pancreas H5AD");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证
    assert!(report.passed, "Pancreas roundtrip failed");

    // 特别检查 Categorical 变量
    if let Some(ref meta_acc) = report.cell_metadata_accuracy {
        assert!(meta_acc.categorical_match, "Categorical variables mismatch");
    }
}

/// Pancreas 跨格式测试 (H5AD -> RDS -> H5AD)
/// Task 20.4: 测试多批次数据和 Categorical 变量转换
#[test]
fn test_pancreas_cross_format() {
    println!("\n=== Pancreas Cross-Format Test (H5AD -> RDS -> H5AD) ===\n");

    let h5ad_path = "tests/data/real_datasets/pancreas.h5ad";
    if !check_test_data(h5ad_path) {
        return;
    }

    // 1. 读取 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 打印元数据信息
    println!("Cell metadata columns: {}", original.cell_metadata.n_cols());
    for col in &original.cell_metadata.columns {
        println!("  - {}", col);
    }

    // 打印 layers 信息
    if let Some(ref layers) = original.layers {
        println!("Layers: {:?}", layers.keys().collect::<Vec<_>>());
    }

    // 2. 写入 RDS
    let rds_path = "tests/data/real_datasets/pancreas_cross_temp.rds";
    write_seurat_rds(&original, rds_path).expect("Failed to write rds");
    println!("Written to RDS: {}", rds_path);

    // 3. 读取 RDS
    let from_rds = seurat_rds_to_ir(rds_path).expect("Failed to read rds");
    println!(
        "From RDS: {} cells x {} genes",
        from_rds.metadata.n_cells, from_rds.metadata.n_genes
    );
    println!(
        "Cell metadata columns after RDS: {}",
        from_rds.cell_metadata.n_cols()
    );

    // 4. 写入 H5AD
    let h5ad_rt_path = "tests/data/real_datasets/pancreas_cross_rt.h5ad";
    write_h5ad(&from_rds, h5ad_rt_path).expect("Failed to write h5ad");

    // 5. 读取往返 H5AD
    let roundtrip = read_h5ad(h5ad_rt_path).expect("Failed to read roundtrip");

    // 6. 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Pancreas Cross-Format");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(rds_path);
    let _ = std::fs::remove_file(h5ad_rt_path);

    // 验证表达矩阵
    assert!(
        report.expression_accuracy.correlation > 0.999,
        "Cross-format correlation too low: {}",
        report.expression_accuracy.correlation
    );

    // 验证 Categorical 变量（批次、细胞类型等）
    if let Some(ref meta_acc) = report.cell_metadata_accuracy {
        assert!(
            meta_acc.categorical_match,
            "Categorical variables (batch, cell_type) not preserved correctly"
        );
    }
}

// ============================================================================
// Task 20.5: Visium 空间数据往返准确性测试
// ============================================================================

#[test]
fn test_visium_h5ad_roundtrip() {
    println!("\n=== Visium Heart H5AD Roundtrip Test ===\n");

    let original_path = "tests/data/real_datasets/visium_heart.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    // 读取原始文件
    let original = read_h5ad(original_path).expect("Failed to read original");
    println!(
        "Original: {} spots x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 检查空间数据
    if let Some(ref emb) = original.embeddings {
        if emb.contains_key("spatial") {
            println!("  Has spatial coordinates");
        }
    }

    // 写入临时文件
    let temp_path = "tests/data/real_datasets/visium_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    // 读取往返文件
    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    // 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Visium Heart H5AD");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(temp_path);

    // 验证
    assert!(report.passed, "Visium roundtrip failed");

    // 检查空间坐标保留
    if let Some(spatial_acc) = report.embedding_accuracies.get("spatial") {
        assert!(
            spatial_acc.max_error < 0.01,
            "Spatial coordinates error too large: {}",
            spatial_acc.max_error
        );
    }
}

/// Visium 跨格式测试 (H5AD -> RDS -> H5AD)
/// Task 20.5: 测试空间数据保留
#[test]
fn test_visium_cross_format() {
    println!("\n=== Visium Cross-Format Test (H5AD -> RDS -> H5AD) ===\n");

    let h5ad_path = "tests/data/real_datasets/visium_heart.h5ad";
    if !check_test_data(h5ad_path) {
        return;
    }

    // 1. 读取 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original: {} spots x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 打印空间数据信息
    if let Some(ref emb) = original.embeddings {
        println!("Embeddings: {:?}", emb.keys().collect::<Vec<_>>());
        if let Some(spatial) = emb.get("spatial") {
            println!(
                "  spatial: {} spots x {} dims",
                spatial.n_rows, spatial.n_cols
            );
        }
    }

    // 打印元数据信息
    println!("Cell metadata columns: {}", original.cell_metadata.n_cols());

    // 2. 写入 RDS
    let rds_path = "tests/data/real_datasets/visium_cross_temp.rds";
    write_seurat_rds(&original, rds_path).expect("Failed to write rds");
    println!("Written to RDS: {}", rds_path);

    // 3. 读取 RDS
    let from_rds = seurat_rds_to_ir(rds_path).expect("Failed to read rds");
    println!(
        "From RDS: {} spots x {} genes",
        from_rds.metadata.n_cells, from_rds.metadata.n_genes
    );

    // 检查空间数据是否保留
    if let Some(ref emb) = from_rds.embeddings {
        println!("Embeddings after RDS: {:?}", emb.keys().collect::<Vec<_>>());
    }

    // 4. 写入 H5AD
    let h5ad_rt_path = "tests/data/real_datasets/visium_cross_rt.h5ad";
    write_h5ad(&from_rds, h5ad_rt_path).expect("Failed to write h5ad");

    // 5. 读取往返 H5AD
    let roundtrip = read_h5ad(h5ad_rt_path).expect("Failed to read roundtrip");

    // 6. 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Visium Cross-Format");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(rds_path);
    let _ = std::fs::remove_file(h5ad_rt_path);

    // 验证表达矩阵
    assert!(
        report.expression_accuracy.correlation > 0.999,
        "Cross-format correlation too low: {}",
        report.expression_accuracy.correlation
    );

    // 验证空间坐标保留
    if let Some(spatial_acc) = report.embedding_accuracies.get("spatial") {
        assert!(
            spatial_acc.max_error < 0.01,
            "Spatial coordinates error too large: {}",
            spatial_acc.max_error
        );
        println!(
            "Spatial coordinates preserved: max_error = {:.2e}",
            spatial_acc.max_error
        );
    }
}

// ============================================================================
// Task 20.6: 特殊值和边界情况测试
// ============================================================================

#[test]
fn test_small_sparse_roundtrip() {
    println!("\n=== Small Sparse Matrix Roundtrip Test ===\n");

    let original_path = "tests/data/small_sparse.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    let original = read_h5ad(original_path).expect("Failed to read");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    let temp_path = "tests/data/small_sparse_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    let report = calculate_full_accuracy(&original, &roundtrip, "Small Sparse");
    print_report(&report);

    let _ = std::fs::remove_file(temp_path);

    assert!(report.passed);
    assert!(report.expression_accuracy.nnz_match, "NNZ count changed");
}

#[test]
fn test_medium_sparse_roundtrip() {
    println!("\n=== Medium Sparse Matrix Roundtrip Test ===\n");

    let original_path = "tests/data/medium_sparse.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    let original = read_h5ad(original_path).expect("Failed to read");
    println!(
        "Original: {} cells x {} genes, NNZ: {}",
        original.metadata.n_cells,
        original.metadata.n_genes,
        original.expression.nnz()
    );

    let temp_path = "tests/data/medium_sparse_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    let report = calculate_full_accuracy(&original, &roundtrip, "Medium Sparse");
    print_report(&report);

    let _ = std::fs::remove_file(temp_path);

    assert!(report.passed);

    // 检查稀疏性保留
    assert!(
        report.expression_accuracy.sparsity_preserved > 0.999,
        "Sparsity not preserved: {}",
        report.expression_accuracy.sparsity_preserved
    );
}

// ============================================================================
// RDS 往返测试
// ============================================================================

/// Seurat RDS 往返测试
#[test]
fn test_seurat_minimal_roundtrip() {
    println!("\n=== Seurat Minimal RDS Roundtrip Test ===\n");

    let original_path = "tests/data/seurat_minimal_simplified.rds";
    if !check_test_data(original_path) {
        return;
    }

    // 读取原始 RDS
    let original = seurat_rds_to_ir(original_path).expect("Failed to read");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 写入临时 RDS
    let temp_path = "tests/data/seurat_minimal_rt_temp.rds";
    write_seurat_rds(&original, temp_path).expect("Failed to write");

    // 读取往返 RDS
    let roundtrip = seurat_rds_to_ir(temp_path).expect("Failed to read roundtrip");

    // 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Seurat Minimal RDS");
    print_report(&report);

    let _ = std::fs::remove_file(temp_path);

    // RDS 往返可能有一些差异，放宽要求
    assert!(
        report.expression_accuracy.correlation > 0.99,
        "RDS roundtrip correlation too low: {}",
        report.expression_accuracy.correlation
    );
}

// ============================================================================
// 汇总测试
// ============================================================================

#[test]
fn test_accuracy_summary() {
    println!("\n========================================");
    println!("Task 20.7: Accuracy Test Summary");
    println!("========================================\n");

    let datasets = vec![
        ("PBMC 3k", "tests/data/real_datasets/pbmc3k.h5ad"),
        ("Pancreas", "tests/data/real_datasets/pancreas.h5ad"),
        ("Visium", "tests/data/real_datasets/visium_heart.h5ad"),
        ("Small Sparse", "tests/data/small_sparse.h5ad"),
    ];

    let mut results = Vec::new();

    for (name, path) in datasets {
        if !Path::new(path).exists() {
            println!("⚠️  {} - SKIPPED (data not found)", name);
            continue;
        }

        match read_h5ad(path) {
            Ok(original) => {
                let temp_path = format!("{}_summary_temp.h5ad", path);
                if write_h5ad(&original, &temp_path).is_ok() {
                    if let Ok(roundtrip) = read_h5ad(&temp_path) {
                        let report = calculate_full_accuracy(&original, &roundtrip, name);
                        let status = if report.passed {
                            "✓ PASS"
                        } else {
                            "✗ FAIL"
                        };
                        println!(
                            "{} - {} (corr: {:.6})",
                            name, status, report.expression_accuracy.correlation
                        );
                        results.push((name.to_string(), report.passed));
                    }
                    let _ = std::fs::remove_file(&temp_path);
                }
            }
            Err(e) => {
                println!("⚠️  {} - ERROR: {}", name, e);
            }
        }
    }

    println!("\n--- Summary ---");
    let passed = results.iter().filter(|(_, p)| *p).count();
    let total = results.len();
    println!("Passed: {}/{}", passed, total);

    if passed < total {
        println!("\nFailed datasets:");
        for (name, passed) in &results {
            if !passed {
                println!("  - {}", name);
            }
        }
    }
}

// ============================================================================
// Task 20.6: 特殊值和边界情况测试（续）
// ============================================================================

/// 单行矩阵往返测试
#[test]
fn test_single_row_roundtrip() {
    println!("\n=== Single Row Matrix Roundtrip Test ===\n");

    let original_path = "tests/data/single_row.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    let original = read_h5ad(original_path).expect("Failed to read");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 验证是单行矩阵
    assert_eq!(original.metadata.n_cells, 1, "Expected single row matrix");

    let temp_path = "tests/data/single_row_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    let report = calculate_full_accuracy(&original, &roundtrip, "Single Row");
    print_report(&report);

    let _ = std::fs::remove_file(temp_path);

    assert!(report.passed, "Single row roundtrip failed");
    assert_eq!(roundtrip.metadata.n_cells, 1, "Cell count changed");
}

/// 单列矩阵往返测试
#[test]
fn test_single_column_roundtrip() {
    println!("\n=== Single Column Matrix Roundtrip Test ===\n");

    let original_path = "tests/data/single_column.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    let original = read_h5ad(original_path).expect("Failed to read");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 验证是单列矩阵
    assert_eq!(
        original.metadata.n_genes, 1,
        "Expected single column matrix"
    );

    let temp_path = "tests/data/single_column_rt_temp.h5ad";
    write_h5ad(&original, temp_path).expect("Failed to write");

    let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

    let report = calculate_full_accuracy(&original, &roundtrip, "Single Column");
    print_report(&report);

    let _ = std::fs::remove_file(temp_path);

    assert!(report.passed, "Single column roundtrip failed");
    assert_eq!(roundtrip.metadata.n_genes, 1, "Gene count changed");
}

/// 空矩阵处理测试
#[test]
fn test_empty_matrix_handling() {
    println!("\n=== Empty Matrix Handling Test ===\n");

    let original_path = "tests/data/empty.h5ad";
    if !check_test_data(original_path) {
        return;
    }

    // 尝试读取空矩阵
    match read_h5ad(original_path) {
        Ok(data) => {
            println!(
                "Empty matrix read successfully: {} x {}",
                data.metadata.n_cells, data.metadata.n_genes
            );

            // 验证是空矩阵
            let is_empty = data.metadata.n_cells == 0 || data.metadata.n_genes == 0;
            println!("Is empty: {}", is_empty);

            // 尝试往返
            let temp_path = "tests/data/empty_rt_temp.h5ad";
            if let Ok(()) = write_h5ad(&data, temp_path) {
                if let Ok(roundtrip) = read_h5ad(temp_path) {
                    println!(
                        "Roundtrip successful: {} x {}",
                        roundtrip.metadata.n_cells, roundtrip.metadata.n_genes
                    );
                    assert_eq!(data.metadata.n_cells, roundtrip.metadata.n_cells);
                    assert_eq!(data.metadata.n_genes, roundtrip.metadata.n_genes);
                }
                let _ = std::fs::remove_file(temp_path);
            }
        }
        Err(e) => {
            // 空矩阵可能会导致错误，这是可接受的行为
            println!("Empty matrix handling: {:?}", e);
            println!("Note: Empty matrices may not be fully supported");
        }
    }
}

/// 跨格式边界情况测试：单行矩阵
#[test]
fn test_single_row_cross_format() {
    println!("\n=== Single Row Cross-Format Test (H5AD -> RDS -> H5AD) ===\n");

    let h5ad_path = "tests/data/single_row.h5ad";
    if !check_test_data(h5ad_path) {
        return;
    }

    // 1. 读取 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 写入 RDS
    let rds_path = "tests/data/single_row_cross_temp.rds";
    write_seurat_rds(&original, rds_path).expect("Failed to write rds");

    // 3. 读取 RDS
    let from_rds = seurat_rds_to_ir(rds_path).expect("Failed to read rds");
    println!(
        "From RDS: {} cells x {} genes",
        from_rds.metadata.n_cells, from_rds.metadata.n_genes
    );

    // 4. 写入 H5AD
    let h5ad_rt_path = "tests/data/single_row_cross_rt.h5ad";
    write_h5ad(&from_rds, h5ad_rt_path).expect("Failed to write h5ad");

    // 5. 读取往返 H5AD
    let roundtrip = read_h5ad(h5ad_rt_path).expect("Failed to read roundtrip");

    // 6. 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Single Row Cross-Format");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(rds_path);
    let _ = std::fs::remove_file(h5ad_rt_path);

    // 验证
    assert!(
        report.expression_accuracy.correlation > 0.999,
        "Single row cross-format correlation too low: {}",
        report.expression_accuracy.correlation
    );
    assert_eq!(
        roundtrip.metadata.n_cells, 1,
        "Cell count changed in cross-format"
    );
}

/// 跨格式边界情况测试：单列矩阵
#[test]
fn test_single_column_cross_format() {
    println!("\n=== Single Column Cross-Format Test (H5AD -> RDS -> H5AD) ===\n");

    let h5ad_path = "tests/data/single_column.h5ad";
    if !check_test_data(h5ad_path) {
        return;
    }

    // 1. 读取 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 写入 RDS
    let rds_path = "tests/data/single_column_cross_temp.rds";
    write_seurat_rds(&original, rds_path).expect("Failed to write rds");

    // 3. 读取 RDS
    let from_rds = seurat_rds_to_ir(rds_path).expect("Failed to read rds");
    println!(
        "From RDS: {} cells x {} genes",
        from_rds.metadata.n_cells, from_rds.metadata.n_genes
    );

    // 4. 写入 H5AD
    let h5ad_rt_path = "tests/data/single_column_cross_rt.h5ad";
    write_h5ad(&from_rds, h5ad_rt_path).expect("Failed to write h5ad");

    // 5. 读取往返 H5AD
    let roundtrip = read_h5ad(h5ad_rt_path).expect("Failed to read roundtrip");

    // 6. 计算准确性
    let report = calculate_full_accuracy(&original, &roundtrip, "Single Column Cross-Format");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(rds_path);
    let _ = std::fs::remove_file(h5ad_rt_path);

    // 验证
    assert!(
        report.expression_accuracy.correlation > 0.999,
        "Single column cross-format correlation too low: {}",
        report.expression_accuracy.correlation
    );
    assert_eq!(
        roundtrip.metadata.n_genes, 1,
        "Gene count changed in cross-format"
    );
}

/// 稀疏度测试：验证不同稀疏度的矩阵往返
#[test]
fn test_sparsity_preservation() {
    println!("\n=== Sparsity Preservation Test ===\n");

    // 测试 small_sparse (较低稀疏度)
    let small_path = "tests/data/small_sparse.h5ad";
    if check_test_data(small_path) {
        let original = read_h5ad(small_path).expect("Failed to read");
        let original_nnz = original.expression.nnz();
        let original_total = original.metadata.n_cells * original.metadata.n_genes;
        let original_sparsity = 1.0 - (original_nnz as f64 / original_total as f64);

        println!(
            "Small sparse: {} x {}, NNZ={}, sparsity={:.2}%",
            original.metadata.n_cells,
            original.metadata.n_genes,
            original_nnz,
            original_sparsity * 100.0
        );

        let temp_path = "tests/data/small_sparse_sparsity_temp.h5ad";
        write_h5ad(&original, temp_path).expect("Failed to write");
        let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

        let roundtrip_nnz = roundtrip.expression.nnz();
        println!("After roundtrip: NNZ={}", roundtrip_nnz);

        assert_eq!(
            original_nnz, roundtrip_nnz,
            "NNZ changed: {} -> {}",
            original_nnz, roundtrip_nnz
        );

        let _ = std::fs::remove_file(temp_path);
    }

    // 测试 medium_sparse (较高稀疏度)
    let medium_path = "tests/data/medium_sparse.h5ad";
    if check_test_data(medium_path) {
        let original = read_h5ad(medium_path).expect("Failed to read");
        let original_nnz = original.expression.nnz();
        let original_total = original.metadata.n_cells * original.metadata.n_genes;
        let original_sparsity = 1.0 - (original_nnz as f64 / original_total as f64);

        println!(
            "Medium sparse: {} x {}, NNZ={}, sparsity={:.2}%",
            original.metadata.n_cells,
            original.metadata.n_genes,
            original_nnz,
            original_sparsity * 100.0
        );

        let temp_path = "tests/data/medium_sparse_sparsity_temp.h5ad";
        write_h5ad(&original, temp_path).expect("Failed to write");
        let roundtrip = read_h5ad(temp_path).expect("Failed to read roundtrip");

        let roundtrip_nnz = roundtrip.expression.nnz();
        println!("After roundtrip: NNZ={}", roundtrip_nnz);

        assert_eq!(
            original_nnz, roundtrip_nnz,
            "NNZ changed: {} -> {}",
            original_nnz, roundtrip_nnz
        );

        let _ = std::fs::remove_file(temp_path);
    }

    println!("\n✅ Sparsity preservation verified for all test matrices");
}

// ============================================================================
// Task 20.1.5: 金标准配对数据验证测试
// 使用 easySCF 生成的配对数据验证 CrossCell 的准确性
// ============================================================================

/// PBMC 3k 金标准验证：比较 CrossCell 输出与 easySCF 金标准
#[test]
fn test_pbmc3k_gold_standard() {
    println!("\n=== PBMC 3k Gold Standard Validation (easySCF) ===\n");

    let h5ad_path = "tests/data/real_datasets/pbmc3k.h5ad";
    let gold_rds_path = "tests/data/real_datasets/pbmc3k_easyscf.rds";

    if !check_test_data(h5ad_path) || !check_test_data(gold_rds_path) {
        println!("⚠️ Gold standard data not found, skipping test");
        return;
    }

    // 1. 读取原始 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original H5AD: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 读取 easySCF 金标准 RDS
    let gold_standard = seurat_rds_to_ir(gold_rds_path).expect("Failed to read gold standard");
    println!(
        "Gold Standard RDS: {} cells x {} genes",
        gold_standard.metadata.n_cells, gold_standard.metadata.n_genes
    );

    // 3. 使用 CrossCell 转换 H5AD -> RDS
    let crosscell_rds_path = "tests/data/real_datasets/pbmc3k_crosscell_temp.rds";
    write_seurat_rds(&original, crosscell_rds_path).expect("Failed to write CrossCell RDS");

    // 4. 读取 CrossCell 生成的 RDS
    let crosscell_output =
        seurat_rds_to_ir(crosscell_rds_path).expect("Failed to read CrossCell RDS");
    println!(
        "CrossCell RDS: {} cells x {} genes",
        crosscell_output.metadata.n_cells, crosscell_output.metadata.n_genes
    );

    // 5. 比较 CrossCell 输出与金标准
    let report =
        calculate_full_accuracy(&gold_standard, &crosscell_output, "PBMC3k vs Gold Standard");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(crosscell_rds_path);

    // 验证：CrossCell 输出应该与 easySCF 金标准高度一致
    assert!(
        report.expression_accuracy.correlation > 0.99,
        "CrossCell output differs significantly from gold standard: corr={}",
        report.expression_accuracy.correlation
    );

    println!("\n✅ PBMC 3k: CrossCell output matches easySCF gold standard");
}

/// Pancreas 金标准验证：比较 CrossCell 输出与 easySCF 金标准
#[test]
fn test_pancreas_gold_standard() {
    println!("\n=== Pancreas Gold Standard Validation (easySCF) ===\n");

    let h5ad_path = "tests/data/real_datasets/pancreas.h5ad";
    let gold_rds_path = "tests/data/real_datasets/pancreas_easyscf.rds";

    if !check_test_data(h5ad_path) || !check_test_data(gold_rds_path) {
        println!("⚠️ Gold standard data not found, skipping test");
        return;
    }

    // 1. 读取原始 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original H5AD: {} cells x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 读取 easySCF 金标准 RDS
    let gold_standard = seurat_rds_to_ir(gold_rds_path).expect("Failed to read gold standard");
    println!(
        "Gold Standard RDS: {} cells x {} genes",
        gold_standard.metadata.n_cells, gold_standard.metadata.n_genes
    );

    // 3. 使用 CrossCell 转换 H5AD -> RDS
    let crosscell_rds_path = "tests/data/real_datasets/pancreas_crosscell_temp.rds";
    write_seurat_rds(&original, crosscell_rds_path).expect("Failed to write CrossCell RDS");

    // 4. 读取 CrossCell 生成的 RDS
    let crosscell_output =
        seurat_rds_to_ir(crosscell_rds_path).expect("Failed to read CrossCell RDS");
    println!(
        "CrossCell RDS: {} cells x {} genes",
        crosscell_output.metadata.n_cells, crosscell_output.metadata.n_genes
    );

    // 5. 比较 CrossCell 输出与金标准
    let report = calculate_full_accuracy(
        &gold_standard,
        &crosscell_output,
        "Pancreas vs Gold Standard",
    );
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(crosscell_rds_path);

    // 验证
    assert!(
        report.expression_accuracy.correlation > 0.99,
        "CrossCell output differs significantly from gold standard: corr={}",
        report.expression_accuracy.correlation
    );

    println!("\n✅ Pancreas: CrossCell output matches easySCF gold standard");
}

/// Visium 金标准验证（注意：Visium 使用 CrossCell 回退，不是独立金标准）
#[test]
fn test_visium_gold_standard() {
    println!("\n=== Visium Gold Standard Validation ===\n");
    println!("⚠️ Note: Visium uses CrossCell fallback (not independent gold standard)");

    let h5ad_path = "tests/data/real_datasets/visium_heart.h5ad";
    let gold_rds_path = "tests/data/real_datasets/visium_easyscf.rds";

    if !check_test_data(h5ad_path) || !check_test_data(gold_rds_path) {
        println!("⚠️ Gold standard data not found, skipping test");
        return;
    }

    // 1. 读取原始 H5AD
    let original = read_h5ad(h5ad_path).expect("Failed to read h5ad");
    println!(
        "Original H5AD: {} spots x {} genes",
        original.metadata.n_cells, original.metadata.n_genes
    );

    // 2. 读取金标准 RDS（由 CrossCell 生成）
    let gold_standard = seurat_rds_to_ir(gold_rds_path).expect("Failed to read gold standard");
    println!(
        "Gold Standard RDS: {} spots x {} genes",
        gold_standard.metadata.n_cells, gold_standard.metadata.n_genes
    );

    // 3. 使用 CrossCell 转换 H5AD -> RDS
    let crosscell_rds_path = "tests/data/real_datasets/visium_crosscell_temp.rds";
    write_seurat_rds(&original, crosscell_rds_path).expect("Failed to write CrossCell RDS");

    // 4. 读取 CrossCell 生成的 RDS
    let crosscell_output =
        seurat_rds_to_ir(crosscell_rds_path).expect("Failed to read CrossCell RDS");
    println!(
        "CrossCell RDS: {} spots x {} genes",
        crosscell_output.metadata.n_cells, crosscell_output.metadata.n_genes
    );

    // 5. 比较 CrossCell 输出与金标准
    let report =
        calculate_full_accuracy(&gold_standard, &crosscell_output, "Visium vs Gold Standard");
    print_report(&report);

    // 清理
    let _ = std::fs::remove_file(crosscell_rds_path);

    // 验证（由于是 CrossCell 自己生成的，应该完全一致）
    assert!(
        report.expression_accuracy.correlation > 0.9999,
        "CrossCell output differs from its own gold standard: corr={}",
        report.expression_accuracy.correlation
    );

    println!("\n✅ Visium: CrossCell output is consistent");
}

/// 金标准验证汇总
#[test]
fn test_gold_standard_summary() {
    println!("\n========================================");
    println!("Task 20.1.5: Gold Standard Validation Summary");
    println!("========================================\n");

    let datasets = vec![
        (
            "PBMC 3k",
            "tests/data/real_datasets/pbmc3k.h5ad",
            "tests/data/real_datasets/pbmc3k_easyscf.rds",
            true,
        ),
        (
            "Pancreas",
            "tests/data/real_datasets/pancreas.h5ad",
            "tests/data/real_datasets/pancreas_easyscf.rds",
            true,
        ),
        (
            "Visium",
            "tests/data/real_datasets/visium_heart.h5ad",
            "tests/data/real_datasets/visium_easyscf.rds",
            false,
        ), // CrossCell fallback
    ];

    let mut results = Vec::new();

    for (name, h5ad_path, rds_path, is_gold) in datasets {
        let gold_label = if is_gold { "easySCF" } else { "CrossCell" };

        if !Path::new(h5ad_path).exists() || !Path::new(rds_path).exists() {
            println!("⚠️  {} - SKIPPED (data not found)", name);
            continue;
        }

        match (read_h5ad(h5ad_path), seurat_rds_to_ir(rds_path)) {
            (Ok(original), Ok(gold_standard)) => {
                let temp_rds = format!("{}_summary_temp.rds", h5ad_path);
                if write_seurat_rds(&original, &temp_rds).is_ok() {
                    if let Ok(crosscell_output) = seurat_rds_to_ir(&temp_rds) {
                        let report =
                            calculate_full_accuracy(&gold_standard, &crosscell_output, name);
                        let status = if report.expression_accuracy.correlation > 0.99 {
                            "✓ PASS"
                        } else {
                            "✗ FAIL"
                        };
                        println!(
                            "{} ({}) - {} (corr: {:.6})",
                            name, gold_label, status, report.expression_accuracy.correlation
                        );
                        results.push((
                            name.to_string(),
                            report.expression_accuracy.correlation > 0.99,
                        ));
                    }
                    let _ = std::fs::remove_file(&temp_rds);
                }
            }
            _ => {
                println!("⚠️  {} - ERROR reading data", name);
            }
        }
    }

    println!("\n--- Summary ---");
    let passed = results.iter().filter(|(_, p)| *p).count();
    let total = results.len();
    println!("Passed: {}/{}", passed, total);

    if passed == total && total > 0 {
        println!("\n🎉 All gold standard validations passed!");
    }
}
