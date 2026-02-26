//! Fidelity Probe ??用真实数据快速验证跨格式往返中的信息损??
//!
//! 方向 1: RDS ??IR ??H5AD ??IR (Seurat ??AnnData 往??
//! 方向 2: H5AD ??IR ??RDS ??IR (AnnData ??Seurat 往??
//!
//! 目的：在投入完整审计框架之前，证明信息损失确实存在且可量化??

use crosscell::anndata::{read_h5ad, write_h5ad};
use crosscell::seurat::{seurat_rds_to_ir, write_seurat_rds, read_seurat_direct};
use crosscell::validation::calculate_full_accuracy;
use crosscell::ir::SingleCellData;
use std::path::Path;

/// 打印 IR 的组件概??
fn print_ir_summary(label: &str, ir: &SingleCellData) {
    let (n_cells, n_genes) = ir.expression.shape();
    println!("\n=== {} ===", label);
    println!("  Cells: {}, Genes: {}", n_cells, n_genes);
    println!("  Expression NNZ: {}", ir.expression.nnz());
    println!("  Cell metadata columns ({}): {:?}", ir.cell_metadata.n_cols(), ir.cell_metadata.columns);
    println!("  Gene metadata columns ({}): {:?}", ir.gene_metadata.n_cols(), ir.gene_metadata.columns);

    if let Some(ref emb) = ir.embeddings {
        let names: Vec<_> = emb.keys().collect();
        println!("  Embeddings ({}): {:?}", emb.len(), names);
        for (name, e) in emb {
            println!("    {}: {} × {}", name, e.n_rows, e.n_cols);
        }
    } else {
        println!("  Embeddings: None");
    }

    if let Some(ref layers) = ir.layers {
        let names: Vec<_> = layers.keys().collect();
        println!("  Layers ({}): {:?}", layers.len(), names);
    } else {
        println!("  Layers: None");
    }

    if let Some(ref pw) = ir.cell_pairwise {
        let names: Vec<_> = pw.keys().collect();
        println!("  Cell pairwise ({}): {:?}", pw.len(), names);
    } else {
        println!("  Cell pairwise: None");
    }

    if let Some(ref gl) = ir.gene_loadings {
        let names: Vec<_> = gl.keys().collect();
        println!("  Gene loadings/varm ({}): {:?}", gl.len(), names);
    } else {
        println!("  Gene loadings/varm: None");
    }
}

/// 打印两个 IR 之间的逐组件差??
fn print_component_diff(original: &SingleCellData, roundtrip: &SingleCellData) {
    println!("\n--- Component Diff ---");

    // Expression matrix
    let orig_nnz = original.expression.nnz();
    let rt_nnz = roundtrip.expression.nnz();
    println!("  Expression NNZ: {} ??{} (delta: {})", orig_nnz, rt_nnz, rt_nnz as i64 - orig_nnz as i64);

    // Cell metadata columns
    let orig_cols: std::collections::HashSet<_> = original.cell_metadata.columns.iter().collect();
    let rt_cols: std::collections::HashSet<_> = roundtrip.cell_metadata.columns.iter().collect();
    let lost_cols: Vec<_> = orig_cols.difference(&rt_cols).collect();
    let gained_cols: Vec<_> = rt_cols.difference(&orig_cols).collect();
    if !lost_cols.is_empty() {
        println!("  Cell metadata LOST columns: {:?}", lost_cols);
    }
    if !gained_cols.is_empty() {
        println!("  Cell metadata GAINED columns: {:?}", gained_cols);
    }

    // Gene metadata columns
    let orig_gcols: std::collections::HashSet<_> = original.gene_metadata.columns.iter().collect();
    let rt_gcols: std::collections::HashSet<_> = roundtrip.gene_metadata.columns.iter().collect();
    let lost_gcols: Vec<_> = orig_gcols.difference(&rt_gcols).collect();
    let gained_gcols: Vec<_> = rt_gcols.difference(&orig_gcols).collect();
    if !lost_gcols.is_empty() {
        println!("  Gene metadata LOST columns: {:?}", lost_gcols);
    }
    if !gained_gcols.is_empty() {
        println!("  Gene metadata GAINED columns: {:?}", gained_gcols);
    }

    // Type changes in shared columns
    print_type_changes("Cell metadata", &original.cell_metadata, &roundtrip.cell_metadata);
    print_type_changes("Gene metadata", &original.gene_metadata, &roundtrip.gene_metadata);

    // Embeddings
    let orig_emb_names: std::collections::HashSet<_> = original.embeddings.as_ref()
        .map(|e| e.keys().collect())
        .unwrap_or_default();
    let rt_emb_names: std::collections::HashSet<_> = roundtrip.embeddings.as_ref()
        .map(|e| e.keys().collect())
        .unwrap_or_default();
    let lost_emb: Vec<_> = orig_emb_names.difference(&rt_emb_names).collect();
    let gained_emb: Vec<_> = rt_emb_names.difference(&orig_emb_names).collect();
    if !lost_emb.is_empty() {
        println!("  Embeddings LOST: {:?}", lost_emb);
    }
    if !gained_emb.is_empty() {
        println!("  Embeddings GAINED: {:?}", gained_emb);
    }

    // Layers
    let orig_layer_names: std::collections::HashSet<_> = original.layers.as_ref()
        .map(|l| l.keys().collect())
        .unwrap_or_default();
    let rt_layer_names: std::collections::HashSet<_> = roundtrip.layers.as_ref()
        .map(|l| l.keys().collect())
        .unwrap_or_default();
    let lost_layers: Vec<_> = orig_layer_names.difference(&rt_layer_names).collect();
    if !lost_layers.is_empty() {
        println!("  Layers LOST: {:?}", lost_layers);
    }
}

/// 打印共享列的类型变化
fn print_type_changes(label: &str, original: &crosscell::ir::DataFrame, roundtrip: &crosscell::ir::DataFrame) {
    let orig_cols: std::collections::HashSet<_> = original.columns.iter().cloned().collect();
    let rt_cols: std::collections::HashSet<_> = roundtrip.columns.iter().cloned().collect();
    let shared: Vec<_> = orig_cols.intersection(&rt_cols).collect();

    let mut type_changes = Vec::new();
    for col in &shared {
        if let (Some(orig_arr), Some(rt_arr)) = (original.column(col), roundtrip.column(col)) {
            let orig_type = orig_arr.data_type().clone();
            let rt_type = rt_arr.data_type().clone();
            if orig_type != rt_type {
                type_changes.push(format!("    {} : {:?} ??{:?}", col, orig_type, rt_type));
            }
        }
    }
    if !type_changes.is_empty() {
        println!("  {} type changes:", label);
        for change in &type_changes {
            println!("{}", change);
        }
    }
}

// ============================================================================
// Direction 1: RDS ??IR ??H5AD ??IR
// ============================================================================

#[test]
fn test_fidelity_probe_rds_to_h5ad_roundtrip() {
    let rds_path = "data/generated/seurat_v4_pbmc3k_processed.rds";
    if !Path::new(rds_path).exists() {
        eprintln!("SKIP: {} not found", rds_path);
        return;
    }

    println!("\n{}", "=".repeat(70));
    println!("FIDELITY PROBE: RDS ??IR ??H5AD ??IR");
    println!("Dataset: {}", rds_path);
    println!("{}", "=".repeat(70));

    // Step 1: RDS ??IR (original) using direct reader for native S4 Seurat objects
    let direct_result = read_seurat_direct(rds_path, false)
        .expect("Failed to load RDS as IR");
    let original_ir = direct_result.data;
    println!("  Seurat version: {:?}", direct_result.version);
    if direct_result.skipped.has_skipped() {
        println!("  Skipped components: {:?}", direct_result.skipped);
    }
    print_ir_summary("Original IR (from RDS)", &original_ir);

    // Step 2: IR ??H5AD (write to temp file)
    let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
    let h5ad_path = temp_dir.path().join("probe_roundtrip.h5ad");
    write_h5ad(&original_ir, &h5ad_path)
        .expect("Failed to write H5AD");
    println!("\n  Wrote H5AD to: {:?}", h5ad_path);

    // Step 3: H5AD ??IR (read back)
    let roundtrip_ir = read_h5ad(&h5ad_path)
        .expect("Failed to read H5AD back as IR");
    print_ir_summary("Roundtrip IR (from H5AD)", &roundtrip_ir);

    // Step 4: Compare
    print_component_diff(&original_ir, &roundtrip_ir);

    // Step 5: calculate_full_accuracy
    let report = calculate_full_accuracy(&original_ir, &roundtrip_ir, "rds_to_h5ad_probe");
    println!("\n{}", report.detailed_report());

    // Print key fidelity numbers
    println!("\n--- Key Fidelity Numbers ---");
    println!("  Expression correlation:    {:.8}", report.expression_accuracy.correlation);
    println!("  Expression max abs error:  {:.2e}", report.expression_accuracy.max_absolute_error);
    println!("  Expression mean abs error: {:.2e}", report.expression_accuracy.mean_absolute_error);
    println!("  Sparsity preserved:        {:.4}%", report.expression_accuracy.sparsity_preserved * 100.0);
    println!("  NNZ match:                 {} ({} vs {})",
        report.expression_accuracy.nnz_match,
        report.expression_accuracy.original_nnz,
        report.expression_accuracy.converted_nnz);

    if let Some(ref meta) = report.cell_metadata_accuracy {
        println!("  Cell metadata cols match:  {}", meta.column_names_match);
        println!("  Cell metadata types match: {}", meta.dtypes_match);
        println!("  Cell metadata categorical: {}", meta.categorical_match);
        if !meta.missing_columns.is_empty() {
            println!("  Cell metadata missing:     {:?}", meta.missing_columns);
        }
    }

    if let Some(ref meta) = report.gene_metadata_accuracy {
        println!("  Gene metadata cols match:  {}", meta.column_names_match);
        println!("  Gene metadata types match: {}", meta.dtypes_match);
        if !meta.missing_columns.is_empty() {
            println!("  Gene metadata missing:     {:?}", meta.missing_columns);
        }
    }

    for (name, acc) in &report.embedding_accuracies {
        println!("  Embedding '{}': corr={:.8}, max_err={:.2e}, mean_err={:.2e}",
            name, acc.correlation, acc.max_error, acc.mean_error);
    }

    // We do NOT assert pass/fail here ??this is a probe to discover what losses exist.
    // The test always passes; the value is in the printed output.
    println!("\n  Overall passed (accuracy check): {}", report.passed);
}

// ============================================================================
// Direction 2: H5AD ??IR ??RDS ??IR
// ============================================================================

#[test]
fn test_fidelity_probe_h5ad_to_rds_roundtrip() {
    let h5ad_path = "data/generated/scanpy_pbmc3k.h5ad";
    if !Path::new(h5ad_path).exists() {
        eprintln!("SKIP: {} not found", h5ad_path);
        return;
    }

    println!("\n{}", "=".repeat(70));
    println!("FIDELITY PROBE: H5AD ??IR ??RDS ??IR");
    println!("Dataset: {}", h5ad_path);
    println!("{}", "=".repeat(70));

    // Step 1: H5AD ??IR (original)
    let original_ir = read_h5ad(h5ad_path)
        .expect("Failed to load H5AD as IR");
    print_ir_summary("Original IR (from H5AD)", &original_ir);

    // Step 2: IR ??RDS (write to temp file)
    let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
    let rds_path = temp_dir.path().join("probe_roundtrip.rds");
    let rds_path_str = rds_path.to_str().expect("Invalid path");
    write_seurat_rds(&original_ir, rds_path_str)
        .expect("Failed to write RDS");
    println!("\n  Wrote RDS to: {:?}", rds_path);

    // Step 3: RDS ??IR (read back)
    let roundtrip_ir = seurat_rds_to_ir(rds_path_str)
        .expect("Failed to read RDS back as IR");
    print_ir_summary("Roundtrip IR (from RDS)", &roundtrip_ir);

    // Step 4: Compare
    print_component_diff(&original_ir, &roundtrip_ir);

    // Step 5: calculate_full_accuracy
    let report = calculate_full_accuracy(&original_ir, &roundtrip_ir, "h5ad_to_rds_probe");
    println!("\n{}", report.detailed_report());

    // Print key fidelity numbers
    println!("\n--- Key Fidelity Numbers ---");
    println!("  Expression correlation:    {:.8}", report.expression_accuracy.correlation);
    println!("  Expression max abs error:  {:.2e}", report.expression_accuracy.max_absolute_error);
    println!("  Expression mean abs error: {:.2e}", report.expression_accuracy.mean_absolute_error);
    println!("  Sparsity preserved:        {:.4}%", report.expression_accuracy.sparsity_preserved * 100.0);
    println!("  NNZ match:                 {} ({} vs {})",
        report.expression_accuracy.nnz_match,
        report.expression_accuracy.original_nnz,
        report.expression_accuracy.converted_nnz);

    if let Some(ref meta) = report.cell_metadata_accuracy {
        println!("  Cell metadata cols match:  {}", meta.column_names_match);
        println!("  Cell metadata types match: {}", meta.dtypes_match);
        println!("  Cell metadata categorical: {}", meta.categorical_match);
        if !meta.missing_columns.is_empty() {
            println!("  Cell metadata missing:     {:?}", meta.missing_columns);
        }
    }

    if let Some(ref meta) = report.gene_metadata_accuracy {
        println!("  Gene metadata cols match:  {}", meta.column_names_match);
        println!("  Gene metadata types match: {}", meta.dtypes_match);
        if !meta.missing_columns.is_empty() {
            println!("  Gene metadata missing:     {:?}", meta.missing_columns);
        }
    }

    for (name, acc) in &report.embedding_accuracies {
        println!("  Embedding '{}': corr={:.8}, max_err={:.2e}, mean_err={:.2e}",
            name, acc.correlation, acc.max_error, acc.mean_error);
    }

    println!("\n  Overall passed (accuracy check): {}", report.passed);
}
