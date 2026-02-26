//! Fidelity Probe ??对比外部工具 (anndataR, zellkonverter, convert2anndata) 转换??H5AD
//!
//! 前置条件：需要先??benchmark 容器中运??R 脚本生成转换文件??
//!   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
//!     Rscript scripts/probe_external_tools.R \
//!     data/generated/seurat_v4_pbmc3k_processed.rds \
//!     /benchmark/results/probe_output
//!
//! 然后??dev 容器中运行本测试??
//!   docker-compose run --rm dev cargo test --test test_fidelity_probe_external -- --nocapture

use crosscell::anndata::read_h5ad;
use crosscell::ir::SingleCellData;
use crosscell::seurat::read_seurat_direct;
use crosscell::validation::calculate_full_accuracy;
use std::path::Path;

/// 对比原始 RDS IR 和某工具转换??H5AD IR
fn probe_tool(tool_name: &str, rds_path: &str, h5ad_path: &str) {
    println!("\n{}", "=".repeat(70));
    println!("TOOL: {}", tool_name);
    println!("{}", "=".repeat(70));

    if !Path::new(h5ad_path).exists() {
        println!("  SKIP: {} not found", h5ad_path);
        return;
    }

    // Load original RDS as IR
    let original = read_seurat_direct(rds_path, false).expect("Failed to load original RDS");
    let original_ir = &original.data;
    let (orig_cells, orig_genes) = original_ir.expression.shape();

    // Load tool-converted H5AD as IR
    let converted_ir = match read_h5ad(h5ad_path) {
        Ok(ir) => ir,
        Err(e) => {
            println!("  ??Failed to read H5AD: {}", e);
            println!("  (This itself is a fidelity finding ??the tool produced an H5AD that CrossCell cannot parse)");
            return;
        }
    };
    let (conv_cells, conv_genes) = converted_ir.expression.shape();

    println!(
        "  Original:  {} cells × {} genes, NNZ={}",
        orig_cells,
        orig_genes,
        original_ir.expression.nnz()
    );
    println!(
        "  Converted: {} cells × {} genes, NNZ={}",
        conv_cells,
        conv_genes,
        converted_ir.expression.nnz()
    );

    // Dimension check
    if orig_cells != conv_cells || orig_genes != conv_genes {
        println!("  ??DIMENSION MISMATCH ??cannot compute full accuracy");
        println!("    Cell diff: {} ??{}", orig_cells, conv_cells);
        println!("    Gene diff: {} ??{}", orig_genes, conv_genes);

        // Still try to report metadata
        println!(
            "  Original cell meta cols ({}): {:?}",
            original_ir.cell_metadata.n_cols(),
            original_ir.cell_metadata.columns
        );
        println!(
            "  Converted cell meta cols ({}): {:?}",
            converted_ir.cell_metadata.n_cols(),
            converted_ir.cell_metadata.columns
        );
        println!(
            "  Original gene meta cols ({}): {:?}",
            original_ir.gene_metadata.n_cols(),
            original_ir.gene_metadata.columns
        );
        println!(
            "  Converted gene meta cols ({}): {:?}",
            converted_ir.gene_metadata.n_cols(),
            converted_ir.gene_metadata.columns
        );

        if let Some(ref emb) = original_ir.embeddings {
            let names: Vec<_> = emb.keys().collect();
            println!("  Original embeddings: {:?}", names);
        }
        if let Some(ref emb) = converted_ir.embeddings {
            let names: Vec<_> = emb.keys().collect();
            println!("  Converted embeddings: {:?}", names);
        }
        return;
    }

    // Full accuracy report
    let report = calculate_full_accuracy(original_ir, &converted_ir, tool_name);

    println!("\n  --- Expression ---");
    println!(
        "    Correlation:     {:.8}",
        report.expression_accuracy.correlation
    );
    println!(
        "    Max abs error:   {:.2e}",
        report.expression_accuracy.max_absolute_error
    );
    println!(
        "    Mean abs error:  {:.2e}",
        report.expression_accuracy.mean_absolute_error
    );
    println!(
        "    Sparsity:        {:.4}%",
        report.expression_accuracy.sparsity_preserved * 100.0
    );
    println!(
        "    NNZ match:       {} ({} vs {})",
        report.expression_accuracy.nnz_match,
        report.expression_accuracy.original_nnz,
        report.expression_accuracy.converted_nnz
    );

    if let Some(ref meta) = report.cell_metadata_accuracy {
        println!("\n  --- Cell Metadata ---");
        println!("    Cols match:      {}", meta.column_names_match);
        println!("    Types match:     {}", meta.dtypes_match);
        println!("    Categorical:     {}", meta.categorical_match);
        println!("    Numeric err:     {:.2e}", meta.numeric_max_error);
        if !meta.missing_columns.is_empty() {
            println!("    Missing cols:    {:?}", meta.missing_columns);
        }
        if !meta.extra_columns.is_empty() {
            println!("    Extra cols:      {:?}", meta.extra_columns);
        }
    }

    if let Some(ref meta) = report.gene_metadata_accuracy {
        println!("\n  --- Gene Metadata ---");
        println!("    Cols match:      {}", meta.column_names_match);
        println!("    Types match:     {}", meta.dtypes_match);
        if !meta.missing_columns.is_empty() {
            println!("    Missing cols:    {:?}", meta.missing_columns);
        }
        if !meta.extra_columns.is_empty() {
            println!("    Extra cols:      {:?}", meta.extra_columns);
        }
    }

    if !report.embedding_accuracies.is_empty() {
        println!("\n  --- Embeddings ---");
        for (name, acc) in &report.embedding_accuracies {
            println!(
                "    {}: corr={:.8}, max_err={:.2e}, mean_err={:.2e}",
                name, acc.correlation, acc.max_error, acc.mean_error
            );
        }
    }

    if !report.layer_accuracies.is_empty() {
        println!("\n  --- Layers ---");
        for (name, acc) in &report.layer_accuracies {
            println!(
                "    {}: corr={:.6}, max_err={:.2e}",
                name, acc.correlation, acc.max_absolute_error
            );
        }
    }

    println!("\n  Overall passed: {}", report.passed);
}

#[test]
fn test_fidelity_probe_external_tools() {
    let rds_path = "data/generated/seurat_v4_pbmc3k_processed.rds";
    if !Path::new(rds_path).exists() {
        eprintln!("SKIP: {} not found", rds_path);
        return;
    }

    let probe_dir = "benchmark/results/probe_output";
    if !Path::new(probe_dir).exists() {
        eprintln!("SKIP: {} not found ??run R probe script first", probe_dir);
        return;
    }

    probe_tool(
        "anndataR",
        rds_path,
        &format!("{}/anndatar_output.h5ad", probe_dir),
    );

    probe_tool(
        "zellkonverter",
        rds_path,
        &format!("{}/zellkonverter_output.h5ad", probe_dir),
    );

    probe_tool(
        "convert2anndata",
        rds_path,
        &format!("{}/convert2anndata_output.h5ad", probe_dir),
    );
}
