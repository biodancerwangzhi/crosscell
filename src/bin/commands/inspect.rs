//! Inspect command implementation - check file information and compatibility

use anyhow::{Context, Result};
use crosscell::{anndata, seurat};
use crosscell::diagnostics::detect_issues;
use crosscell::ir::SingleCellData;
use log::info;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

pub fn run(input: PathBuf, detailed: bool, output: Option<PathBuf>) -> Result<()> {
    info!("Inspecting file: {}", input.display());

    // Validate input file exists
    if !input.exists() {
        anyhow::bail!("Input file does not exist: {}", input.display());
    }

    // Detect file format
    let format = detect_format(&input)?;
    info!("File format: {}", format);

    // Get file size
    let file_size = fs::metadata(&input)?.len();
    let size_mb = file_size as f64 / (1024.0 * 1024.0);

    // Read and analyze file
    let data = match format.as_str() {
        "AnnData (.h5ad)" => {
            let input_str = input.to_str().context("Invalid input path")?;
            anndata::reader::read_h5ad(input_str)
                .context("Failed to read H5AD file")?
        }
        "Seurat (.rds)" => {
            let input_str = input.to_str().context("Invalid input path")?;
            // Use read_seurat_direct for raw Seurat RDS files
            let result = seurat::read_seurat_direct(std::path::Path::new(input_str), false)
                .map_err(|e| anyhow::anyhow!("Failed to read Seurat RDS file: {}", e))?;
            result.data
        }
        _ => anyhow::bail!("Unsupported file format: {}", format),
    };

    // Extract basic info
    let n_cells = data.metadata.n_cells;
    let n_genes = data.metadata.n_genes;
    let has_embeddings = data.embeddings.as_ref().map_or(false, |e| !e.is_empty());
    let has_layers = data.layers.as_ref().map_or(false, |l| !l.is_empty());
    let sparsity = calculate_sparsity(&data);

    // Run diagnostics
    let diagnostic_report = detect_issues(&data);

    // Display basic information
    println!("📊 File: {}", input.display());
    println!("Format: {} | Size: {:.2} MB", format, size_mb);
    println!("Cells: {} | Genes: {} | Sparsity: {:.1}%", n_cells, n_genes, sparsity);
    println!();

    // Display data components
    println!("✅ Data Components");
    println!("  ├─ Expression matrix: {}", if sparsity > 50.0 { "Sparse" } else { "Dense" });
    println!("  ├─ Cell metadata (obs): {} columns", data.cell_metadata.n_cols());
    println!("  ├─ Gene metadata (var): {} columns", data.gene_metadata.n_cols());
    
    if has_embeddings {
        let emb_count = data.embeddings.as_ref().map_or(0, |e| e.len());
        let emb_names: Vec<_> = data.embeddings.as_ref()
            .map_or(vec![], |e| e.keys().cloned().collect());
        println!("  ├─ Embeddings (obsm): {} ({:?})", emb_count, emb_names);
    }
    
    if has_layers {
        let layer_count = data.layers.as_ref().map_or(0, |l| l.len());
        let layer_names: Vec<_> = data.layers.as_ref()
            .map_or(vec![], |l| l.keys().cloned().collect());
        println!("  ├─ Layers: {} ({:?})", layer_count, layer_names);
    }
    
    if data.spatial.is_some() {
        println!("  └─ Spatial data: Yes");
    } else {
        println!("  └─ Spatial data: No");
    }

    // Display compatibility score and issues
    println!();
    println!("🎯 Conversion Compatibility: {}/100", diagnostic_report.compatibility_score);
    
    if diagnostic_report.has_issues() {
        println!();
        let auto_fix_text = if diagnostic_report.auto_fixable_count > 0 {
            format!(" ({} auto-fixable)", diagnostic_report.auto_fixable_count)
        } else {
            String::new()
        };
        
        println!("⚠️  Detected {} issue(s){}:", diagnostic_report.total_issues, auto_fix_text);
        
        for issue in &diagnostic_report.issues {
            let icon = match issue.severity {
                crosscell::diagnostics::IssueSeverity::Critical => "❌",
                crosscell::diagnostics::IssueSeverity::Warning => "⚠️",
                crosscell::diagnostics::IssueSeverity::Info => "ℹ️",
            };
            let fix_hint = if issue.auto_fixable { " → will auto-fix" } else { "" };
            println!("  {} {}{}", icon, issue.issue_type, fix_hint);
            
            if detailed {
                println!("     💡 {}", issue.suggestion);
            }
        }
        
        println!();
        println!("💡 Tip: Run 'crosscell convert' to automatically fix these issues");
    } else {
        println!("  ✓ No issues detected - ready for conversion");
    }

    // Detailed information
    if detailed {
        println!();
        println!("📋 Detailed Information");
        println!("  ├─ Total elements: {}", n_cells * n_genes);
        println!("  ├─ Memory estimate (dense): {:.2} MB", 
            (n_cells * n_genes * 8) as f64 / (1024.0 * 1024.0));
        println!("  ├─ Source format: {}", data.metadata.source_format);
        if let Some(ref version) = data.metadata.source_version {
            println!("  ├─ Source version: {}", version);
        }
        if let Some(ref active_assay) = data.metadata.active_assay {
            println!("  └─ Active assay: {}", active_assay);
        }
        
        // Show column names if detailed
        if !data.cell_metadata.columns.is_empty() {
            println!();
            println!("  Cell metadata columns:");
            for col in &data.cell_metadata.columns {
                println!("    - {}", col);
            }
        }
        
        if !data.gene_metadata.columns.is_empty() {
            println!();
            println!("  Gene metadata columns:");
            for col in &data.gene_metadata.columns {
                println!("    - {}", col);
            }
        }
    }

    // Write JSON report if requested
    if let Some(output_path) = output {
        info!("Writing report to: {}", output_path.display());
        
        let issues_json: Vec<_> = diagnostic_report.issues.iter().map(|issue| {
            json!({
                "severity": format!("{}", issue.severity),
                "type": format!("{}", issue.issue_type),
                "auto_fixable": issue.auto_fixable,
                "suggestion": issue.suggestion,
            })
        }).collect();
        
        let report = json!({
            "format": format,
            "path": input.to_str(),
            "size_bytes": file_size,
            "size_mb": size_mb,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "sparsity": sparsity,
            "has_embeddings": has_embeddings,
            "has_layers": has_layers,
            "has_spatial": data.spatial.is_some(),
            "compatibility_score": diagnostic_report.compatibility_score,
            "issues": {
                "total": diagnostic_report.total_issues,
                "critical": diagnostic_report.critical_count,
                "warning": diagnostic_report.warning_count,
                "info": diagnostic_report.info_count,
                "auto_fixable": diagnostic_report.auto_fixable_count,
                "details": issues_json,
            },
            "cell_metadata_columns": data.cell_metadata.columns,
            "gene_metadata_columns": data.gene_metadata.columns,
        });
        
        let json_str = serde_json::to_string_pretty(&report)?;
        fs::write(&output_path, json_str)?;
        
        println!();
        println!("📝 Report saved to: {}", output_path.display());
    }

    Ok(())
}

/// Calculate sparsity of expression matrix
fn calculate_sparsity(data: &SingleCellData) -> f64 {
    let n_cells = data.metadata.n_cells;
    let n_genes = data.metadata.n_genes;
    let total = n_cells * n_genes;
    
    if total == 0 {
        return 0.0;
    }
    
    match &data.expression {
        crosscell::ir::ExpressionMatrix::SparseCSR(m) => {
            let nnz = m.data.len();
            100.0 * (1.0 - (nnz as f64 / total as f64))
        }
        crosscell::ir::ExpressionMatrix::SparseCSC(m) => {
            let nnz = m.data.len();
            100.0 * (1.0 - (nnz as f64 / total as f64))
        }
        crosscell::ir::ExpressionMatrix::Dense(_) => 0.0,
        crosscell::ir::ExpressionMatrix::Lazy(m) => {
            let nnz = m.nnz.unwrap_or(0);
            100.0 * (1.0 - (nnz as f64 / total as f64))
        }
    }
}

/// Detect file format from extension
fn detect_format(path: &PathBuf) -> Result<String> {
    let extension = path
        .extension()
        .and_then(|s| s.to_str())
        .context("Failed to get file extension")?;

    match extension.to_lowercase().as_str() {
        "h5ad" | "h5" => Ok("AnnData (.h5ad)".to_string()),
        "rds" => Ok("Seurat (.rds)".to_string()),
        _ => anyhow::bail!("Unsupported file format: {}", extension),
    }
}
