//! Validate command implementation - verify roundtrip conversion consistency

use anyhow::{Context, Result};
use crosscell::validation;
use log::info;
use serde_json::json;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

pub fn run(
    original: PathBuf,
    converted: PathBuf,
    tolerance: f64,
    output: Option<PathBuf>,
    strict: bool,
    cluster_column: Option<String>,
    validate_schema: Option<String>,
) -> Result<()> {
    let start_time = Instant::now();

    info!("Validating conversion");
    info!("Original: {}", original.display());
    info!("Converted: {}", converted.display());
    info!("Tolerance: {}", tolerance);
    info!("Strict mode: {}", strict);

    // Validate files exist
    if !original.exists() {
        anyhow::bail!("Original file does not exist: {}", original.display());
    }
    if !converted.exists() {
        anyhow::bail!("Converted file does not exist: {}", converted.display());
    }

    // Detect file formats
    let original_format = detect_format(&original)?;
    let converted_format = detect_format(&converted)?;

    println!("🔍 Validating Conversion");
    println!("  Original:  {} ({})", original.display(), original_format);
    println!(
        "  Converted: {} ({})",
        converted.display(),
        converted_format
    );
    println!();

    // Read files and perform validation
    let report = if original_format == converted_format {
        // Same format - roundtrip validation
        match original_format.as_str() {
            "h5ad" => {
                let original_str = original.to_str().context("Invalid original path")?;
                let converted_str = converted.to_str().context("Invalid converted path")?;
                validation::roundtrip::validate_h5ad_roundtrip(
                    original_str,
                    converted_str,
                    tolerance,
                )
                .map_err(|e| anyhow::anyhow!(e))?
            }
            "rds" => {
                let original_str = original.to_str().context("Invalid original path")?;
                let converted_str = converted.to_str().context("Invalid converted path")?;
                validation::roundtrip::validate_rds_roundtrip(
                    original_str,
                    converted_str,
                    tolerance,
                )
                .map_err(|e| anyhow::anyhow!(e))?
            }
            _ => anyhow::bail!("Unsupported format for validation"),
        }
    } else {
        // Cross-format validation
        let (h5ad_path, rds_path) = if original_format == "h5ad" {
            (
                original.to_str().context("Invalid path")?,
                converted.to_str().context("Invalid path")?,
            )
        } else {
            (
                converted.to_str().context("Invalid path")?,
                original.to_str().context("Invalid path")?,
            )
        };

        validation::roundtrip::validate_cross_format(h5ad_path, rds_path, tolerance)
            .map_err(|e| anyhow::anyhow!(e))?
    };

    let elapsed = start_time.elapsed();

    // Calculate ARI/NMI if cluster_column is specified
    let (ari_result, nmi_result) = if let Some(ref col_name) = cluster_column {
        info!("Calculating ARI/NMI for cluster column: {}", col_name);

        // Load data to extract cluster labels
        let (orig_data, conv_data) = match (original_format.as_str(), converted_format.as_str()) {
            ("h5ad", "h5ad") => {
                let orig = crosscell::anndata::reader::read_h5ad(
                    original.to_str().context("Invalid path")?,
                )
                .context("Failed to read original H5AD")?;
                let conv = crosscell::anndata::reader::read_h5ad(
                    converted.to_str().context("Invalid path")?,
                )
                .context("Failed to read converted H5AD")?;
                (orig, conv)
            }
            ("rds", "rds") => {
                let orig = crosscell::seurat::seurat_to_ir::seurat_rds_to_ir(
                    original.to_str().context("Invalid path")?,
                )
                .context("Failed to read original RDS")?;
                let conv = crosscell::seurat::seurat_to_ir::seurat_rds_to_ir(
                    converted.to_str().context("Invalid path")?,
                )
                .context("Failed to read converted RDS")?;
                (orig, conv)
            }
            ("h5ad", "rds") | ("rds", "h5ad") => {
                let (h5ad_p, rds_p) = if original_format == "h5ad" {
                    (&original, &converted)
                } else {
                    (&converted, &original)
                };
                let h5ad_data =
                    crosscell::anndata::reader::read_h5ad(h5ad_p.to_str().context("Invalid path")?)
                        .context("Failed to read H5AD")?;
                let rds_data = crosscell::seurat::seurat_to_ir::seurat_rds_to_ir(
                    rds_p.to_str().context("Invalid path")?,
                )
                .context("Failed to read RDS")?;
                if original_format == "h5ad" {
                    (h5ad_data, rds_data)
                } else {
                    (rds_data, h5ad_data)
                }
            }
            _ => anyhow::bail!("Unsupported format combination"),
        };

        // Extract labels from obs
        let orig_labels = extract_string_labels(&orig_data.cell_metadata, col_name)
            .with_context(|| format!("Failed to extract '{}' from original obs", col_name))?;
        let conv_labels = extract_string_labels(&conv_data.cell_metadata, col_name)
            .with_context(|| format!("Failed to extract '{}' from converted obs", col_name))?;

        let ari = crosscell::validation::calculate_ari(&orig_labels, &conv_labels);
        let nmi = crosscell::validation::calculate_nmi(&orig_labels, &conv_labels);
        (Some(ari), Some(nmi))
    } else {
        (None, None)
    };

    // Schema validation if requested
    let schema_result = if let Some(ref schema_version) = validate_schema {
        info!("Validating against schema: {}", schema_version);

        // Load data for schema validation (reuse if already loaded for ARI/NMI, otherwise load fresh)
        let data = load_single_file_data(&original, &original_format)?;
        let result = crosscell::validation::validate_cellxgene_schema(&data, schema_version);
        Some(result)
    } else {
        None
    };

    // Display results
    if report.passed() {
        println!("✅ Validation Passed");
    } else {
        if strict {
            println!("❌ Validation Failed (Strict Mode)");
        } else {
            println!("⚠️  Validation Completed with Warnings");
        }
    }

    println!();
    println!("{}", report.summary());

    if !report.passed() || strict {
        println!();
        println!("📋 Detailed Report:");
        println!("{}", report.detailed_report());
    }

    println!();
    println!("⚡ Performance Metrics");
    println!("  ├─ Validation time: {:.2}s", elapsed.as_secs_f64());
    println!("  └─ Tolerance: {}", tolerance);

    // Display ARI/NMI if calculated
    if let (Some(ari), Some(nmi)) = (ari_result, nmi_result) {
        println!();
        println!(
            "🔬 Cluster Metrics (column: {})",
            cluster_column.as_deref().unwrap_or("")
        );
        println!("  ├─ ARI: {:.6}", ari);
        println!("  └─ NMI: {:.6}", nmi);
    }

    // Display schema validation result
    if let Some(ref sr) = schema_result {
        println!();
        println!("{}", sr.summary());
    }

    // Write validation report if requested
    if let Some(output_path) = output {
        info!("Writing validation report to: {}", output_path.display());

        let json_report = json!({
            "passed": report.passed(),
            "original_file": original.to_str(),
            "converted_file": converted.to_str(),
            "tolerance": tolerance,
            "strict_mode": strict,
            "validation_time_seconds": elapsed.as_secs_f64(),
            "summary": report.summary(),
            "detailed_report": report.detailed_report(),
            "ari": ari_result,
            "nmi": nmi_result,
            "schema_validation": schema_result.as_ref().map(|sr| json!({
                "passed": sr.passed,
                "missing_obs_fields": sr.missing_obs_fields,
                "missing_var_fields": sr.missing_var_fields,
            })),
        });

        let json_str = serde_json::to_string_pretty(&json_report)?;
        fs::write(&output_path, json_str)?;

        println!();
        println!("📝 Validation report saved to: {}", output_path.display());
    }

    // Exit with error if validation failed in strict mode
    if strict && !report.passed() {
        anyhow::bail!("Validation failed in strict mode");
    }

    Ok(())
}

/// Detect file format from extension
fn detect_format(path: &PathBuf) -> Result<String> {
    let extension = path
        .extension()
        .and_then(|s| s.to_str())
        .context("Failed to get file extension")?;

    match extension.to_lowercase().as_str() {
        "h5ad" | "h5" => Ok("h5ad".to_string()),
        "rds" => Ok("rds".to_string()),
        _ => anyhow::bail!("Unsupported file format: {}", extension),
    }
}

/// Load a single file into SingleCellData IR
fn load_single_file_data(path: &PathBuf, format: &str) -> Result<crosscell::ir::SingleCellData> {
    let path_str = path.to_str().context("Invalid file path")?;
    match format {
        "h5ad" => {
            crosscell::anndata::reader::read_h5ad(path_str).context("Failed to read H5AD file")
        }
        "rds" => crosscell::seurat::seurat_to_ir::seurat_rds_to_ir(path_str)
            .context("Failed to read RDS file"),
        _ => anyhow::bail!("Unsupported format: {}", format),
    }
}

/// Extract string labels from a DataFrame column (supports String, Dictionary/Categorical)
fn extract_string_labels(df: &crosscell::ir::DataFrame, column_name: &str) -> Result<Vec<String>> {
    use arrow::array::*;
    use arrow::datatypes::DataType;

    let array = df
        .column(column_name)
        .ok_or_else(|| anyhow::anyhow!("Column '{}' not found in obs metadata", column_name))?;

    match array.data_type() {
        DataType::Utf8 => {
            let arr = array
                .as_any()
                .downcast_ref::<StringArray>()
                .context("Failed to downcast to StringArray")?;
            Ok((0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        String::new()
                    } else {
                        arr.value(i).to_string()
                    }
                })
                .collect())
        }
        DataType::LargeUtf8 => {
            let arr = array
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .context("Failed to downcast to LargeStringArray")?;
            Ok((0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        String::new()
                    } else {
                        arr.value(i).to_string()
                    }
                })
                .collect())
        }
        DataType::Dictionary(_, _) => {
            // Try Int32 dictionary first (most common for categorical)
            if let Some(dict) = array
                .as_any()
                .downcast_ref::<DictionaryArray<arrow::datatypes::Int32Type>>()
            {
                let values = dict.values();
                if let Some(str_values) = values.as_any().downcast_ref::<StringArray>() {
                    return Ok((0..dict.len())
                        .map(|i| {
                            if dict.is_null(i) {
                                String::new()
                            } else {
                                let key = dict.keys().value(i) as usize;
                                str_values.value(key).to_string()
                            }
                        })
                        .collect());
                }
            }
            anyhow::bail!("Unsupported dictionary type for column '{}'", column_name)
        }
        dt => anyhow::bail!(
            "Unsupported data type {:?} for cluster column '{}'",
            dt,
            column_name
        ),
    }
}
