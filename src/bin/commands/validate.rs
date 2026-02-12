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
    println!("  Converted: {} ({})", converted.display(), converted_format);
    println!();

    // Read files and perform validation
    let report = if original_format == converted_format {
        // Same format - roundtrip validation
        match original_format.as_str() {
            "h5ad" => {
                let original_str = original.to_str().context("Invalid original path")?;
                let converted_str = converted.to_str().context("Invalid converted path")?;
                validation::roundtrip::validate_h5ad_roundtrip(original_str, converted_str, tolerance)
                    .map_err(|e| anyhow::anyhow!(e))?
            }
            "rds" => {
                let original_str = original.to_str().context("Invalid original path")?;
                let converted_str = converted.to_str().context("Invalid converted path")?;
                validation::roundtrip::validate_rds_roundtrip(original_str, converted_str, tolerance)
                    .map_err(|e| anyhow::anyhow!(e))?
            }
            _ => anyhow::bail!("Unsupported format for validation"),
        }
    } else {
        // Cross-format validation
        let (h5ad_path, rds_path) = if original_format == "h5ad" {
            (original.to_str().context("Invalid path")?, converted.to_str().context("Invalid path")?)
        } else {
            (converted.to_str().context("Invalid path")?, original.to_str().context("Invalid path")?)
        };
        
        validation::roundtrip::validate_cross_format(h5ad_path, rds_path, tolerance)
            .map_err(|e| anyhow::anyhow!(e))?
    };

    let elapsed = start_time.elapsed();

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
