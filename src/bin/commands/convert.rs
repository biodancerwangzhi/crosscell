//! Convert command implementation - format conversion
//!
//! Uses stage-based progress display instead of percentage bars,
//! which is more honest for operations where we can't predict exact progress.
//!
//! Supports multiple formats through the plugin system:
//! - AnnData (.h5ad)
//! - Seurat (.rds)
//! - Loom (.loom)

use anyhow::{Context, Result};
use crosscell::anndata::{
    streaming_convert, streaming_h5ad_to_rds, streaming_rds_to_h5ad, StreamingH5adReader,
};
use crosscell::diagnostics::{detect_issues, DataCleaner};
use crosscell::formats::create_default_registry;
use crosscell::ir::UnstructuredValue;
use crosscell::seurat::direct_read::{read_seurat_direct, SeuratVersion};
use crosscell::storage::{DiskBackedConfig, DiskBackedMatrix};
use crosscell::validation;
use log::{info, warn};
use std::path::PathBuf;
use std::time::Instant;

use super::progress::{format_bytes, StageProgress};

/// List all supported formats
pub fn list_formats() {
    let registry = create_default_registry();

    println!();
    println!("📋 Supported Formats");
    println!("─────────────────────────────────────────");
    println!();

    for info in registry.list_info() {
        let caps = match (info.can_read, info.can_write) {
            (true, true) => "read/write",
            (true, false) => "read only",
            (false, true) => "write only",
            (false, false) => "none",
        };
        let exts = info.extensions.join(", ");
        println!("  {:<12} ({})  [{}]", info.name, exts, caps);
        if !info.description.is_empty() {
            println!("               {}", info.description);
        }
    }
    println!();
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    input: PathBuf,
    output: PathBuf,
    format: String,
    validate: bool,
    _sparse_format: String,
    _compression: String,
    _compression_level: u32,
    lazy: bool,
    chunk_size: usize,
    streaming: bool,
    estimate_memory: bool,
    dry_run: bool,
    disk_backed: bool,
    temp_dir: Option<PathBuf>,
    _direct: bool,
    simplify_first: bool,
    normalize: bool,
    top_genes: Option<usize>,
    gene_id_column: Option<String>,
    keep_layers: bool,
) -> Result<()> {
    let start_time = Instant::now();
    let registry = create_default_registry();

    // Validate input file exists
    if !input.exists() {
        anyhow::bail!("Input file does not exist: {}", input.display());
    }

    // Detect input format using registry
    let input_converter = registry
        .detect(&input)
        .context(format!("Unsupported input format: {}", input.display()))?;
    let input_format = input_converter.name().to_string();

    // Get target format converter
    let target_format = format.to_lowercase();
    let output_converter = registry.get(&target_format).context(format!(
        "Unknown target format: {}. Use 'crosscell formats' to list supported formats.",
        format
    ))?;

    // Validate converters support required operations
    if !input_converter.can_read() {
        anyhow::bail!("Format '{}' does not support reading", input_format);
    }
    if !output_converter.can_write() {
        anyhow::bail!("Format '{}' does not support writing", target_format);
    }

    // Check if conversion is needed
    if input_format == target_format {
        warn!("Input and output formats are the same. No conversion needed.");
        return Ok(());
    }

    // Auto-detect large files and suggest streaming mode
    let file_size = std::fs::metadata(&input)?.len();
    let should_use_streaming =
        streaming || should_auto_enable_streaming(&input, file_size, &input_format);

    // If streaming mode is enabled, use streaming conversion
    if should_use_streaming {
        return run_streaming_conversion(
            input,
            output,
            &input_format,
            &target_format,
            chunk_size,
            validate,
            dry_run,
            file_size,
            start_time,
        );
    }

    // Determine if we should use direct reading for Seurat files
    // Default is direct reading (unless --simplify-first is specified)
    let use_direct_read = input_format == "seurat" && !simplify_first;

    // Print header
    println!();
    println!("🔄 CrossCell Format Conversion");
    println!("   {} → {}", input.display(), output.display());
    println!(
        "   Format: {} → {}",
        input_format.to_uppercase(),
        target_format.to_uppercase()
    );

    if use_direct_read {
        println!("   Mode: Direct read (no R preprocessing)");
    } else if input_format == "seurat" && simplify_first {
        println!("   Mode: Simplify-first (legacy, requires R preprocessing)");
    }

    if disk_backed {
        let temp = temp_dir
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| std::env::temp_dir().display().to_string());
        println!("   Disk-backed: {} (temp: {})", "enabled", temp);
        println!("   Chunk size: {} rows", chunk_size);
    } else if lazy {
        println!("   Lazy loading: enabled (chunk size: {})", chunk_size);
    }

    if dry_run {
        println!("   Dry run: enabled (no files will be modified)");
    }
    println!();

    // Estimate memory if requested
    if estimate_memory {
        print_memory_estimate(&input, lazy, disk_backed, chunk_size)?;
    }

    info!("Converting {} to {}", input.display(), output.display());
    info!("Target format: {}", format);

    if use_direct_read {
        info!("Using direct read mode for Seurat file");
    }

    if disk_backed {
        info!("Disk-backed mode enabled with chunk size: {}", chunk_size);
    } else if lazy {
        info!("Lazy loading enabled with chunk size: {}", chunk_size);
    }

    // Define stages based on conversion type
    let stages = if use_direct_read {
        if validate {
            vec![
                "Reading Seurat object directly",
                "Running diagnostics",
                "Converting format",
                "Writing output file",
                "Validating result",
            ]
        } else {
            vec![
                "Reading Seurat object directly",
                "Running diagnostics",
                "Converting format",
                "Writing output file",
            ]
        }
    } else if disk_backed {
        if validate {
            vec![
                "Reading input file",
                "Running diagnostics",
                "Writing to disk cache",
                "Converting format",
                "Writing output file",
                "Validating result",
                "Cleaning up",
            ]
        } else {
            vec![
                "Reading input file",
                "Running diagnostics",
                "Writing to disk cache",
                "Converting format",
                "Writing output file",
                "Cleaning up",
            ]
        }
    } else if validate {
        vec![
            "Reading input file",
            "Running diagnostics",
            "Converting format",
            "Writing output file",
            "Validating result",
        ]
    } else {
        vec![
            "Reading input file",
            "Running diagnostics",
            "Converting format",
            "Writing output file",
        ]
    };

    let mut progress = StageProgress::new(stages, true);

    // Stage 1: Read input file
    progress.next_stage();

    let mut ir_data = if use_direct_read {
        // Try direct reading for Seurat files first
        progress.update_message("parsing Seurat object structure");

        match read_seurat_direct(&input, keep_layers) {
            Ok(result) => {
                // Show version and skipped components
                let version_str = match result.version {
                    SeuratVersion::V3 => "V3",
                    SeuratVersion::V4 => "V4",
                    SeuratVersion::V5 => "V5",
                    SeuratVersion::Unknown => "Unknown",
                };

                let detail = format!(
                    "Seurat {} - {} cells × {} genes",
                    version_str, result.data.metadata.n_cells, result.data.metadata.n_genes
                );
                progress.complete_stage(Some(&detail));

                // Report skipped components
                if result.skipped.has_skipped() {
                    let skipped = &result.skipped;
                    if !skipped.closures.is_empty() {
                        progress
                            .verbose(&format!("⚠ Skipped {} closure(s)", skipped.closures.len()));
                    }
                    if !skipped.environments.is_empty() {
                        progress.verbose(&format!(
                            "⚠ Skipped {} environment(s)",
                            skipped.environments.len()
                        ));
                    }
                    if !skipped.external_pointers.is_empty() {
                        progress.verbose(&format!(
                            "⚠ Skipped {} external pointer(s)",
                            skipped.external_pointers.len()
                        ));
                    }
                    if !skipped.bytecodes.is_empty() {
                        progress.verbose(&format!(
                            "⚠ Skipped {} bytecode(s)",
                            skipped.bytecodes.len()
                        ));
                    }
                }

                result.data
            }
            Err(direct_err) => {
                // Fall back to plugin system for simplified files
                eprintln!("Direct read failed: {}", direct_err);
                progress.verbose("Direct read failed, trying simplified format...");

                let data = input_converter
                    .read(&input)
                    .map_err(|e| anyhow::anyhow!("Failed to read {} file: {}", input_format, e))?;

                let detail = format!(
                    "{} cells × {} genes (simplified)",
                    data.metadata.n_cells, data.metadata.n_genes
                );
                progress.complete_stage(Some(&detail));

                data
            }
        }
    } else {
        // Use plugin system for other formats
        progress.update_message(&format!(
            "parsing {} structure",
            input_format.to_uppercase()
        ));

        let data = input_converter
            .read(&input)
            .map_err(|e| anyhow::anyhow!("Failed to read {} file: {}", input_format, e))?;

        let detail = format!(
            "{} cells × {} genes",
            data.metadata.n_cells, data.metadata.n_genes
        );
        progress.complete_stage(Some(&detail));

        data
    };

    // Stage 2: Run diagnostics
    progress.next_stage();
    progress.update_message("analyzing data structure");

    let diagnostic_report = detect_issues(&ir_data);
    let data_cleaner = DataCleaner::with_defaults();
    let cleaning_preview = data_cleaner.preview(&ir_data);

    if diagnostic_report.has_issues() {
        let detail = format!("{} issue(s) found", diagnostic_report.total_issues);
        progress.complete_stage(Some(&detail));

        // Show issues
        for issue in &diagnostic_report.issues {
            let icon = match issue.severity {
                crosscell::diagnostics::IssueSeverity::Critical => "❌",
                crosscell::diagnostics::IssueSeverity::Warning => "⚠️",
                crosscell::diagnostics::IssueSeverity::Info => "ℹ️",
            };
            progress.verbose(&format!("{} {}", icon, issue.issue_type));
        }
    } else {
        progress.complete_stage(Some("no issues"));
    }

    // Handle dry run mode
    if dry_run {
        progress.finish();
        print_dry_run_summary(
            &input,
            &output,
            &input_format,
            &target_format,
            &ir_data,
            &diagnostic_report,
            &cleaning_preview,
        )?;
        return Ok(());
    }

    // Auto-fix issues if needed
    if diagnostic_report.auto_fixable_count > 0 {
        progress.update_message("auto-fixing issues");
        let cleaning_report = data_cleaner.clean(&mut ir_data);

        if cleaning_report.has_modifications() {
            progress.verbose(&format!(
                "Fixed {} issue(s)",
                cleaning_report.total_modifications
            ));

            // Store audit log
            let audit_log = cleaning_report.to_audit_log();
            if let Some(ref mut uns) = ir_data.unstructured {
                if let Ok(json_str) = audit_log.to_json() {
                    uns.insert(
                        "crosscell_audit_log".to_string(),
                        UnstructuredValue::String(json_str),
                    );
                }
            }
        }
    }

    // Disk-backed mode: write to disk cache first
    let _disk_matrix: Option<DiskBackedMatrix> = if disk_backed {
        progress.next_stage();
        progress.update_message("caching matrix to disk");

        let config = DiskBackedConfig::new()
            .with_chunk_size(chunk_size)
            .with_compression(true)
            .with_cleanup(true);

        let config = if let Some(ref dir) = temp_dir {
            config.with_temp_dir(dir)
        } else {
            config
        };

        // Check if expression matrix is sparse
        if ir_data.expression.is_sparse() {
            match DiskBackedMatrix::from_matrix(&ir_data.expression, config) {
                Ok(dm) => {
                    let disk_usage = dm.disk_usage();
                    let num_chunks = dm.num_chunks();
                    progress.complete_stage(Some(&format!(
                        "{} chunks, {}",
                        num_chunks,
                        format_bytes(disk_usage as usize)
                    )));
                    Some(dm)
                }
                Err(e) => {
                    progress.complete_stage(Some("skipped (not sparse)"));
                    warn!("Disk-backed mode skipped: {}", e);
                    None
                }
            }
        } else {
            progress.complete_stage(Some("skipped (dense matrix)"));
            warn!("Disk-backed mode only supports sparse matrices");
            None
        }
    } else {
        None
    };

    // Apply AI-Ready transforms (normalize → top-genes → gene-id-column)
    if normalize {
        info!("Applying library-size normalization + log1p");
        ir_data.expression = crosscell::transform::normalize_library_size(&ir_data.expression)
            .map_err(|e| anyhow::anyhow!("Normalization failed: {}", e))?;
        // Also normalize layers if present
        if let Some(ref mut layers) = ir_data.layers {
            for (name, matrix) in layers.iter_mut() {
                match crosscell::transform::normalize_library_size(matrix) {
                    Ok(normalized) => *matrix = normalized,
                    Err(e) => warn!("Skipping normalization for layer '{}': {}", name, e),
                }
            }
        }
    }

    if let Some(n) = top_genes {
        info!("Selecting top {} variable genes", n);
        crosscell::transform::select_top_variable_genes(&mut ir_data, n)
            .map_err(|e| anyhow::anyhow!("Gene selection failed: {}", e))?;
    }

    if let Some(ref col) = gene_id_column {
        info!("Applying gene ID column: {}", col);
        crosscell::transform::apply_gene_id_column(&mut ir_data, col)
            .map_err(|e| anyhow::anyhow!("Gene ID replacement failed: {}", e))?;
    }

    // Stage 3: Convert format (prepare for writing)
    progress.next_stage();
    progress.update_message(&format!(
        "preparing {} output",
        target_format.to_uppercase()
    ));
    progress.complete_stage(Some(&format!("{} format", target_format.to_uppercase())));

    // Stage 4: Write output file using plugin system
    progress.next_stage();
    progress.update_message(&format!("writing {} file", target_format.to_uppercase()));

    output_converter
        .write(&ir_data, &output)
        .map_err(|e| anyhow::anyhow!("Failed to write {} file: {}", target_format, e))?;

    let file_size = std::fs::metadata(&output)
        .map(|m| format_bytes(m.len() as usize))
        .unwrap_or_else(|_| "unknown".to_string());
    progress.complete_stage(Some(&file_size));

    // Stage 5: Validate (if requested)
    if validate {
        progress.next_stage();
        progress.update_message("comparing input and output");

        // Re-read the output file using the plugin system
        let converted = output_converter
            .read(&output)
            .map_err(|e| anyhow::anyhow!("Failed to read converted file for validation: {}", e))?;

        let validation_result =
            validation::roundtrip::validate_roundtrip(&ir_data, &converted, 1e-7);

        if validation_result.passed() {
            progress.complete_stage(Some("passed"));
        } else {
            progress.complete_stage(Some("warnings"));
            progress.println(&format!("\n{}", validation_result.summary()));
        }
    }

    // Cleanup stage for disk-backed mode
    if disk_backed {
        progress.next_stage();
        progress.update_message("removing temp files");
        // DiskBackedMatrix will cleanup on drop
        drop(_disk_matrix);
        progress.complete_stage(Some("done"));
    }

    // Finish and show summary
    progress.finish();

    let elapsed = start_time.elapsed();
    info!("Conversion completed in {:.2}s", elapsed.as_secs_f64());

    // Final summary
    println!();
    println!("📁 Output: {}", output.display());

    Ok(())
}

/// Print memory estimation
fn print_memory_estimate(
    input: &PathBuf,
    lazy: bool,
    disk_backed: bool,
    chunk_size: usize,
) -> Result<()> {
    let file_size = std::fs::metadata(input)?.len();
    let estimated_full = file_size as f64 * 3.0;
    let estimated_lazy = file_size as f64 * 0.5;
    let estimated_disk_backed =
        (chunk_size as f64 * 30000.0 * 8.0 * 2.0) / (1024.0 * 1024.0 * 1024.0); // ~2 chunks in memory

    println!("📊 Memory Estimation:");
    println!("   File size: {}", format_bytes(file_size as usize));
    println!("   Full load: ~{}", format_bytes(estimated_full as usize));

    if disk_backed {
        println!("   Disk-backed: ~{:.1} GB (peak)", estimated_disk_backed);
        println!("   Chunk size: {} rows", chunk_size);
        println!(
            "   Disk usage: ~{}",
            format_bytes((file_size as f64 * 0.7) as usize)
        );
    } else if lazy {
        println!(
            "   Lazy load: ~{} (peak)",
            format_bytes(estimated_lazy as usize)
        );
        println!("   Chunk size: {} rows", chunk_size);
    }

    if estimated_full > 8.0 * 1024.0 * 1024.0 * 1024.0 && !lazy && !disk_backed {
        println!();
        println!("💡 Tip: Use --disk-backed flag for ultra-large datasets (10M+ cells)");
        println!("        Use --lazy flag to reduce memory usage for large datasets");
    }
    println!();

    Ok(())
}

/// Print dry run summary
fn print_dry_run_summary(
    input: &PathBuf,
    output: &PathBuf,
    input_format: &str,
    target_format: &str,
    ir_data: &crosscell::ir::SingleCellData,
    diagnostic_report: &crosscell::diagnostics::DiagnosticReport,
    cleaning_preview: &crosscell::diagnostics::CleaningReport,
) -> Result<()> {
    println!();
    println!("📋 Dry Run Summary");
    println!("─────────────────────────────────────────");
    println!();
    println!(
        "Input:  {} ({})",
        input.display(),
        input_format.to_uppercase()
    );
    println!(
        "Output: {} ({})",
        output.display(),
        target_format.to_uppercase()
    );
    println!("Cells:  {}", ir_data.metadata.n_cells);
    println!("Genes:  {}", ir_data.metadata.n_genes);
    println!();

    if diagnostic_report.has_issues() {
        println!("🔍 Issues detected: {}", diagnostic_report.total_issues);
        for issue in &diagnostic_report.issues {
            let icon = match issue.severity {
                crosscell::diagnostics::IssueSeverity::Critical => "❌",
                crosscell::diagnostics::IssueSeverity::Warning => "⚠️",
                crosscell::diagnostics::IssueSeverity::Info => "ℹ️",
            };
            println!("   {} {}", icon, issue.issue_type);
        }
        println!();
    }

    if cleaning_preview.has_modifications() {
        println!(
            "🔧 Planned modifications: {}",
            cleaning_preview.total_modifications
        );

        if !cleaning_preview.columns_renamed.is_empty() {
            println!(
                "   Rename {} column(s)",
                cleaning_preview.columns_renamed.len()
            );
        }
        if !cleaning_preview.duplicates_resolved.is_empty() {
            println!(
                "   Deduplicate {} name(s)",
                cleaning_preview.duplicates_resolved.len()
            );
        }
        if !cleaning_preview.names_truncated.is_empty() {
            println!(
                "   Truncate {} name(s)",
                cleaning_preview.names_truncated.len()
            );
        }
        if !cleaning_preview.columns_dropped.is_empty() {
            println!(
                "   Drop {} column(s)",
                cleaning_preview.columns_dropped.len()
            );
        }
        println!();
    } else {
        println!("✅ No modifications needed - data is clean");
        println!();
    }

    // Estimate output size
    let input_size = std::fs::metadata(input)?.len();
    let estimated_output = (input_size as f64 * 0.95) as u64;

    println!("📊 Estimates:");
    println!("   Input size:  {}", format_bytes(input_size as usize));
    println!(
        "   Output size: ~{}",
        format_bytes(estimated_output as usize)
    );
    println!(
        "   Compatibility: {}/100",
        diagnostic_report.compatibility_score
    );
    println!();

    println!("💡 To proceed with conversion, run without --dry-run:");
    println!(
        "   crosscell convert -i {} -o {} -f {}",
        input.display(),
        output.display(),
        target_format
    );

    Ok(())
}

/// Check if streaming mode should be auto-enabled based on file size
fn should_auto_enable_streaming(_input: &PathBuf, file_size: u64, input_format: &str) -> bool {
    // Auto-enable streaming for files > 5GB
    const STREAMING_THRESHOLD: u64 = 5 * 1024 * 1024 * 1024; // 5 GB

    if file_size > STREAMING_THRESHOLD {
        // Only auto-enable for supported formats
        matches!(input_format, "anndata" | "seurat")
    } else {
        false
    }
}

/// Estimate cell count from file size (rough approximation)
fn estimate_cell_count(file_size: u64, input_format: &str) -> usize {
    // Rough estimate: ~3KB per cell for sparse data
    let bytes_per_cell = match input_format {
        "anndata" => 3000,
        "seurat" => 4000,
        _ => 3500,
    };
    (file_size / bytes_per_cell) as usize
}

/// Estimate memory usage for streaming mode
fn estimate_streaming_memory(n_genes: usize, chunk_size: usize) -> usize {
    // Memory = chunk_size × n_genes × 8 bytes (f64) × 2 (read + write buffers)
    chunk_size * n_genes * 8 * 2
}

/// Run streaming conversion
fn run_streaming_conversion(
    input: PathBuf,
    output: PathBuf,
    input_format: &str,
    target_format: &str,
    chunk_size: usize,
    validate: bool,
    dry_run: bool,
    file_size: u64,
    start_time: Instant,
) -> Result<()> {
    let estimated_cells = estimate_cell_count(file_size, input_format);

    // Print header
    println!();
    println!("🔄 CrossCell Streaming Conversion");
    println!("   {} → {}", input.display(), output.display());
    println!(
        "   Format: {} → {}",
        input_format.to_uppercase(),
        target_format.to_uppercase()
    );
    println!("   Mode: Streaming (constant memory)");
    println!("   Chunk size: {} cells", chunk_size);
    println!();

    // Show file info and estimates
    println!("📊 File Information:");
    println!("   File size: {}", format_bytes(file_size as usize));
    println!("   Estimated cells: ~{}", format_number(estimated_cells));

    // Estimate memory based on typical gene count (30000)
    let estimated_memory = estimate_streaming_memory(30000, chunk_size);
    println!("   Estimated memory: ~{}", format_bytes(estimated_memory));
    println!();

    if dry_run {
        println!("📋 Dry Run - No files will be modified");
        println!();
        println!("💡 To proceed with streaming conversion:");
        println!(
            "   crosscell convert -i {} -o {} -f {} --streaming --chunk-size {}",
            input.display(),
            output.display(),
            target_format,
            chunk_size
        );
        return Ok(());
    }

    // Determine conversion type and run
    match (input_format, target_format) {
        ("anndata", "anndata") => {
            // H5AD → H5AD (re-chunking or compression)
            run_h5ad_to_h5ad_streaming(&input, &output, chunk_size)?;
        }
        ("anndata", "seurat") => {
            // H5AD → RDS
            run_h5ad_to_rds_streaming(&input, &output, chunk_size)?;
        }
        ("seurat", "anndata") => {
            // RDS → H5AD
            run_rds_to_h5ad_streaming(&input, &output, chunk_size)?;
        }
        _ => {
            anyhow::bail!(
                "Streaming mode not supported for {} → {} conversion. \
                 Supported: anndata↔seurat, anndata→anndata",
                input_format,
                target_format
            );
        }
    }

    let elapsed = start_time.elapsed();

    // Show completion summary
    println!();
    println!("✅ Streaming conversion completed!");
    println!("   Time: {:.2}s", elapsed.as_secs_f64());

    if let Ok(output_meta) = std::fs::metadata(&output) {
        println!(
            "   Output size: {}",
            format_bytes(output_meta.len() as usize)
        );
    }

    println!("   📁 Output: {}", output.display());

    // Validation in streaming mode is limited
    if validate {
        println!();
        println!("⚠️  Note: Full validation not available in streaming mode.");
        println!("   Use 'crosscell validate' command for detailed validation.");
    }

    Ok(())
}

/// Run H5AD → H5AD streaming conversion
fn run_h5ad_to_h5ad_streaming(input: &PathBuf, output: &PathBuf, chunk_size: usize) -> Result<()> {
    println!("🔄 Starting H5AD → H5AD streaming conversion...");
    println!();

    // Get metadata first
    let reader = StreamingH5adReader::new(input, chunk_size)
        .map_err(|e| anyhow::anyhow!("Failed to open input file: {}", e))?;

    let metadata = reader.metadata();
    let total_cells = metadata.n_cells;
    let total_chunks = reader.total_chunks();

    println!("   Cells: {}", format_number(total_cells));
    println!("   Genes: {}", format_number(metadata.n_genes));
    println!("   Chunks: {}", total_chunks);
    println!();

    // Run streaming conversion with progress
    let mut processed = 0usize;
    let progress_callback = |current: usize, total: usize| {
        if current > processed {
            processed = current;
            let pct = (current as f64 / total as f64 * 100.0) as usize;
            print!(
                "\r   Progress: {}% ({}/{} cells)",
                pct,
                format_number(current),
                format_number(total)
            );
            std::io::Write::flush(&mut std::io::stdout()).ok();
        }
    };

    streaming_convert(input, output, chunk_size, Some(progress_callback))
        .map_err(|e| anyhow::anyhow!("Streaming conversion failed: {}", e))?;

    println!();
    Ok(())
}

/// Run H5AD → RDS streaming conversion
fn run_h5ad_to_rds_streaming(input: &PathBuf, output: &PathBuf, chunk_size: usize) -> Result<()> {
    println!("🔄 Starting H5AD → RDS streaming conversion...");
    println!();

    // Get metadata first
    let reader = StreamingH5adReader::new(input, chunk_size)
        .map_err(|e| anyhow::anyhow!("Failed to open input file: {}", e))?;

    let metadata = reader.metadata();
    let total_cells = metadata.n_cells;

    println!("   Cells: {}", format_number(total_cells));
    println!("   Genes: {}", format_number(metadata.n_genes));
    println!();

    // Run streaming conversion with progress
    let mut processed = 0usize;
    let progress_callback = |current: usize, total: usize| {
        if current > processed {
            processed = current;
            let pct = (current as f64 / total as f64 * 100.0) as usize;
            print!(
                "\r   Progress: {}% ({}/{} cells)",
                pct,
                format_number(current),
                format_number(total)
            );
            std::io::Write::flush(&mut std::io::stdout()).ok();
        }
    };

    streaming_h5ad_to_rds(input, output, chunk_size, Some(progress_callback))
        .map_err(|e| anyhow::anyhow!("Streaming conversion failed: {}", e))?;

    println!();
    Ok(())
}

/// Run RDS → H5AD streaming conversion
fn run_rds_to_h5ad_streaming(input: &PathBuf, output: &PathBuf, chunk_size: usize) -> Result<()> {
    println!("🔄 Starting RDS → H5AD streaming conversion...");
    println!();
    println!("   Note: RDS format requires full file load, streaming applies to output.");
    println!();

    // Run streaming conversion with progress
    let mut processed = 0usize;
    let progress_callback = |current: usize, total: usize| {
        if current > processed {
            processed = current;
            let pct = (current as f64 / total as f64 * 100.0) as usize;
            print!(
                "\r   Progress: {}% ({}/{} cells)",
                pct,
                format_number(current),
                format_number(total)
            );
            std::io::Write::flush(&mut std::io::stdout()).ok();
        }
    };

    streaming_rds_to_h5ad(input, output, chunk_size, Some(progress_callback))
        .map_err(|e| anyhow::anyhow!("Streaming conversion failed: {}", e))?;

    println!();
    Ok(())
}

/// Format a number with thousands separators
fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    let chars: Vec<char> = s.chars().collect();

    for (i, c) in chars.iter().enumerate() {
        if i > 0 && (chars.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(*c);
    }

    result
}
