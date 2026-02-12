//! CrossCell CLI - Single-cell transcriptomics data format conversion tool
//!
//! Supports conversion between AnnData (.h5ad) and Seurat (.rds) formats

use clap::{Parser, Subcommand};
use env_logger::Env;
use log::{error, info};
use std::path::PathBuf;
use std::process;

mod commands;

use commands::{convert, inspect, simplify_seurat, validate};

/// CrossCell - Single-cell transcriptomics data format conversion tool
#[derive(Parser)]
#[command(name = "crosscell")]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    /// Subcommand to execute
    #[command(subcommand)]
    command: Commands,

    /// Enable verbose output mode
    #[arg(short, long, global = true)]
    verbose: bool,

    /// Enable quiet mode (errors only)
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,

    /// Number of threads (defaults to all CPU cores)
    #[arg(short, long, global = true)]
    threads: Option<usize>,

    /// Log file path
    #[arg(long, global = true)]
    log_file: Option<PathBuf>,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert file format (.h5ad ↔ .rds ↔ .loom)
    Convert {
        /// Input file path
        #[arg(short, long)]
        input: PathBuf,

        /// Output file path
        #[arg(short, long)]
        output: PathBuf,

        /// Target format (anndata, seurat, or loom)
        #[arg(short, long)]
        format: String,

        /// Validate data consistency after conversion
        #[arg(long)]
        validate: bool,

        /// Sparse matrix format (csr, csc, auto)
        #[arg(long, default_value = "auto")]
        sparse_format: String,

        /// Compression algorithm (none, gzip, lz4)
        #[arg(short, long, default_value = "gzip")]
        compression: String,

        /// Compression level (0-9)
        #[arg(long, default_value = "6")]
        compression_level: u32,

        /// Enable lazy loading for large datasets (reduces memory usage)
        #[arg(long)]
        lazy: bool,

        /// Chunk size for lazy loading (number of rows per chunk)
        #[arg(long, default_value = "10000")]
        chunk_size: usize,

        /// Enable streaming mode for large files (constant memory usage)
        #[arg(long)]
        streaming: bool,

        /// Estimate memory usage before conversion
        #[arg(long)]
        estimate_memory: bool,

        /// Preview changes without actually converting (dry run mode)
        #[arg(long)]
        dry_run: bool,

        /// Enable disk-backed mode for ultra-large datasets (10M+ cells)
        #[arg(long)]
        disk_backed: bool,

        /// Temporary directory for disk-backed mode
        #[arg(long)]
        temp_dir: Option<PathBuf>,

        /// Use direct reading for Seurat files (default, no R preprocessing needed)
        #[arg(long)]
        direct: bool,

        /// Use legacy simplify-first mode for Seurat files (requires R preprocessing)
        #[arg(long, conflicts_with = "direct")]
        simplify_first: bool,
    },

    /// List all supported formats
    Formats,

    /// Inspect file information and compatibility
    Inspect {
        /// Input file path
        #[arg(short, long)]
        input: PathBuf,

        /// Show detailed information
        #[arg(short, long)]
        detailed: bool,

        /// Output report path (JSON format)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Validate roundtrip conversion consistency
    Validate {
        /// Original file path
        #[arg(long)]
        original: PathBuf,

        /// Converted file path
        #[arg(long)]
        converted: PathBuf,

        /// Numerical tolerance
        #[arg(long, default_value = "1e-7")]
        tolerance: f64,

        /// Output validation report path
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Strict mode (fail on any difference)
        #[arg(long)]
        strict: bool,
    },

    /// Simplify Seurat object for conversion
    SimplifySeurat {
        /// Input Seurat RDS file
        #[arg(short, long)]
        input: PathBuf,

        /// Output simplified RDS file
        #[arg(short, long)]
        output: PathBuf,

        /// Keep graph objects
        #[arg(long)]
        keep_graphs: bool,

        /// Keep command history
        #[arg(long)]
        keep_commands: bool,

        /// Keep tool functions
        #[arg(long)]
        keep_tools: bool,
    },
}

fn main() {
    let cli = Cli::parse();

    // 设置日志级别
    let log_level = if cli.verbose {
        "debug"
    } else if cli.quiet {
        "error"
    } else {
        "info"
    };

    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();

    // 设置线程数
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or_else(|e| {
                error!("Failed to set thread count: {}", e);
                process::exit(1);
            });
        info!("Using {} threads", threads);
    }

    // 执行命令
    let result = match cli.command {
        Commands::Convert {
            input,
            output,
            format,
            validate,
            sparse_format,
            compression,
            compression_level,
            lazy,
            chunk_size,
            streaming,
            estimate_memory,
            dry_run,
            disk_backed,
            temp_dir,
            direct,
            simplify_first,
        } => convert::run(
            input,
            output,
            format,
            validate,
            sparse_format,
            compression,
            compression_level,
            lazy,
            chunk_size,
            streaming,
            estimate_memory,
            dry_run,
            disk_backed,
            temp_dir,
            direct,
            simplify_first,
        ),

        Commands::Formats => {
            convert::list_formats();
            Ok(())
        }

        Commands::Inspect {
            input,
            detailed,
            output,
        } => inspect::run(input, detailed, output),

        Commands::Validate {
            original,
            converted,
            tolerance,
            output,
            strict,
        } => validate::run(original, converted, tolerance, output, strict),

        Commands::SimplifySeurat {
            input,
            output,
            keep_graphs,
            keep_commands,
            keep_tools,
        } => simplify_seurat::run(input, output, keep_graphs, keep_commands, keep_tools),
    };

    // 处理结果
    if let Err(e) = result {
        error!("Error: {}", e);
        process::exit(1);
    }
}
