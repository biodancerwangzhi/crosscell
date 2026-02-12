//! Seurat 对象支持
//!
//! 提供从 Seurat RDS 文件提取数据到 CrossCell IR 格式的功能。
//!
//! ## 支持的 Seurat 版本
//!
//! - Seurat v5 (Assay5 结构)
//! - Seurat v3/v4 (传统 Assay 结构)
//!
//! ## 示例
//!
//! ```rust,no_run
//! use crosscell::seurat::seurat_rds_to_ir;
//!
//! let ir = seurat_rds_to_ir("data/seurat_object.rds")?;
//! println!("Loaded {} cells × {} genes", ir.metadata.n_cells, ir.metadata.n_genes);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod error;
pub mod extract;
pub mod seurat_to_ir;
pub mod ir_to_seurat;
pub mod assay5;
pub mod direct_read;

// 重新导出核心类型和函数
pub use error::SeuratError;
pub use seurat_to_ir::seurat_rds_to_ir;
pub use ir_to_seurat::{ir_to_seurat_rds, write_seurat_rds};
pub use assay5::{AssayType, detect_assay_type, extract_assay5_layers, extract_assay5_main_matrix};
pub use direct_read::{read_seurat_direct, DirectReadResult, SeuratVersion, SkippedComponents};
