//! SingleCellExperiment (SCE) 对象支持
//!
//! 提供从 SingleCellExperiment RDS 文件提取数据到 CrossCell IR 格式的功能。
//!
//! ## 支持的 SCE 版本
//!
//! - SingleCellExperiment (Bioconductor)
//! - SummarizedExperiment (基类)
//!
//! ## 示例
//!
//! ```rust,no_run
//! use crosscell::sce::sce_rds_to_ir;
//!
//! let ir = sce_rds_to_ir("data/sce_object.rds")?;
//! println!("Loaded {} cells × {} genes", ir.metadata.n_cells, ir.metadata.n_genes);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod error;
pub mod extract;
pub mod sce_to_ir;
pub mod ir_to_sce;

// 重新导出核心类型和函数
pub use error::SceError;
pub use sce_to_ir::sce_rds_to_ir;
pub use ir_to_sce::{ir_to_sce_rds, write_sce_rds};
pub use extract::is_sce_object;
