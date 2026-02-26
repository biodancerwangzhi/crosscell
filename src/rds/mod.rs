//! RDS - R Data Serialization 格式支持
//!
//! 本模块是 rds2cpp (C++) 库的 Rust 移植版本，提供完整的 RDS 文件解析和写入功能。
//!
//! ## 架构
//!
//! 模块结构与 rds2cpp 一一对应：
//! - `sexp_type` - R 数据类型枚举 (对应 SEXPType.hpp)
//! - `string_encoding` - 字符串编码类型 (对应 StringEncoding.hpp)
//! - `r_object` - R 对象数据结构 (对应 RObject.hpp)
//! - `symbol` - R 符号 (对应 Symbol.hpp)
//! - `environment` - R 环境 (对应 Environment.hpp)
//! - `external_pointer` - 外部指针 (对应 ExternalPointer.hpp)
//! - `rds_file` - RDS 文件结构 (对应 RdsFile.hpp)
//! - `parse/` - 解析模块
//! - `write/` - 写入模块
//!
//! ## 示例
//!
//! ```rust,no_run
//! use crosscell::rds::{parse_rds, write_rds, RdsFile};
//! use std::path::Path;
//!
//! // 读取 RDS 文件
//! let file = parse_rds(Path::new("data.rds"))?;
//! println!("Format version: {}", file.format_version);
//!
//! // 写入 RDS 文件
//! write_rds(&file, Path::new("output.rds"))?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

// 核心类型定义模块
pub mod environment;
pub mod error;
pub mod external_pointer;
pub mod r_object;
pub mod rds_file;
pub mod sexp_type;
pub mod string_encoding;
pub mod symbol;

// 解析和写入子模块
pub mod parse;
pub mod write;

// 属性测试模块（仅测试时编译）
#[cfg(test)]
mod proptest;

// 重新导出核心类型，便于外部使用
pub use environment::Environment;
pub use error::{RdsError, Result};
pub use external_pointer::ExternalPointer;
pub use r_object::{
    Attributes, BuiltInFunction, ComplexVector, DoubleVector, EnvironmentIndex, ExpressionVector,
    ExternalPointerIndex, GenericVector, IntegerVector, LanguageObject, LogicalVector, PairList,
    RObject, RawVector, S4Object, StringVector, SymbolIndex,
};
pub use rds_file::RdsFile;
pub use sexp_type::SEXPType;
pub use string_encoding::StringEncoding;
pub use symbol::Symbol;

// 便捷函数
use std::path::Path;

/// 读取 RDS 文件（便捷函数）
pub fn parse_rds<P: AsRef<Path>>(path: P) -> Result<RdsFile> {
    parse::rds::parse_rds_file(path, &parse::rds::ParseRdsOptions::default())
}

/// 写入 RDS 文件（便捷函数）
pub fn write_rds<P: AsRef<Path>>(file: &RdsFile, path: P) -> Result<()> {
    write::rds::write_rds_file_default(file, path)
}
