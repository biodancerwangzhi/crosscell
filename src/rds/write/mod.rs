//! RDS 写入模块
//!
//! 提供 RDS 文件的写入功能。
//! 对应 rds2cpp 的写入相关文件。

// 工具函数
pub mod utils;
pub mod shared_info;

// 写入入口
pub mod rds;
pub mod object;

// 基础类型写入
pub mod single_string;
pub mod atomic;
pub mod attributes;

// 复合类型写入
pub mod pairlist;
pub mod list;
pub mod builtin;
pub mod expression;
pub mod language;
pub mod s4;

// 重新导出公共写入函数
pub use object::write_object;
pub use rds::{write_rds, write_rds_file};
pub use shared_info::SharedWriteInfo;
pub use single_string::write_single_string;
pub use atomic::{
    write_integer_body, write_logical_body, write_double_body,
    write_raw_body, write_complex_body, write_string_body,
};
pub use attributes::write_attributes;
pub use pairlist::write_pairlist;
pub use list::write_list_body;
pub use builtin::write_builtin_body;
pub use expression::write_expression_body;
pub use language::write_language;
pub use s4::write_s4;
