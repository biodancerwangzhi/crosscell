//! RDS 写入模块
//!
//! 提供 RDS 文件的写入功能。
//! 对应 rds2cpp 的写入相关文件。

// 工具函数
pub mod shared_info;
pub mod utils;

// 写入入口
pub mod object;
pub mod rds;

// 基础类型写入
pub mod atomic;
pub mod attributes;
pub mod single_string;

// 复合类型写入
pub mod builtin;
pub mod expression;
pub mod language;
pub mod list;
pub mod pairlist;
pub mod s4;

// 重新导出公共写入函数
pub use atomic::{
    write_complex_body, write_double_body, write_integer_body, write_logical_body, write_raw_body,
    write_string_body,
};
pub use attributes::write_attributes;
pub use builtin::write_builtin_body;
pub use expression::write_expression_body;
pub use language::write_language;
pub use list::write_list_body;
pub use object::write_object;
pub use pairlist::write_pairlist;
pub use rds::{write_rds, write_rds_file};
pub use s4::write_s4;
pub use shared_info::SharedWriteInfo;
pub use single_string::write_single_string;
