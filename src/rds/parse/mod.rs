//! RDS 解析模块
//!
//! 提供 RDS 文件的解析功能。
//! 对应 rds2cpp 的解析相关文件。

// 工具函数
pub mod header;
pub mod shared_info;
pub mod utils;

// 解析入口
pub mod object;
pub mod rds;

// 基础类型解析
pub mod atomic;
pub mod single_string;

// 复合类型解析
pub mod altrep;
pub mod attributes;
pub mod builtin;
pub mod environment;
pub mod expression;
pub mod external_pointer;
pub mod language;
pub mod list;
pub mod pairlist;
pub mod s4;
pub mod symbol;

// 重新导出公共解析函数
pub use utils::{
    as_char, get_flags, get_gp_high, get_gp_low, get_sexp_type, is_little_endian, quick_extract,
    read_f64, read_f64_be, read_header, read_i32, read_i32_be, read_length, Header,
};

pub use header::{
    get_reference_index, has_attributes, has_tag, is_object, parse_header, ParsedHeader,
};

pub use shared_info::{ReferenceType, SharedParseInfo, StoredString};

pub use single_string::{parse_single_string, parse_single_string_with_header, StringInfo};

pub use atomic::{
    parse_complex_body, parse_double_body, parse_integer_body, parse_logical_body, parse_raw_body,
    parse_string_body, parse_string_body_with_ref,
};

pub use symbol::parse_symbol_body;

pub use attributes::{parse_attributes, parse_attributes_body};

pub use pairlist::parse_pairlist_body;

pub use list::parse_list_body;

pub use builtin::parse_builtin_body;

pub use expression::parse_expression_body;

pub use language::parse_language_body;

pub use environment::{
    parse_base_environment_body, parse_empty_environment_body, parse_global_environment_body,
    parse_new_environment_body,
};

pub use external_pointer::parse_external_pointer_body;

pub use s4::parse_s4_body;

pub use altrep::parse_altrep_body;

pub use object::parse_object;

pub use rds::{parse_rds, parse_rds_file, ParseRdsOptions};
