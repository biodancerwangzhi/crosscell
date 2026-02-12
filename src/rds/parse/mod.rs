//! RDS 解析模块
//!
//! 提供 RDS 文件的解析功能。
//! 对应 rds2cpp 的解析相关文件。

// 工具函数
pub mod utils;
pub mod shared_info;
pub mod header;

// 解析入口
pub mod rds;
pub mod object;

// 基础类型解析
pub mod single_string;
pub mod atomic;

// 复合类型解析
pub mod symbol;
pub mod attributes;
pub mod pairlist;
pub mod list;
pub mod builtin;
pub mod expression;
pub mod language;
pub mod environment;
pub mod external_pointer;
pub mod s4;
pub mod altrep;

// 重新导出公共解析函数
pub use utils::{
    Header, is_little_endian, read_i32_be, read_f64_be, read_i32, read_f64,
    read_length, quick_extract, as_char, read_header, get_sexp_type, get_flags,
    get_gp_high, get_gp_low,
};

pub use header::{
    ParsedHeader, parse_header, has_attributes, has_tag, is_object, get_reference_index,
};

pub use shared_info::{SharedParseInfo, ReferenceType, StoredString};

pub use single_string::{StringInfo, parse_single_string, parse_single_string_with_header};

pub use atomic::{
    parse_integer_body, parse_logical_body, parse_double_body,
    parse_raw_body, parse_complex_body, parse_string_body, parse_string_body_with_ref,
};

pub use symbol::parse_symbol_body;

pub use attributes::{parse_attributes, parse_attributes_body};

pub use pairlist::parse_pairlist_body;

pub use list::parse_list_body;

pub use builtin::parse_builtin_body;

pub use expression::parse_expression_body;

pub use language::parse_language_body;

pub use environment::{
    parse_new_environment_body, parse_global_environment_body,
    parse_base_environment_body, parse_empty_environment_body,
};

pub use external_pointer::parse_external_pointer_body;

pub use s4::parse_s4_body;

pub use altrep::parse_altrep_body;

pub use object::parse_object;

pub use rds::{parse_rds, parse_rds_file, ParseRdsOptions};
