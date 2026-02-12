//! 头部解析
//!
//! 解析 RDS 对象的 4 字节头部。
//! 对应 rds2cpp 的头部解析逻辑。

use std::io::Read;
use crate::rds::error::Result;
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use super::utils::{Header, read_header as read_header_bytes, get_sexp_type, get_flags, get_gp_high};

/// 头部标志位常量
pub mod flags {
    /// 对象标志位 (bit 0)
    pub const IS_OBJECT: u8 = 0x01;
    /// 属性标志位 (bit 1)
    pub const HAS_ATTRIBUTES: u8 = 0x02;
    /// 标签标志位 (bit 2)
    pub const HAS_TAG: u8 = 0x04;
}

/// 编码标志位常量（在 gp 字段的高字节）
pub mod encoding_flags {
    /// Latin1 编码标志
    pub const LATIN1: u8 = 0x04;
    /// UTF-8 编码标志
    pub const UTF8: u8 = 0x08;
    /// ASCII 编码标志
    pub const ASCII: u8 = 0x40;
}

/// 解析后的头部信息
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParsedHeader {
    /// 原始 4 字节头部
    pub raw: Header,
    /// SEXP 类型
    pub sexp_type: SEXPType,
    /// 是否有属性
    pub has_attributes: bool,
    /// 是否有标签
    pub has_tag: bool,
    /// 是否是对象
    pub is_object: bool,
}

impl ParsedHeader {
    /// 从原始头部字节创建解析后的头部
    pub fn from_raw(raw: Header) -> Option<Self> {
        let sexp_type = SEXPType::from_u8(get_sexp_type(&raw))?;
        let flags = get_flags(&raw);
        
        Some(Self {
            raw,
            sexp_type,
            has_attributes: (flags & flags::HAS_ATTRIBUTES) != 0,
            has_tag: (flags & flags::HAS_TAG) != 0,
            is_object: (flags & flags::IS_OBJECT) != 0,
        })
    }

    /// 获取字符串编码（用于 CHAR 类型）
    pub fn get_encoding(&self) -> StringEncoding {
        let gp = get_gp_high(&self.raw);
        
        if (gp & encoding_flags::UTF8) != 0 {
            StringEncoding::Utf8
        } else if (gp & encoding_flags::LATIN1) != 0 {
            StringEncoding::Latin1
        } else if (gp & encoding_flags::ASCII) != 0 {
            StringEncoding::Ascii
        } else {
            StringEncoding::None
        }
    }
}


/// 从 reader 解析头部
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的头部信息
///
/// # 错误
/// 如果读取失败或类型无效返回错误
pub fn parse_header<R: Read>(reader: &mut R) -> Result<ParsedHeader> {
    let raw = read_header_bytes(reader)?;
    ParsedHeader::from_raw(raw).ok_or_else(|| {
        crate::rds::error::RdsError::UnsupportedType(get_sexp_type(&raw))
    })
}

/// 检查头部是否有属性标记
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// 如果有属性标记返回 true
#[inline]
pub fn has_attributes(header: &Header) -> bool {
    (get_flags(header) & flags::HAS_ATTRIBUTES) != 0
}

/// 检查头部是否有标签标记
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// 如果有标签标记返回 true
#[inline]
pub fn has_tag(header: &Header) -> bool {
    (get_flags(header) & flags::HAS_TAG) != 0
}

/// 检查头部是否是对象
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// 如果是对象返回 true
#[inline]
pub fn is_object(header: &Header) -> bool {
    (get_flags(header) & flags::IS_OBJECT) != 0
}

/// 从头部获取引用索引
///
/// 用于 REF 类型，引用索引存储在头部的前 3 字节（大端序）
/// 参考 rds2cpp SharedParseInfo.hpp compute_reference_index
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// 引用索引
pub fn get_reference_index(header: &Header) -> usize {
    let b0 = header[0] as usize;
    let b1 = header[1] as usize;
    let b2 = header[2] as usize;
    // 大端序：header[0] 是高字节，header[2] 是低字节
    (b0 << 16) | (b1 << 8) | b2
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parsed_header_from_raw() {
        // 测试整数类型头部
        let raw: Header = [0x00, 0x00, 0x00, 0x0D]; // Int type = 13
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.sexp_type, SEXPType::Int);
        assert!(!parsed.has_attributes);
        assert!(!parsed.has_tag);
        assert!(!parsed.is_object);
    }

    #[test]
    fn test_parsed_header_with_attributes() {
        // 测试带属性的头部
        let raw: Header = [0x00, 0x00, 0x02, 0x0D]; // Int with has_attributes
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.sexp_type, SEXPType::Int);
        assert!(parsed.has_attributes);
        assert!(!parsed.has_tag);
    }

    #[test]
    fn test_parsed_header_with_tag() {
        // 测试带标签的头部
        let raw: Header = [0x00, 0x00, 0x04, 0x02]; // List with has_tag
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.sexp_type, SEXPType::List);
        assert!(parsed.has_tag);
    }

    #[test]
    fn test_parsed_header_is_object() {
        // 测试对象标志
        let raw: Header = [0x00, 0x00, 0x01, 0x19]; // S4 with is_object
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.sexp_type, SEXPType::S4);
        assert!(parsed.is_object);
    }

    #[test]
    fn test_parsed_header_invalid_type() {
        // 测试无效类型
        let raw: Header = [0x00, 0x00, 0x00, 0x0B]; // 11 is invalid
        let parsed = ParsedHeader::from_raw(raw);
        assert!(parsed.is_none());
    }

    #[test]
    fn test_get_encoding_utf8() {
        let raw: Header = [0x00, 0x08, 0x00, 0x09]; // CHAR with UTF-8
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.get_encoding(), StringEncoding::Utf8);
    }

    #[test]
    fn test_get_encoding_latin1() {
        let raw: Header = [0x00, 0x04, 0x00, 0x09]; // CHAR with Latin1
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.get_encoding(), StringEncoding::Latin1);
    }

    #[test]
    fn test_get_encoding_ascii() {
        let raw: Header = [0x00, 0x40, 0x00, 0x09]; // CHAR with ASCII
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.get_encoding(), StringEncoding::Ascii);
    }

    #[test]
    fn test_get_encoding_none() {
        let raw: Header = [0x00, 0x00, 0x00, 0x09]; // CHAR with no encoding
        let parsed = ParsedHeader::from_raw(raw).unwrap();
        assert_eq!(parsed.get_encoding(), StringEncoding::None);
    }

    #[test]
    fn test_parse_header_from_reader() {
        let data = vec![0x00, 0x00, 0x00, 0x0D]; // Int type
        let mut cursor = Cursor::new(data);
        let parsed = parse_header(&mut cursor).unwrap();
        assert_eq!(parsed.sexp_type, SEXPType::Int);
    }

    #[test]
    fn test_parse_header_invalid_type() {
        let data = vec![0x00, 0x00, 0x00, 0x0B]; // Invalid type 11
        let mut cursor = Cursor::new(data);
        let result = parse_header(&mut cursor);
        assert!(result.is_err());
    }

    #[test]
    fn test_has_attributes() {
        let header: Header = [0x00, 0x00, 0x02, 0x00];
        assert!(has_attributes(&header));
        
        let header: Header = [0x00, 0x00, 0x00, 0x00];
        assert!(!has_attributes(&header));
    }

    #[test]
    fn test_has_tag() {
        let header: Header = [0x00, 0x00, 0x04, 0x00];
        assert!(has_tag(&header));
        
        let header: Header = [0x00, 0x00, 0x00, 0x00];
        assert!(!has_tag(&header));
    }

    #[test]
    fn test_is_object() {
        let header: Header = [0x00, 0x00, 0x01, 0x00];
        assert!(is_object(&header));
        
        let header: Header = [0x00, 0x00, 0x00, 0x00];
        assert!(!is_object(&header));
    }

    #[test]
    fn test_get_reference_index() {
        // 测试引用索引提取（大端序：header[0] 是高字节）
        let header: Header = [0x00, 0x00, 0x01, 0xFF]; // index = 1
        assert_eq!(get_reference_index(&header), 1);
        
        let header: Header = [0x00, 0x01, 0x00, 0xFF]; // index = 256
        assert_eq!(get_reference_index(&header), 256);
        
        let header: Header = [0x01, 0x00, 0x00, 0xFF]; // index = 65536
        assert_eq!(get_reference_index(&header), 65536);
        
        let header: Header = [0x56, 0x34, 0x12, 0xFF]; // index = 0x563412
        assert_eq!(get_reference_index(&header), 0x563412);
    }
}
