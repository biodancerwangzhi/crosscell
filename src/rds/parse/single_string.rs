//! 单字符串解析
//!
//! 解析 R 的单个字符串（CHARSXP）。

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use super::utils::{read_header, get_sexp_type, quick_extract};
use super::header::{ParsedHeader, encoding_flags};

/// 解析后的字符串信息
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringInfo {
    pub value: String,
    pub encoding: StringEncoding,
    pub missing: bool,
}

impl Default for StringInfo {
    fn default() -> Self {
        Self { value: String::new(), encoding: StringEncoding::Utf8, missing: false }
    }
}

impl StringInfo {
    pub fn na() -> Self {
        Self { value: String::new(), encoding: StringEncoding::None, missing: true }
    }

    pub fn new(value: String, encoding: StringEncoding) -> Self {
        Self { value, encoding, missing: false }
    }
}

fn get_string_encoding(header: &[u8; 4]) -> StringEncoding {
    let gp = header[1];
    if (gp & encoding_flags::UTF8) != 0 { StringEncoding::Utf8 }
    else if (gp & encoding_flags::LATIN1) != 0 { StringEncoding::Latin1 }
    else if (gp & encoding_flags::ASCII) != 0 { StringEncoding::Ascii }
    else { StringEncoding::Utf8 }
}

/// 解析单个字符串（CHARSXP）
pub fn parse_single_string<R: Read>(reader: &mut R) -> Result<StringInfo> {
    let header = read_header(reader)?;
    let sexp_type = get_sexp_type(&header);
    if sexp_type != SEXPType::Char as u8 {
        return Err(RdsError::ParseError {
            context: "single_string".to_string(),
            message: format!("Expected CHAR type (9), got {}", sexp_type),
        });
    }
    let encoding = get_string_encoding(&header);
    
    // 读取长度（直接读取 i32 以检测 NA）
    let mut len_buf = [0u8; 4];
    reader.read_exact(&mut len_buf)?;
    let length_i32 = i32::from_be_bytes(len_buf);
    
    // -1 表示 NA 字符串
    if length_i32 == -1 {
        return Ok(StringInfo::na());
    }
    
    let length = length_i32 as usize;
    let bytes = quick_extract(reader, length)?;
    let value = String::from_utf8_lossy(&bytes).into_owned();
    Ok(StringInfo::new(value, encoding))
}

/// 解析单个字符串（带已解析的头部）
pub fn parse_single_string_with_header<R: Read>(
    reader: &mut R, header: &ParsedHeader,
) -> Result<StringInfo> {
    if header.sexp_type != SEXPType::Char {
        return Err(RdsError::ParseError {
            context: "single_string".to_string(),
            message: format!("Expected CHAR type, got {:?}", header.sexp_type),
        });
    }
    let encoding = header.get_encoding();
    
    // 读取长度（直接读取 i32 以检测 NA）
    let mut len_buf = [0u8; 4];
    reader.read_exact(&mut len_buf)?;
    let length_i32 = i32::from_be_bytes(len_buf);
    
    // -1 表示 NA 字符串
    if length_i32 == -1 {
        return Ok(StringInfo::na());
    }
    
    let length = length_i32 as usize;
    let bytes = quick_extract(reader, length)?;
    let value = String::from_utf8_lossy(&bytes).into_owned();
    Ok(StringInfo::new(value, encoding))
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_char_header(enc: StringEncoding) -> [u8; 4] {
        let gp = match enc {
            StringEncoding::Utf8 => encoding_flags::UTF8,
            StringEncoding::Latin1 => encoding_flags::LATIN1,
            StringEncoding::Ascii => encoding_flags::ASCII,
            StringEncoding::None => 0,
        };
        [0x00, gp, 0x00, SEXPType::Char as u8]
    }

    fn make_length(len: i32) -> [u8; 4] { len.to_be_bytes() }

    #[test]
    fn test_string_info_default() {
        let info = StringInfo::default();
        assert_eq!(info.value, "");
        assert!(!info.missing);
    }

    #[test]
    fn test_string_info_na() {
        let info = StringInfo::na();
        assert!(info.missing);
    }

    #[test]
    fn test_parse_simple_string() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(5));
        data.extend_from_slice(b"hello");
        let mut cursor = Cursor::new(data);
        let info = parse_single_string(&mut cursor).unwrap();
        assert_eq!(info.value, "hello");
        assert!(!info.missing);
    }

    #[test]
    fn test_parse_empty_string() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(0));
        let mut cursor = Cursor::new(data);
        let info = parse_single_string(&mut cursor).unwrap();
        assert_eq!(info.value, "");
        assert!(!info.missing);
    }

    #[test]
    fn test_parse_na_string() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::None));
        data.extend_from_slice(&make_length(-1));
        let mut cursor = Cursor::new(data);
        let info = parse_single_string(&mut cursor).unwrap();
        assert!(info.missing);
    }

    #[test]
    fn test_parse_unicode_string() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        let s = "你好";
        data.extend_from_slice(&make_length(s.len() as i32));
        data.extend_from_slice(s.as_bytes());
        let mut cursor = Cursor::new(data);
        let info = parse_single_string(&mut cursor).unwrap();
        assert_eq!(info.value, "你好");
    }

    #[test]
    fn test_parse_wrong_type_error() {
        let mut data = Vec::new();
        data.extend_from_slice(&[0x00, 0x00, 0x00, SEXPType::Int as u8]);
        data.extend_from_slice(&make_length(5));
        let mut cursor = Cursor::new(data);
        assert!(parse_single_string(&mut cursor).is_err());
    }

    #[test]
    fn test_get_string_encoding_variants() {
        assert_eq!(get_string_encoding(&[0, encoding_flags::UTF8, 0, 9]), StringEncoding::Utf8);
        assert_eq!(get_string_encoding(&[0, encoding_flags::LATIN1, 0, 9]), StringEncoding::Latin1);
        assert_eq!(get_string_encoding(&[0, encoding_flags::ASCII, 0, 9]), StringEncoding::Ascii);
        assert_eq!(get_string_encoding(&[0, 0, 0, 9]), StringEncoding::Utf8);
    }
}
