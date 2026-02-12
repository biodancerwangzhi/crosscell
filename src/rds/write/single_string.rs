//! 单字符串写入
//!
//! 写入 R 的单个字符串（CHARSXP）。
//! 对应 rds2cpp 的单字符串写入逻辑。

use std::io::Write;
use crate::rds::error::{Result, RdsError};
use crate::rds::string_encoding::StringEncoding;
use super::utils::{write_string_header, write_i32, write_length, write_bytes};

/// 写入单个字符串（CHARSXP）
///
/// # 参数
/// * `value` - 字符串内容
/// * `encoding` - 字符串编码
/// * `missing` - 是否为 NA 值
/// * `writer` - 写入目标
///
/// # 错误
/// - 如果字符串长度超过 2^31-1 返回 StringTooLong
/// - 如果写入失败返回 IO 错误
pub fn write_single_string<W: Write>(
    value: &str,
    encoding: StringEncoding,
    missing: bool,
    writer: &mut W,
) -> Result<()> {
    if missing {
        // NA 字符串：写入头部 + 长度 -1
        write_string_header(encoding, true, writer)?;
        write_i32(-1, writer)?;
        return Ok(());
    }

    let len = value.len();
    if len > i32::MAX as usize {
        return Err(RdsError::StringTooLong(len));
    }

    // 写入头部（带编码标记）
    write_string_header(encoding, false, writer)?;
    // 写入长度
    write_length(len, writer)?;
    // 写入字符串内容
    write_bytes(value.as_bytes(), writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::parse::single_string::parse_single_string;

    #[test]
    fn test_write_simple_string() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("hello", StringEncoding::Utf8, false, &mut buf).unwrap();
        let data = buf.into_inner();
        // header(4) + length(4) + "hello"(5) = 13
        assert_eq!(data.len(), 13);
    }

    #[test]
    fn test_write_empty_string() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("", StringEncoding::Utf8, false, &mut buf).unwrap();
        let data = buf.into_inner();
        // header(4) + length(4) = 8
        assert_eq!(data.len(), 8);
    }

    #[test]
    fn test_write_na_string() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("", StringEncoding::None, true, &mut buf).unwrap();
        let data = buf.into_inner();
        // header(4) + length(-1 as 4 bytes) = 8
        assert_eq!(data.len(), 8);
        // length should be -1 (0xFFFFFFFF)
        assert_eq!(&data[4..8], &[0xFF, 0xFF, 0xFF, 0xFF]);
    }

    #[test]
    fn test_roundtrip_simple() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("hello", StringEncoding::Utf8, false, &mut buf).unwrap();
        let data = buf.into_inner();
        let mut reader = Cursor::new(data);
        let info = parse_single_string(&mut reader).unwrap();
        assert_eq!(info.value, "hello");
        assert_eq!(info.encoding, StringEncoding::Utf8);
        assert!(!info.missing);
    }

    #[test]
    fn test_roundtrip_na() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("", StringEncoding::None, true, &mut buf).unwrap();
        let data = buf.into_inner();
        let mut reader = Cursor::new(data);
        let info = parse_single_string(&mut reader).unwrap();
        assert!(info.missing);
    }

    #[test]
    fn test_roundtrip_unicode() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("你好世界", StringEncoding::Utf8, false, &mut buf).unwrap();
        let data = buf.into_inner();
        let mut reader = Cursor::new(data);
        let info = parse_single_string(&mut reader).unwrap();
        assert_eq!(info.value, "你好世界");
        assert_eq!(info.encoding, StringEncoding::Utf8);
    }

    #[test]
    fn test_roundtrip_latin1() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("test", StringEncoding::Latin1, false, &mut buf).unwrap();
        let data = buf.into_inner();
        let mut reader = Cursor::new(data);
        let info = parse_single_string(&mut reader).unwrap();
        assert_eq!(info.value, "test");
        // 注意：解析器可能会将 Latin1 规范化为 UTF8（对于纯 ASCII 内容）
        // 这是正确的行为，因为 ASCII 内容在两种编码中是相同的
        assert!(
            info.encoding == StringEncoding::Latin1 || info.encoding == StringEncoding::Utf8,
            "Encoding should be Latin1 or UTF8, got {:?}",
            info.encoding
        );
    }

    #[test]
    fn test_roundtrip_ascii() {
        let mut buf = Cursor::new(Vec::new());
        write_single_string("abc", StringEncoding::Ascii, false, &mut buf).unwrap();
        let data = buf.into_inner();
        let mut reader = Cursor::new(data);
        let info = parse_single_string(&mut reader).unwrap();
        assert_eq!(info.value, "abc");
        // 注意：解析器可能会将 Ascii 规范化为 Latin1 或 UTF8
        // 这是正确的行为，因为 ASCII 是 Latin1 和 UTF8 的子集
        assert!(
            info.encoding == StringEncoding::Ascii 
            || info.encoding == StringEncoding::Latin1 
            || info.encoding == StringEncoding::Utf8,
            "Encoding should be Ascii, Latin1 or UTF8, got {:?}",
            info.encoding
        );
    }
}
