//! 内置函数写入
//!
//! 写入 R 的内置函数引用。

use std::io::Write;
use crate::rds::error::Result;
use crate::rds::r_object::BuiltInFunction;
use super::utils::{write_length, write_bytes};

/// 写入内置函数体（不含头部）
pub fn write_builtin_body<W: Write>(func: &BuiltInFunction, writer: &mut W) -> Result<()> {
    write_length(func.name.len(), writer)?;
    write_bytes(func.name.as_bytes(), writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::parse::builtin::parse_builtin_body;

    #[test]
    fn test_roundtrip() {
        let f = BuiltInFunction { name: "sum".into() };
        let mut buf = Cursor::new(Vec::new());
        write_builtin_body(&f, &mut buf).unwrap();
        let mut r = Cursor::new(buf.into_inner());
        let parsed = parse_builtin_body(&mut r).unwrap();
        assert_eq!(parsed.name, "sum");
    }

    #[test]
    fn test_empty_name() {
        let f = BuiltInFunction { name: String::new() };
        let mut buf = Cursor::new(Vec::new());
        write_builtin_body(&f, &mut buf).unwrap();
        let mut r = Cursor::new(buf.into_inner());
        let parsed = parse_builtin_body(&mut r).unwrap();
        assert_eq!(parsed.name, "");
    }
}
