//! 内置函数解析
use std::io::Read;
use crate::rds::error::Result;
use crate::rds::r_object::BuiltInFunction;
use super::utils::{read_length, quick_extract};

/// 解析内置函数体
pub fn parse_builtin_body<R: Read>(reader: &mut R) -> Result<BuiltInFunction> {
    let length = read_length(reader)?;
    let bytes = quick_extract(reader, length)?;
    let name = String::from_utf8_lossy(&bytes).into_owned();
    Ok(BuiltInFunction { name })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_builtin() {
        let mut d = Vec::new();
        d.extend_from_slice(&3i32.to_be_bytes());
        d.extend_from_slice(b"sum");
        let mut c = Cursor::new(d);
        let f = parse_builtin_body(&mut c).unwrap();
        assert_eq!(f.name, "sum");
    }

    #[test]
    fn test_parse_builtin_empty() {
        let d = 0i32.to_be_bytes().to_vec();
        let mut c = Cursor::new(d);
        let f = parse_builtin_body(&mut c).unwrap();
        assert_eq!(f.name, "");
    }
}
