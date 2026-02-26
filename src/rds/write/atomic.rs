//! 原子向量写入
//!
//! 写入所有原子向量类型（整数、逻辑、浮点、原始、复数、字符串）。
//! 对应 rds2cpp 的原子向量写入逻辑。

use super::single_string::write_single_string;
use super::utils::{write_bytes, write_f64, write_i32, write_length};
use crate::rds::error::Result;
use crate::rds::r_object::{
    ComplexVector, DoubleVector, IntegerVector, LogicalVector, RawVector, StringVector,
};
use std::io::Write;

/// 写入整数向量体（不含头部）
pub fn write_integer_body<W: Write>(vec: &IntegerVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    for &val in &vec.data {
        write_i32(val, writer)?;
    }
    Ok(())
}

/// 写入逻辑向量体（不含头部）
pub fn write_logical_body<W: Write>(vec: &LogicalVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    for &val in &vec.data {
        write_i32(val, writer)?;
    }
    Ok(())
}

/// 写入浮点向量体（不含头部）
pub fn write_double_body<W: Write>(vec: &DoubleVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    for &val in &vec.data {
        write_f64(val, writer)?;
    }
    Ok(())
}

/// 写入原始向量体（不含头部）
pub fn write_raw_body<W: Write>(vec: &RawVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    write_bytes(&vec.data, writer)?;
    Ok(())
}

/// 写入复数向量体（不含头部）
pub fn write_complex_body<W: Write>(vec: &ComplexVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    for c in &vec.data {
        write_f64(c.re, writer)?;
        write_f64(c.im, writer)?;
    }
    Ok(())
}

/// 写入字符串向量体（不含头部）
pub fn write_string_body<W: Write>(vec: &StringVector, writer: &mut W) -> Result<()> {
    write_length(vec.data.len(), writer)?;
    for i in 0..vec.data.len() {
        let missing = vec.missing.get(i).copied().unwrap_or(false);
        let encoding = vec.encodings.get(i).copied().unwrap_or_default();
        write_single_string(&vec.data[i], encoding, missing, writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::parse::atomic;
    use crate::rds::string_encoding::StringEncoding;
    use num_complex::Complex64;
    use std::io::Cursor;

    fn rt_int(d: Vec<i32>) -> Vec<i32> {
        let v = IntegerVector {
            data: d,
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_integer_body(&v, &mut b).unwrap();
        atomic::parse_integer_body(&mut Cursor::new(b.into_inner()))
            .unwrap()
            .data
    }

    #[test]
    fn test_integer_roundtrip() {
        assert_eq!(rt_int(vec![]), Vec::<i32>::new());
        assert_eq!(
            rt_int(vec![1, -1, 0, i32::MAX, i32::MIN]),
            vec![1, -1, 0, i32::MAX, i32::MIN]
        );
    }

    #[test]
    fn test_logical_roundtrip() {
        let v = LogicalVector {
            data: vec![0, 1, i32::MIN],
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_logical_body(&v, &mut b).unwrap();
        let p = atomic::parse_logical_body(&mut Cursor::new(b.into_inner())).unwrap();
        assert_eq!(p.data, v.data);
    }

    #[test]
    fn test_double_roundtrip() {
        let v = DoubleVector {
            data: vec![0.0, 1.5, -3.14, f64::INFINITY],
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_double_body(&v, &mut b).unwrap();
        let p = atomic::parse_double_body(&mut Cursor::new(b.into_inner())).unwrap();
        assert_eq!(p.data, v.data);
    }

    #[test]
    fn test_raw_roundtrip() {
        let v = RawVector {
            data: vec![0x00, 0xFF, 0x42],
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_raw_body(&v, &mut b).unwrap();
        let p = atomic::parse_raw_body(&mut Cursor::new(b.into_inner())).unwrap();
        assert_eq!(p.data, v.data);
    }

    #[test]
    fn test_complex_roundtrip() {
        let v = ComplexVector {
            data: vec![Complex64::new(1.0, 2.0), Complex64::new(-3.0, 0.0)],
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_complex_body(&v, &mut b).unwrap();
        let p = atomic::parse_complex_body(&mut Cursor::new(b.into_inner())).unwrap();
        assert_eq!(p.data, v.data);
    }

    #[test]
    fn test_string_roundtrip() {
        let v = StringVector {
            data: vec!["hello".into(), "".into(), "world".into()],
            encodings: vec![
                StringEncoding::Utf8,
                StringEncoding::Utf8,
                StringEncoding::Latin1,
            ],
            missing: vec![false, true, false],
            attributes: Default::default(),
        };
        let mut b = Cursor::new(Vec::new());
        write_string_body(&v, &mut b).unwrap();
        let p = atomic::parse_string_body(&mut Cursor::new(b.into_inner())).unwrap();
        assert_eq!(p.data[0], "hello");
        assert!(p.missing[1]);
        assert_eq!(p.data[2], "world");
        // 注意：解析器可能会将 Latin1 规范化为 UTF8（对于纯 ASCII 内容）
        // 这是正确的行为，因为 ASCII 内容在两种编码中是相同的
        assert!(
            p.encodings[2] == StringEncoding::Latin1 || p.encodings[2] == StringEncoding::Utf8,
            "Encoding should be Latin1 or UTF8, got {:?}",
            p.encodings[2]
        );
    }
}
