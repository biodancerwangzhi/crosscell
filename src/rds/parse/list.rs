//! 列表解析
use super::shared_info::SharedParseInfo;
use super::utils::read_length;
use crate::rds::error::Result;
use crate::rds::r_object::{GenericVector, RObject};
use std::io::Read;

/// 解析列表体（VEC 类型）
pub fn parse_list_body<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<GenericVector>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    for _ in 0..length {
        data.push(parse_fn(reader, shared)?);
    }
    Ok(GenericVector {
        data,
        attributes: Default::default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::parse::utils::read_header;
    use crate::rds::sexp_type::SEXPType;
    use std::io::Cursor;

    #[test]
    fn test_empty_list() {
        let d = 0i32.to_be_bytes().to_vec();
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |_: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| Ok(RObject::Null);
        let list = parse_list_body(&mut c, &mut s, &mut pf).unwrap();
        assert!(list.data.is_empty());
    }

    #[test]
    fn test_list_with_nulls() {
        let mut d = Vec::new();
        d.extend_from_slice(&2i32.to_be_bytes());
        d.extend_from_slice(&[0, 0, 0, SEXPType::NilValue as u8]);
        d.extend_from_slice(&[0, 0, 0, SEXPType::NilValue as u8]);
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |r: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| {
            let _ = read_header(r)?;
            Ok(RObject::Null)
        };
        let list = parse_list_body(&mut c, &mut s, &mut pf).unwrap();
        assert_eq!(list.data.len(), 2);
    }
}
