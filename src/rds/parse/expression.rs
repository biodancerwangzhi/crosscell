//! 表达式解析
use std::io::Read;
use crate::rds::error::Result;
use crate::rds::r_object::{ExpressionVector, RObject};
use super::utils::read_length;
use super::shared_info::SharedParseInfo;

/// 解析表达式向量体
pub fn parse_expression_body<R: Read, F>(
    reader: &mut R, shared: &mut SharedParseInfo, parse_fn: &mut F,
) -> Result<ExpressionVector>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    for _ in 0..length {
        data.push(parse_fn(reader, shared)?);
    }
    Ok(ExpressionVector { data, attributes: Default::default() })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_empty_expression() {
        let d = 0i32.to_be_bytes().to_vec();
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |_: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| Ok(RObject::Null);
        let expr = parse_expression_body(&mut c, &mut s, &mut pf).unwrap();
        assert!(expr.data.is_empty());
    }
}
