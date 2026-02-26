//! 外部指针解析
use super::attributes::parse_attributes;
use super::header::{has_attributes, ParsedHeader};
use super::shared_info::SharedParseInfo;
use crate::rds::error::Result;
use crate::rds::external_pointer::ExternalPointer;
use crate::rds::r_object::{ExternalPointerIndex, RObject};
use std::io::Read;

/// 解析外部指针体
pub fn parse_external_pointer_body<R: Read, F>(
    reader: &mut R,
    header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<ExternalPointerIndex>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    let index = shared.request_external_pointer();
    let protection = Box::new(parse_fn(reader, shared)?);
    let tag = Box::new(parse_fn(reader, shared)?);
    let ptr = ExternalPointer::with_values(*protection, *tag);
    if has_attributes(&header.raw) {
        let _attrs = parse_attributes(reader, shared, parse_fn)?;
    }
    shared.update_external_pointer(index, ptr);
    Ok(ExternalPointerIndex { index })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::parse::utils::read_header;
    use crate::rds::sexp_type::SEXPType;
    use std::io::Cursor;

    #[test]
    fn test_external_pointer() {
        let mut d = Vec::new();
        d.extend_from_slice(&[0, 0, 0, SEXPType::NilValue as u8]); // protection
        d.extend_from_slice(&[0, 0, 0, SEXPType::NilValue as u8]); // tag
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |r: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| {
            let _ = read_header(r)?;
            Ok(RObject::Null)
        };
        let hdr = ParsedHeader::from_raw([0, 0, 0, SEXPType::Extptr as u8]).unwrap();
        let idx = parse_external_pointer_body(&mut c, &hdr, &mut s, &mut pf).unwrap();
        assert_eq!(idx.index, 0);
    }
}
