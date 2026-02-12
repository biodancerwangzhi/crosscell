//! 列表写入
//!
//! 写入 R 的通用列表（VEC 类型）。

use std::io::Write;
use crate::rds::error::Result;
use crate::rds::r_object::{GenericVector, RObject};
use super::utils::write_length;
use super::shared_info::SharedWriteInfo;

/// 写入列表体（不含头部）
pub fn write_list_body<W, F>(
    vec: &GenericVector,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
    write_obj: &mut F,
) -> Result<()>
where
    W: Write,
    F: FnMut(&RObject, &mut W, &mut SharedWriteInfo) -> Result<()>,
{
    write_length(vec.data.len(), writer)?;
    for item in &vec.data {
        write_obj(item, writer, shared)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::symbol::Symbol;
    use crate::rds::environment::Environment;
    use crate::rds::external_pointer::ExternalPointer;
    use super::super::utils::write_header;

    fn mk<'a>(s: &'a [Symbol], e: &'a [Environment], p: &'a [ExternalPointer]) -> SharedWriteInfo<'a> {
        SharedWriteInfo::new(s, e, p)
    }
    fn nw(_: &RObject, w: &mut Cursor<Vec<u8>>, _: &mut SharedWriteInfo) -> Result<()> {
        write_header(SEXPType::NilValue, w)
    }

    #[test]
    fn test_empty_list() {
        let v = GenericVector { data: vec![], attributes: Default::default() };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_list_body(&v, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        // Just length(4)
        assert_eq!(data.len(), 4);
    }

    #[test]
    fn test_list_with_items() {
        let v = GenericVector {
            data: vec![RObject::Null, RObject::Null],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_list_body(&v, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        // length(4) + 2 * NilValue(4) = 12
        assert_eq!(data.len(), 12);
    }
}
