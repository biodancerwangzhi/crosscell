//! 表达式向量写入
//!
//! 写入 R 的表达式向量。

use super::shared_info::SharedWriteInfo;
use super::utils::write_length;
use crate::rds::error::Result;
use crate::rds::r_object::{ExpressionVector, RObject};
use std::io::Write;

/// 写入表达式向量体（不含头部）
pub fn write_expression_body<W, F>(
    vec: &ExpressionVector,
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
    use super::super::utils::write_header;
    use super::*;
    use crate::rds::environment::Environment;
    use crate::rds::external_pointer::ExternalPointer;
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::symbol::Symbol;
    use std::io::Cursor;

    fn mk<'a>(
        s: &'a [Symbol],
        e: &'a [Environment],
        p: &'a [ExternalPointer],
    ) -> SharedWriteInfo<'a> {
        SharedWriteInfo::new(s, e, p)
    }
    fn nw(_: &RObject, w: &mut Cursor<Vec<u8>>, _: &mut SharedWriteInfo) -> Result<()> {
        write_header(SEXPType::NilValue, w)
    }

    #[test]
    fn test_empty_expression() {
        let v = ExpressionVector {
            data: vec![],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_expression_body(&v, &mut buf, &mut sh, &mut nw).unwrap();
        assert_eq!(buf.into_inner().len(), 4); // just length
    }

    #[test]
    fn test_expression_with_items() {
        let v = ExpressionVector {
            data: vec![RObject::Null, RObject::Null],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_expression_body(&v, &mut buf, &mut sh, &mut nw).unwrap();
        assert_eq!(buf.into_inner().len(), 12); // length(4) + 2*NilValue(4)
    }
}
