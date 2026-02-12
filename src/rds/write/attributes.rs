//! 属性写入
//!
//! 写入 R 对象的属性（作为配对列表）。
//! 对应 rds2cpp 的属性写入逻辑。

use std::io::Write;
use crate::rds::error::Result;
use crate::rds::r_object::{Attributes, RObject};
use crate::rds::string_encoding::StringEncoding;
use super::utils::{write_header, write_pairlist_header};
use super::shared_info::SharedWriteInfo;
use crate::rds::sexp_type::SEXPType;

/// 写入属性（作为配对列表）
///
/// # 返回
/// `true` 如果写入了属性，`false` 如果属性为空
pub fn write_attributes<W, F>(
    attrs: &Attributes,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
    write_obj: &mut F,
) -> Result<bool>
where
    W: Write,
    F: FnMut(&RObject, &mut W, &mut SharedWriteInfo) -> Result<()>,
{
    if attrs.is_empty() {
        return Ok(false);
    }
    for i in 0..attrs.names.len() {
        write_pairlist_header(true, writer)?;
        let enc = attrs.encodings.get(i).copied().unwrap_or(StringEncoding::Utf8);
        shared.write_symbol(&attrs.names[i], enc, writer)?;
        write_obj(&attrs.values[i], writer, shared)?;
    }
    write_header(SEXPType::NilValue, writer)?;
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::symbol::Symbol;
    use crate::rds::environment::Environment;
    use crate::rds::external_pointer::ExternalPointer;

    fn mk<'a>(s: &'a [Symbol], e: &'a [Environment], p: &'a [ExternalPointer]) -> SharedWriteInfo<'a> {
        SharedWriteInfo::new(s, e, p)
    }
    fn nw(_: &RObject, w: &mut Cursor<Vec<u8>>, _: &mut SharedWriteInfo) -> Result<()> {
        write_header(SEXPType::NilValue, w)
    }

    #[test]
    fn test_empty_returns_false() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        assert!(!write_attributes(&Attributes::new(), &mut buf, &mut sh, &mut nw).unwrap());
    }

    #[test]
    fn test_single_attr() {
        let mut attrs = Attributes::new();
        attrs.add("names".into(), RObject::Null, StringEncoding::Utf8);
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        assert!(write_attributes(&attrs, &mut buf, &mut sh, &mut nw).unwrap());
        let data = buf.into_inner();
        assert_eq!(data[data.len() - 1], SEXPType::NilValue as u8);
    }

    #[test]
    fn test_symbol_dedup() {
        let mut attrs = Attributes::new();
        attrs.add("x".into(), RObject::Null, StringEncoding::Utf8);
        attrs.add("y".into(), RObject::Null, StringEncoding::Utf8);
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_attributes(&attrs, &mut buf, &mut sh, &mut nw).unwrap();
        // reference_count 从 1 开始，两个符号分配 2 和 3
        assert_eq!(sh.reference_count(), 3);
    }
}
