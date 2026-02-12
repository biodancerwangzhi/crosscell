//! 配对列表写入
//!
//! 写入 R 的配对列表。每个元素是一个 (tag, value) 对，
//! 以 NILVALUE 终止。

use std::io::Write;
use crate::rds::error::Result;
use crate::rds::r_object::{PairList, RObject};
use crate::rds::sexp_type::SEXPType;
use super::utils::{write_header, write_pairlist_header};
use super::shared_info::SharedWriteInfo;
use super::attributes::write_attributes;

/// 写入配对列表
///
/// # 参数
/// * `pairlist` - 配对列表数据
/// * `writer` - 写入目标
/// * `shared` - 共享写入信息
/// * `write_obj` - 写入对象的闭包
pub fn write_pairlist<W, F>(
    pairlist: &PairList,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
    write_obj: &mut F,
) -> Result<()>
where
    W: Write,
    F: FnMut(&RObject, &mut W, &mut SharedWriteInfo) -> Result<()>,
{
    for i in 0..pairlist.data.len() {
        let has_tag = pairlist.has_tag.get(i).copied().unwrap_or(false);
        let is_first = i == 0;

        // 写入 LIST 头部
        if is_first && !pairlist.attributes.is_empty() {
            // 第一个节点带属性标记
            let flags = 0x02 | if has_tag { 0x04 } else { 0x00 };
            let header = [0, 0, flags, SEXPType::List as u8];
            writer.write_all(&header)?;
            // 写入属性
            write_attributes(&pairlist.attributes, writer, shared, write_obj)?;
        } else {
            write_pairlist_header(has_tag, writer)?;
        }

        // 写入标签（符号）
        if has_tag {
            let name = &pairlist.tag_names[i];
            let enc = pairlist.tag_encodings[i];
            shared.write_symbol(name, enc, writer)?;
        }

        // 写入值
        write_obj(&pairlist.data[i], writer, shared)?;
    }

    // 终止
    write_header(SEXPType::NilValue, writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::string_encoding::StringEncoding;
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
    fn test_empty_pairlist() {
        let pl = PairList::default();
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_pairlist(&pl, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        // Only NilValue terminator
        assert_eq!(data.len(), 4);
        assert_eq!(data[3], SEXPType::NilValue as u8);
    }

    #[test]
    fn test_tagged_pairlist() {
        let pl = PairList {
            data: vec![RObject::Null],
            has_tag: vec![true],
            tag_names: vec!["x".into()],
            tag_encodings: vec![StringEncoding::Utf8],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_pairlist(&pl, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        // LIST header + symbol + value + NilValue
        assert!(!data.is_empty());
        assert_eq!(data[3], SEXPType::List as u8);
    }

    #[test]
    fn test_untagged_pairlist() {
        let pl = PairList {
            data: vec![RObject::Null, RObject::Null],
            has_tag: vec![false, false],
            tag_names: vec!["".into(), "".into()],
            tag_encodings: vec![StringEncoding::Utf8, StringEncoding::Utf8],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_pairlist(&pl, &mut buf, &mut sh, &mut nw).unwrap();
        // reference_count 从 1 开始，没有符号写入
        assert_eq!(sh.reference_count(), 1);
    }
}
