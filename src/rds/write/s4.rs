//! S4 对象写入
//!
//! 写入 R 的 S4 对象。S4 对象以配对列表形式写入槽位，
//! 第一个槽位是 class 属性（包含 package 子属性）。

use super::shared_info::SharedWriteInfo;
use super::utils::{write_header, write_pairlist_header};
use crate::rds::error::Result;
use crate::rds::r_object::{Attributes, RObject, S4Object, StringVector};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use std::io::Write;

/// 写入 S4 对象体
///
/// S4 对象格式：
/// 1. class 槽位（字符串向量，带 package 属性）
/// 2. 其他槽位（作为属性配对列表）
/// 3. NILVALUE 终止
pub fn write_s4<W, F>(
    s4: &S4Object,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
    write_obj: &mut F,
) -> Result<()>
where
    W: Write,
    F: FnMut(&RObject, &mut W, &mut SharedWriteInfo) -> Result<()>,
{
    // 写入 class 槽位
    write_pairlist_header(true, writer)?;
    shared.write_symbol("class", StringEncoding::Utf8, writer)?;

    // class 值是一个字符串向量，带 package 属性
    let mut class_attrs = Attributes::new();
    if !s4.package_name.is_empty() {
        let pkg_vec = StringVector {
            data: vec![s4.package_name.clone()],
            encodings: vec![s4.package_encoding],
            missing: vec![false],
            attributes: Default::default(),
        };
        class_attrs.add(
            "package".into(),
            RObject::StringVector(pkg_vec),
            StringEncoding::Utf8,
        );
    }

    let class_vec = StringVector {
        data: vec![s4.class_name.clone()],
        encodings: vec![s4.class_encoding],
        missing: vec![false],
        attributes: class_attrs,
    };
    write_obj(&RObject::StringVector(class_vec), writer, shared)?;

    // 写入其他槽位（作为属性）
    for i in 0..s4.attributes.names.len() {
        write_pairlist_header(true, writer)?;
        let enc = s4
            .attributes
            .encodings
            .get(i)
            .copied()
            .unwrap_or(StringEncoding::Utf8);
        shared.write_symbol(&s4.attributes.names[i], enc, writer)?;
        write_obj(&s4.attributes.values[i], writer, shared)?;
    }

    // 终止
    write_header(SEXPType::NilValue, writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::environment::Environment;
    use crate::rds::external_pointer::ExternalPointer;
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
    fn test_simple_s4() {
        let s4 = S4Object {
            class_name: "MyClass".into(),
            class_encoding: StringEncoding::Utf8,
            package_name: "mypackage".into(),
            package_encoding: StringEncoding::Utf8,
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_s4(&s4, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        assert!(!data.is_empty());
        assert_eq!(data[data.len() - 1], SEXPType::NilValue as u8);
    }

    #[test]
    fn test_s4_with_slots() {
        let mut attrs = Attributes::new();
        attrs.add("slot1".into(), RObject::Null, StringEncoding::Utf8);
        let s4 = S4Object {
            class_name: "Foo".into(),
            class_encoding: StringEncoding::Utf8,
            package_name: "bar".into(),
            package_encoding: StringEncoding::Utf8,
            attributes: attrs,
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_s4(&s4, &mut buf, &mut sh, &mut nw).unwrap();
        // "class" + "slot1" = 2 symbols, reference_count 从 1 开始
        assert_eq!(sh.reference_count(), 3);
    }
}
