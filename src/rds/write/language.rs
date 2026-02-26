//! 语言对象写入
//!
//! 写入 R 的语言对象（函数调用表达式）。
//! 语言对象以配对列表形式写入：第一个元素是函数名符号，
//! 后续元素是参数。

use super::attributes::write_attributes;
use super::shared_info::SharedWriteInfo;
use super::utils::{write_header, write_pairlist_header};
use crate::rds::error::Result;
use crate::rds::r_object::{LanguageObject, RObject};
use crate::rds::sexp_type::SEXPType;
use std::io::Write;

/// 写入语言对象体
///
/// 语言对象格式：
/// 1. LANG 头部（带属性标记，如果有）
/// 2. 属性（如果有）
/// 3. 函数名符号
/// 4. 参数配对列表
/// 5. NILVALUE 终止
pub fn write_language<W, F>(
    lang: &LanguageObject,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
    write_obj: &mut F,
) -> Result<()>
where
    W: Write,
    F: FnMut(&RObject, &mut W, &mut SharedWriteInfo) -> Result<()>,
{
    // 写入属性（如果有）
    if !lang.attributes.is_empty() {
        write_attributes(&lang.attributes, writer, shared, write_obj)?;
    }

    // 写入函数名符号
    shared.write_symbol(&lang.function_name, lang.function_encoding, writer)?;

    // 写入参数
    for i in 0..lang.argument_values.len() {
        let has_name = lang.argument_has_name.get(i).copied().unwrap_or(false);

        // 写入 LIST 头部
        write_pairlist_header(has_name, writer)?;

        // 写入参数名（如果有）
        if has_name {
            let name = &lang.argument_names[i];
            let enc = lang.argument_encodings[i];
            shared.write_symbol(name, enc, writer)?;
        }

        // 写入参数值
        write_obj(&lang.argument_values[i], writer, shared)?;
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
    use crate::rds::string_encoding::StringEncoding;
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
    fn test_no_args() {
        let lang = LanguageObject {
            function_name: "print".into(),
            function_encoding: StringEncoding::Utf8,
            argument_values: vec![],
            argument_names: vec![],
            argument_has_name: vec![],
            argument_encodings: vec![],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_language(&lang, &mut buf, &mut sh, &mut nw).unwrap();
        let data = buf.into_inner();
        assert!(!data.is_empty());
        // Last 4 bytes = NilValue terminator
        assert_eq!(data[data.len() - 1], SEXPType::NilValue as u8);
    }

    #[test]
    fn test_with_named_args() {
        let lang = LanguageObject {
            function_name: "seq".into(),
            function_encoding: StringEncoding::Utf8,
            argument_values: vec![RObject::Null, RObject::Null],
            argument_names: vec!["from".into(), "to".into()],
            argument_has_name: vec![true, true],
            argument_encodings: vec![StringEncoding::Utf8, StringEncoding::Utf8],
            attributes: Default::default(),
        };
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut sh = mk(&s, &e, &p);
        let mut buf = Cursor::new(Vec::new());
        write_language(&lang, &mut buf, &mut sh, &mut nw).unwrap();
        // "seq", "from", "to" = 3 symbols, reference_count 从 1 开始
        assert_eq!(sh.reference_count(), 4);
    }
}
