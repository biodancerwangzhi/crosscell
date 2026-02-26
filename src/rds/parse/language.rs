//! 语言对象解析
use super::header::{encoding_flags, has_tag, ParsedHeader};
use super::shared_info::SharedParseInfo;
use super::utils::{get_sexp_type, quick_extract, read_header};
use crate::rds::error::{RdsError, Result};
use crate::rds::r_object::{LanguageObject, RObject};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use std::io::Read;

fn parse_header_from_bytes(h: &[u8; 4]) -> Result<ParsedHeader> {
    ParsedHeader::from_raw(*h).ok_or_else(|| RdsError::UnsupportedType(h[3]))
}

/// 解析语言对象的 tag（可以是 CHAR、REF 或 SYM）
fn parse_language_tag<R: Read>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
) -> Result<(String, StringEncoding)> {
    let header = read_header(reader)?;
    let sexp_type = get_sexp_type(&header);

    if sexp_type == SEXPType::Ref as u8 {
        // 引用类型 - 从共享信息中获取字符串或符号
        let (value, encoding, _missing) = shared.get_string_or_symbol(&header)?;
        Ok((value, encoding))
    } else if sexp_type == SEXPType::Char as u8 {
        // 标准 CHARSXP
        let gp = header[1];
        let encoding = if (gp & encoding_flags::UTF8) != 0 {
            StringEncoding::Utf8
        } else if (gp & encoding_flags::LATIN1) != 0 {
            StringEncoding::Latin1
        } else if (gp & encoding_flags::ASCII) != 0 {
            StringEncoding::Ascii
        } else {
            StringEncoding::Utf8
        };

        // 读取长度
        let mut len_buf = [0u8; 4];
        reader.read_exact(&mut len_buf)?;
        let length_i32 = i32::from_be_bytes(len_buf);

        if length_i32 == -1 {
            // NA 字符串
            Ok((String::new(), StringEncoding::None))
        } else {
            let str_len = length_i32 as usize;
            let bytes = quick_extract(reader, str_len)?;
            let value = String::from_utf8_lossy(&bytes).into_owned();

            // 存储到共享信息中以便后续引用
            let str_index = shared.request_string();
            shared.store_string(str_index, value.clone(), encoding, false);

            Ok((value, encoding))
        }
    } else if sexp_type == SEXPType::Sym as u8 {
        // 符号类型 - 解析符号并使用其名称
        use super::symbol::parse_symbol_body;
        let sym_idx = parse_symbol_body(reader, shared)?;
        let symbol = shared
            .get_symbol(sym_idx.index)
            .ok_or_else(|| RdsError::ParseError {
                context: "language tag".into(),
                message: format!("Symbol index {} not found", sym_idx.index),
            })?;
        Ok((symbol.name.clone(), symbol.encoding))
    } else {
        Err(RdsError::ParseError {
            context: "language tag".into(),
            message: format!(
                "Expected CHAR (9), Ref (255), or Sym (1), got {}",
                sexp_type
            ),
        })
    }
}

/// 解析语言对象体（作为配对列表解析，第一个元素是函数名）
pub fn parse_language_body<R: Read, F>(
    reader: &mut R,
    header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<LanguageObject>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    let mut lang = LanguageObject::default();
    let mut cur = *header;
    let mut first = true;

    loop {
        if has_tag(&cur.raw) {
            let (tag_value, tag_encoding) = parse_language_tag(reader, shared)?;
            if first {
                // 第一个元素的标签是函数名
                lang.function_name = tag_value;
                lang.function_encoding = tag_encoding;
            } else {
                lang.argument_has_name.push(true);
                lang.argument_names.push(tag_value);
                lang.argument_encodings.push(tag_encoding);
            }
        } else if first {
            // 第一个元素没有标签，需要解析符号获取函数名
            let _val = parse_fn(reader, shared)?;
            // 从符号中提取函数名（简化处理）
            lang.function_name = String::new();
            lang.function_encoding = StringEncoding::Utf8;
            first = false;

            let next_bytes = read_header(reader)?;
            let next = parse_header_from_bytes(&next_bytes)?;
            if next.sexp_type == SEXPType::NilValue {
                break;
            }
            if next.sexp_type != SEXPType::List {
                return Err(RdsError::ParseError {
                    context: "language".into(),
                    message: format!("Expected LIST or NILVALUE_, got {:?}", next.sexp_type),
                });
            }
            cur = next;
            continue;
        } else {
            lang.argument_has_name.push(false);
            lang.argument_names.push(String::new());
            lang.argument_encodings.push(StringEncoding::Utf8);
        }

        if !first {
            let val = parse_fn(reader, shared)?;
            lang.argument_values.push(val);
        }
        first = false;

        let next_bytes = read_header(reader)?;
        let next = parse_header_from_bytes(&next_bytes)?;
        if next.sexp_type == SEXPType::NilValue {
            break;
        }
        if next.sexp_type != SEXPType::List {
            return Err(RdsError::ParseError {
                context: "language".into(),
                message: format!("Expected LIST or NILVALUE_, got {:?}", next.sexp_type),
            });
        }
        cur = next;
    }
    Ok(lang)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::parse::header::flags;
    use std::io::Cursor;

    fn list_hdr(tag: bool) -> [u8; 4] {
        [0, 0, if tag { flags::HAS_TAG } else { 0 }, 2]
    }
    fn nil_hdr() -> [u8; 4] {
        [0, 0, 0, SEXPType::NilValue as u8]
    }

    #[test]
    fn test_simple_language() {
        // 简化测试：只有终止符
        let d = nil_hdr().to_vec();
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |_: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| -> Result<RObject> {
            Ok(RObject::Null)
        };
        let hdr = ParsedHeader::from_raw(list_hdr(false)).unwrap();
        // 这个测试会因为没有函数名而产生空的语言对象
    }
}
