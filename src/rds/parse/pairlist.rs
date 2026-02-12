//! 配对列表解析
//!
//! 参考 rds2cpp 的 parse_pairlist.hpp 实现。
//!
//! ## Pairlist 结构
//!
//! Pairlist 是 R 中的链表结构，每个节点包含：
//! - 可选的属性（第一个节点）
//! - 可选的标签（符号类型）
//! - 值（任意 R 对象）
//! - CDR（下一个节点或 NILVALUE_ 终止符）

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use crate::rds::r_object::{PairList, RObject};
use super::utils::read_header;
use super::header::{ParsedHeader, has_tag, has_attributes};
use super::shared_info::SharedParseInfo;
use super::symbol::parse_symbol_body;

fn parse_header_from_bytes(h: &[u8; 4]) -> Result<ParsedHeader> {
    ParsedHeader::from_raw(*h).ok_or_else(|| RdsError::UnsupportedType(h[3]))
}

/// 解析配对列表标签（SYM 或 REF 类型）
/// 
/// 与 rds2cpp 保持一致：标签是符号类型，不是字符串类型。
fn parse_pairlist_tag<R: Read>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
) -> Result<(String, StringEncoding)> {
    let header = read_header(reader)?;
    let sexp_type = SEXPType::from_u8(header[3]);
    
    match sexp_type {
        Some(SEXPType::Sym) => {
            // 解析符号体
            let sym_idx = parse_symbol_body(reader, shared)?;
            let symbol = shared.get_symbol(sym_idx.index)
                .ok_or_else(|| RdsError::ParseError {
                    context: "pairlist tag".into(),
                    message: format!("Symbol index {} not found", sym_idx.index),
                })?;
            Ok((symbol.name.clone(), symbol.encoding))
        }
        Some(SEXPType::Ref) => {
            // 从引用获取符号索引
            let sym_index = shared.get_symbol_index(&header)?;
            let symbol = shared.get_symbol(sym_index)
                .ok_or_else(|| RdsError::ParseError {
                    context: "pairlist tag".into(),
                    message: format!("Symbol index {} not found", sym_index),
                })?;
            Ok((symbol.name.clone(), symbol.encoding))
        }
        _ => {
            Err(RdsError::ParseError {
                context: "pairlist tag".into(),
                message: format!("Expected SYM or REF for pairlist tag, got {:?}", sexp_type),
            })
        }
    }
}

/// 递归解析 pairlist 节点
/// 
/// 参考 R 的 ReadItem_Iterative 实现。
/// R 的 pairlist CDR 可以是任何类型，不仅仅是 LIST 或 NILVALUE。
fn recursive_parse<R: Read, F>(
    reader: &mut R,
    output: &mut PairList,
    header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
    is_first: bool,
) -> Result<()>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 处理第一个节点的属性
    if is_first && has_attributes(&header.raw) {
        use super::attributes::parse_attributes;
        output.attributes = parse_attributes(reader, shared, parse_fn)?;
    }
    
    // 解析标签（如果有）
    let has_tag_flag = has_tag(&header.raw);
    output.has_tag.push(has_tag_flag);
    
    if has_tag_flag {
        let (tag_name, tag_encoding) = parse_pairlist_tag(reader, shared)?;
        output.tag_names.push(tag_name);
        output.tag_encodings.push(tag_encoding);
    } else {
        output.tag_names.push(String::new());
        output.tag_encodings.push(StringEncoding::default());
    }
    
    // 解析值 (CAR)
    let value = parse_fn(reader, shared)?;
    output.data.push(value);
    
    // 读取下一个节点的头部 (CDR)
    // R 的 pairlist CDR 可以是任何类型，不仅仅是 LIST 或 NILVALUE
    let next_header = read_header(reader)?;
    let next_type = next_header[3];
    
    // 检查是否是终止符 (NILVALUE_ = 254)
    if next_type == SEXPType::NilValue as u8 {
        return Ok(());
    }
    
    // 检查是否是 NULL (NILSXP = 0)
    if next_type == SEXPType::Nil as u8 {
        return Ok(());
    }
    
    // 检查是否是下一个 pairlist 节点 (LIST = 2)
    if next_type == SEXPType::List as u8 {
        let next_parsed = parse_header_from_bytes(&next_header)?;
        return recursive_parse(reader, output, &next_parsed, shared, parse_fn, false);
    }
    
    // CDR 是其他类型 - 解析它但不添加到 data
    // 这种情况在 R 中是合法的，虽然不常见
    // CDR 不是 pairlist 的一部分，它是"剩余"部分
    let next_parsed = parse_header_from_bytes(&next_header)?;
    use super::object::parse_object_with_header;
    let _cdr_value = parse_object_with_header(reader, &next_parsed, shared)?;
    // 不添加到 output.data - CDR 不是 pairlist 元素
    
    Ok(())
}

/// 解析配对列表体
/// 
/// 参考 rds2cpp 的 parse_pairlist_body
pub fn parse_pairlist_body<R: Read, F>(
    reader: &mut R,
    header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<PairList>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    let mut output = PairList::default();
    recursive_parse(reader, &mut output, header, shared, parse_fn, true)?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::parse::header::flags;

    fn list_hdr(tag: bool) -> [u8; 4] {
        [0, 0, if tag { flags::HAS_TAG } else { 0 }, SEXPType::List as u8]
    }
    
    fn nil_hdr() -> [u8; 4] {
        [0, 0, 0, SEXPType::NilValue as u8]
    }

    #[test]
    fn test_single_untagged() {
        let mut d = Vec::new();
        d.extend_from_slice(&[0, 0, 0, SEXPType::NilValue as u8]); // value = NULL
        d.extend_from_slice(&nil_hdr()); // terminator
        
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |r: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| -> Result<RObject> {
            let mut h = [0u8; 4];
            r.read_exact(&mut h)?;
            Ok(RObject::Null)
        };
        let hdr = ParsedHeader::from_raw(list_hdr(false)).unwrap();
        let pl = parse_pairlist_body(&mut c, &hdr, &mut s, &mut pf).unwrap();
        assert_eq!(pl.data.len(), 1);
        assert!(!pl.has_tag[0]);
    }
}
