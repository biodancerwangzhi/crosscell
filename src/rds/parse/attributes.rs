//! 属性解析
use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use crate::rds::r_object::{Attributes, RObject};
use super::utils::read_header;
use super::header::{ParsedHeader, has_tag};
use super::shared_info::SharedParseInfo;
use super::symbol::parse_symbol_body;

fn parse_header_from_bytes(h: &[u8; 4]) -> Result<ParsedHeader> {
    ParsedHeader::from_raw(*h).ok_or_else(|| RdsError::UnsupportedType(h[3]))
}

/// 解析属性标签（SYM 或 REF 类型）
/// 
/// 与 rds2cpp 保持一致：标签是符号类型，不是字符串类型。
fn parse_attribute_tag<R: Read>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
) -> Result<(String, StringEncoding, bool)> {
    let header = read_header(reader)?;
    let sexp_type = SEXPType::from_u8(header[3]);
    
    match sexp_type {
        Some(SEXPType::Sym) => {
            // 解析符号体
            let sym_idx = parse_symbol_body(reader, shared)?;
            let symbol = shared.get_symbol(sym_idx.index)
                .ok_or_else(|| RdsError::ParseError {
                    context: "attribute tag".into(),
                    message: format!("Symbol index {} not found", sym_idx.index),
                })?;
            Ok((symbol.name.clone(), symbol.encoding, false))
        }
        Some(SEXPType::Ref) => {
            // 从引用获取符号索引
            let sym_index = shared.get_symbol_index(&header)?;
            let symbol = shared.get_symbol(sym_index)
                .ok_or_else(|| RdsError::ParseError {
                    context: "attribute tag".into(),
                    message: format!("Symbol index {} not found", sym_index),
                })?;
            Ok((symbol.name.clone(), symbol.encoding, false))
        }
        _ => {
            Err(RdsError::ParseError {
                context: "attribute tag".into(),
                message: format!("Expected SYM or REF for attribute tag, got {:?}", sexp_type),
            })
        }
    }
}

pub fn parse_attributes<R: Read, F>(r: &mut R, s: &mut SharedParseInfo, f: &mut F) -> Result<Attributes>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    let h = read_header(r)?;
    let p = parse_header_from_bytes(&h)?;
    if p.sexp_type == SEXPType::NilValue || p.sexp_type == SEXPType::Nil {
        // NULL - 无属性
        return Ok(Attributes::new());
    }
    if p.sexp_type == SEXPType::Char {
        // CHAR 类型 - 解析并忽略
        // 这种情况可能发生在某些特殊的 Seurat 对象中
        use super::single_string::parse_single_string_with_header;
        let _ = parse_single_string_with_header(r, &p)?;
        return Ok(Attributes::new());
    }
    if p.sexp_type != SEXPType::List {
        return Err(RdsError::ParseError { context: "attr".into(), message: format!("Expected LIST, got {:?}", p.sexp_type) });
    }
    parse_attributes_body(r, &p, s, f)
}

pub fn parse_attributes_body<R: Read, F>(r: &mut R, h: &ParsedHeader, s: &mut SharedParseInfo, f: &mut F) -> Result<Attributes>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    let mut attrs = Attributes::new();
    let mut cur = *h;
    loop {
        if !has_tag(&cur.raw) { return Err(RdsError::ParseError { context: "attr".into(), message: "No tag".into() }); }
        let (tag_name, tag_encoding, missing) = parse_attribute_tag(r, s)?;
        if missing { return Err(RdsError::ParseError { context: "attr".into(), message: "NA name".into() }); }
        let val = f(r, s)?;
        attrs.add(tag_name, val, tag_encoding);
        let nb = read_header(r)?;
        let np = parse_header_from_bytes(&nb)?;
        if np.sexp_type == SEXPType::NilValue { 
            break; 
        }
        if np.sexp_type != SEXPType::List { return Err(RdsError::ParseError { context: "attr".into(), message: "Bad type".into() }); }
        cur = np;
    }
    Ok(attrs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::parse::header::flags;

    #[test]
    fn test_no_tag_error() {
        let d = vec![0u8; 10];
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |_: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| Ok(RObject::Null);
        let hdr = ParsedHeader::from_raw([0, 0, 0, 2]).unwrap();
        assert!(parse_attributes_body(&mut c, &hdr, &mut s, &mut pf).is_err());
    }
}
