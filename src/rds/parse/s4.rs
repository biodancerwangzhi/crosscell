//! S4 对象解析
use super::header::ParsedHeader;
use super::shared_info::SharedParseInfo;
use super::utils::read_header;
use crate::rds::error::Result;
use crate::rds::r_object::{RObject, S4Object};
use crate::rds::sexp_type::SEXPType;
use std::io::Read;

fn parse_header_from_bytes(h: &[u8; 4]) -> Result<ParsedHeader> {
    ParsedHeader::from_raw(*h).ok_or_else(|| crate::rds::error::RdsError::UnsupportedType(h[3]))
}

/// 解析 S4 对象体
///
/// 参考 rds2cpp parse_s4.hpp
/// S4 对象的槽位存储为 pairlist
pub fn parse_s4_body<R: Read, F>(
    reader: &mut R,
    _header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<S4Object>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    let mut s4 = S4Object::default();

    // 读取槽位 pairlist 的头部
    let slot_header_bytes = read_header(reader)?;
    let slot_header = parse_header_from_bytes(&slot_header_bytes)?;

    // 检查是否为 NilValue（空 S4 对象）
    if slot_header.sexp_type == SEXPType::NilValue || slot_header.sexp_type == SEXPType::Nil {
        return Ok(s4);
    }

    // 槽位应该是 pairlist
    if slot_header.sexp_type != SEXPType::List {
        // 非 LIST 类型 - 可能是特殊情况，返回空 S4
        return Ok(s4);
    }

    // 使用 pairlist 解析器解析槽位
    use super::pairlist::parse_pairlist_body;
    let pairlist = parse_pairlist_body(reader, &slot_header, shared, parse_fn)?;

    // 处理槽位
    for i in 0..pairlist.data.len() {
        let tag_name = &pairlist.tag_names[i];
        let tag_encoding = pairlist.tag_encodings[i];
        let val = &pairlist.data[i];

        if tag_name == "class" {
            if let RObject::StringVector(sv) = val {
                if !sv.data.is_empty() {
                    s4.class_name = sv.data[0].clone();
                    s4.class_encoding = sv.encodings.get(0).copied().unwrap_or_default();
                }
                // 从 class 字符串向量的属性中提取 package 名称
                if let Some(pkg_idx) = sv.attributes.names.iter().position(|n| n == "package") {
                    if let Some(RObject::StringVector(pkg_sv)) = sv.attributes.values.get(pkg_idx) {
                        if !pkg_sv.data.is_empty() {
                            s4.package_name = pkg_sv.data[0].clone();
                            s4.package_encoding =
                                pkg_sv.encodings.get(0).copied().unwrap_or_default();
                        }
                    }
                }
            }
        } else {
            s4.attributes
                .add(tag_name.clone(), val.clone(), tag_encoding);
        }
    }

    Ok(s4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_empty_s4() {
        let d = [0u8, 0, 0, SEXPType::NilValue as u8].to_vec();
        let mut c = Cursor::new(d);
        let mut s = SharedParseInfo::new();
        let mut pf = |_: &mut Cursor<Vec<u8>>, _: &mut SharedParseInfo| Ok(RObject::Null);
        let hdr = ParsedHeader::from_raw([0, 0, 1, SEXPType::S4 as u8]).unwrap();
        let s4 = parse_s4_body(&mut c, &hdr, &mut s, &mut pf).unwrap();
        assert!(s4.class_name.is_empty());
    }
}
