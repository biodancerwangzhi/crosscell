//! 符号解析
//!
//! 解析 R 符号（SYMSXP）。
//! 对应 rds2cpp 的符号解析逻辑。

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::symbol::Symbol;
use crate::rds::r_object::SymbolIndex;
use super::shared_info::SharedParseInfo;
use super::utils::{read_header, get_sexp_type, quick_extract};
use super::header::{encoding_flags, get_reference_index};

/// 解析符号体
///
/// # Requirements
/// - 15.1: 解析单字符串作为符号名
/// - 15.2: 符号名为缺失值时返回错误
/// - 15.3: 在共享信息中注册符号并返回索引
pub fn parse_symbol_body<R: Read>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
) -> Result<SymbolIndex> {
    let index = shared.request_symbol();
    
    // 读取头部，检查是 CHAR 还是 REF
    let header = read_header(reader)?;
    let sexp_type = get_sexp_type(&header);
    
    let (name, encoding, missing) = if sexp_type == SEXPType::Ref as u8 {
        // 引用类型 - 从共享信息中获取字符串或符号
        let _ref_index = get_reference_index(&header);
        shared.get_string_or_symbol(&header)?
    } else if sexp_type == SEXPType::Char as u8 {
        // 标准 CHARSXP
        let gp = header[1];
        let encoding = if (gp & encoding_flags::UTF8) != 0 { 
            crate::rds::string_encoding::StringEncoding::Utf8 
        } else if (gp & encoding_flags::LATIN1) != 0 { 
            crate::rds::string_encoding::StringEncoding::Latin1 
        } else if (gp & encoding_flags::ASCII) != 0 { 
            crate::rds::string_encoding::StringEncoding::Ascii 
        } else { 
            crate::rds::string_encoding::StringEncoding::Utf8 
        };
        
        // 读取长度
        let mut len_buf = [0u8; 4];
        reader.read_exact(&mut len_buf)?;
        let length_i32 = i32::from_be_bytes(len_buf);
        
        if length_i32 == -1 {
            // NA 字符串
            (String::new(), crate::rds::string_encoding::StringEncoding::None, true)
        } else {
            let str_len = length_i32 as usize;
            let bytes = quick_extract(reader, str_len)?;
            let value = String::from_utf8_lossy(&bytes).into_owned();
            
            // 存储到共享信息中以便后续引用
            let str_index = shared.request_string();
            shared.store_string(str_index, value.clone(), encoding, false);
            
            (value, encoding, false)
        }
    } else {
        return Err(RdsError::ParseError {
            context: "symbol".to_string(),
            message: format!("Expected CHAR (9) or Ref (255), got {}", sexp_type),
        });
    };
    
    if missing {
        return Err(RdsError::ParseError {
            context: "symbol".to_string(),
            message: "Symbol name cannot be NA".to_string(),
        });
    }
    
    let symbol = Symbol::new(name, encoding);
    shared.update_symbol(index, symbol);
    
    Ok(SymbolIndex { index })
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::string_encoding::StringEncoding;
    use crate::rds::parse::header::encoding_flags;

    fn make_char_header(encoding: StringEncoding) -> [u8; 4] {
        let gp = match encoding {
            StringEncoding::Utf8 => encoding_flags::UTF8,
            StringEncoding::Latin1 => encoding_flags::LATIN1,
            StringEncoding::Ascii => encoding_flags::ASCII,
            StringEncoding::None => 0,
        };
        [0x00, gp, 0x00, SEXPType::Char as u8]
    }

    fn make_length(len: i32) -> [u8; 4] { len.to_be_bytes() }

    #[test]
    fn test_parse_symbol_simple() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(b"x");
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let sym_idx = parse_symbol_body(&mut cursor, &mut shared).unwrap();
        assert_eq!(sym_idx.index, 0);
        
        let symbol = shared.get_symbol(0).unwrap();
        assert_eq!(symbol.name, "x");
        assert_eq!(symbol.encoding, StringEncoding::Utf8);
    }

    #[test]
    fn test_parse_symbol_longer_name() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(11));
        data.extend_from_slice(b"my_variable");
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let sym_idx = parse_symbol_body(&mut cursor, &mut shared).unwrap();
        let symbol = shared.get_symbol(sym_idx.index).unwrap();
        assert_eq!(symbol.name, "my_variable");
    }

    #[test]
    fn test_parse_symbol_na_error() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_char_header(StringEncoding::None));
        data.extend_from_slice(&make_length(-1));
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let result = parse_symbol_body(&mut cursor, &mut shared);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_multiple_symbols() {
        let mut shared = SharedParseInfo::new();
        
        // 第一个符号
        let mut data1 = Vec::new();
        data1.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data1.extend_from_slice(&make_length(1));
        data1.extend_from_slice(b"a");
        let mut cursor1 = Cursor::new(data1);
        let idx1 = parse_symbol_body(&mut cursor1, &mut shared).unwrap();
        
        // 第二个符号
        let mut data2 = Vec::new();
        data2.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data2.extend_from_slice(&make_length(1));
        data2.extend_from_slice(b"b");
        let mut cursor2 = Cursor::new(data2);
        let idx2 = parse_symbol_body(&mut cursor2, &mut shared).unwrap();
        
        assert_eq!(idx1.index, 0);
        assert_eq!(idx2.index, 1);
        assert_eq!(shared.get_symbol(0).unwrap().name, "a");
        assert_eq!(shared.get_symbol(1).unwrap().name, "b");
    }
}
