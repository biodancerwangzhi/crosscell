//! 原子向量解析
//!
//! 解析所有原子向量类型（整数、逻辑、浮点、原始、复数、字符串）。
//! 对应 rds2cpp 的原子向量解析逻辑。
//!
//! ## RDS 原子向量格式
//!
//! 原子向量在 RDS 中的格式：
//! 1. 4 字节头部（包含类型和属性标记）
//! 2. 4 字节长度（或 8 字节大长度）
//! 3. 数据内容（根据类型不同）
//!
//! 数据使用大端序（XDR 格式）存储。

use std::io::Read;
use num_complex::Complex64;
use crate::rds::error::Result;
use crate::rds::r_object::{
    IntegerVector, LogicalVector, DoubleVector, RawVector, ComplexVector, StringVector,
};
use super::utils::{read_length, quick_extract, read_i32_be, read_f64_be, read_header, get_sexp_type};
use super::single_string::parse_single_string;

/// 解析整数向量体
///
/// 从 reader 读取整数向量的数据部分（不包括头部）。
/// 整数使用 4 字节大端序格式存储。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的整数向量
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn parse_integer_body<R: Read>(reader: &mut R) -> Result<IntegerVector> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    
    // 读取所有字节
    let bytes = quick_extract(reader, length * 4)?;
    
    // 转换为 i32（大端序）
    for chunk in bytes.chunks_exact(4) {
        let arr: [u8; 4] = chunk.try_into().unwrap();
        data.push(read_i32_be(&arr));
    }
    
    Ok(IntegerVector {
        data,
        attributes: Default::default(),
    })
}


/// 解析逻辑向量体
///
/// 从 reader 读取逻辑向量的数据部分。
/// R 的逻辑值使用 i32 表示：0 = FALSE, 1 = TRUE, NA_INTEGER = NA
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的逻辑向量
pub fn parse_logical_body<R: Read>(reader: &mut R) -> Result<LogicalVector> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    
    let bytes = quick_extract(reader, length * 4)?;
    
    for chunk in bytes.chunks_exact(4) {
        let arr: [u8; 4] = chunk.try_into().unwrap();
        data.push(read_i32_be(&arr));
    }
    
    Ok(LogicalVector {
        data,
        attributes: Default::default(),
    })
}

/// 解析浮点向量体
///
/// 从 reader 读取浮点向量的数据部分。
/// 浮点数使用 8 字节大端序 IEEE 754 格式存储。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的浮点向量
pub fn parse_double_body<R: Read>(reader: &mut R) -> Result<DoubleVector> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    
    let bytes = quick_extract(reader, length * 8)?;
    
    for chunk in bytes.chunks_exact(8) {
        let arr: [u8; 8] = chunk.try_into().unwrap();
        data.push(read_f64_be(&arr));
    }
    
    Ok(DoubleVector {
        data,
        attributes: Default::default(),
    })
}

/// 解析原始向量体
///
/// 从 reader 读取原始向量的数据部分。
/// 原始向量直接存储字节，无需字节序转换。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的原始向量
pub fn parse_raw_body<R: Read>(reader: &mut R) -> Result<RawVector> {
    let length = read_length(reader)?;
    let data = quick_extract(reader, length)?;
    
    Ok(RawVector {
        data,
        attributes: Default::default(),
    })
}


/// 解析复数向量体
///
/// 从 reader 读取复数向量的数据部分。
/// 每个复数由两个 8 字节浮点数组成（实部 + 虚部）。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的复数向量
pub fn parse_complex_body<R: Read>(reader: &mut R) -> Result<ComplexVector> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    
    // 每个复数 16 字节（8 字节实部 + 8 字节虚部）
    let bytes = quick_extract(reader, length * 16)?;
    
    for chunk in bytes.chunks_exact(16) {
        let re_arr: [u8; 8] = chunk[0..8].try_into().unwrap();
        let im_arr: [u8; 8] = chunk[8..16].try_into().unwrap();
        let re = read_f64_be(&re_arr);
        let im = read_f64_be(&im_arr);
        data.push(Complex64::new(re, im));
    }
    
    Ok(ComplexVector {
        data,
        attributes: Default::default(),
    })
}

/// 解析字符串向量体
///
/// 从 reader 读取字符串向量的数据部分。
/// 字符串向量由多个单字符串（CHARSXP）组成。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 解析后的字符串向量
pub fn parse_string_body<R: Read>(reader: &mut R) -> Result<StringVector> {
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    let mut encodings = Vec::with_capacity(length);
    let mut missing = Vec::with_capacity(length);
    
    for _ in 0..length {
        let info = parse_single_string(reader)?;
        data.push(info.value);
        encodings.push(info.encoding);
        missing.push(info.missing);
    }
    
    Ok(StringVector {
        data,
        encodings,
        missing,
        attributes: Default::default(),
    })
}

/// 解析字符串向量体（支持引用）
///
/// 从 reader 读取字符串向量的数据部分。
/// 字符串向量由多个单字符串（CHARSXP）或引用（REFSXP）组成。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源

/// 辅助函数：使用已解析的头部解析对象（用于字符串向量中的回退处理）
fn parse_object_with_header_fallback<R: Read, F>(
    reader: &mut R,
    header: &crate::rds::parse::header::ParsedHeader,
    shared: &mut crate::rds::parse::shared_info::SharedParseInfo,
    parse_fn: &mut F,
) -> crate::rds::error::Result<crate::rds::r_object::RObject>
where F: FnMut(&mut R, &mut crate::rds::parse::shared_info::SharedParseInfo) -> crate::rds::error::Result<crate::rds::r_object::RObject> {
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::r_object::RObject;
    use crate::rds::parse::header::has_attributes;
    
    match header.sexp_type {
        SEXPType::Nil | SEXPType::NilValue => Ok(RObject::Null),
        SEXPType::MissingArg | SEXPType::UnboundValue => Ok(RObject::Null),
        SEXPType::GlobalEnv | SEXPType::BaseEnv | SEXPType::EmptyEnv => Ok(RObject::Null),
        SEXPType::Char => {
            use super::single_string::parse_single_string_with_header;
            let _ = parse_single_string_with_header(reader, header)?;
            Ok(RObject::Null)
        }
        SEXPType::Sym => {
            use super::symbol::parse_symbol_body;
            let _ = parse_symbol_body(reader, shared)?;
            Ok(RObject::Null)
        }
        SEXPType::Int => {
            let mut vec = parse_integer_body(reader)?;
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                vec.attributes = parse_attributes(reader, shared, parse_fn)?;
            }
            Ok(RObject::IntegerVector(vec))
        }
        SEXPType::Real => {
            let mut vec = parse_double_body(reader)?;
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                vec.attributes = parse_attributes(reader, shared, parse_fn)?;
            }
            Ok(RObject::DoubleVector(vec))
        }
        SEXPType::Lgl => {
            let mut vec = parse_logical_body(reader)?;
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                vec.attributes = parse_attributes(reader, shared, parse_fn)?;
            }
            Ok(RObject::LogicalVector(vec))
        }
        SEXPType::Str => {
            let mut vec = parse_string_body_with_ref(reader, shared)?;
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                vec.attributes = parse_attributes(reader, shared, parse_fn)?;
            }
            Ok(RObject::StringVector(vec))
        }
        SEXPType::Vec => {
            use super::list::parse_list_body;
            let mut vec = parse_list_body(reader, shared, parse_fn)?;
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                vec.attributes = parse_attributes(reader, shared, parse_fn)?;
            }
            Ok(RObject::GenericVector(vec))
        }
        SEXPType::List => {
            use super::pairlist::parse_pairlist_body;
            let _ = parse_pairlist_body(reader, header, shared, parse_fn)?;
            Ok(RObject::Null)
        }
        SEXPType::Env => {
            use super::environment::parse_new_environment_body;
            let _ = parse_new_environment_body(reader, header, shared, parse_fn)?;
            Ok(RObject::Null)
        }
        SEXPType::Clo => {
            // Closure: formals, body, environment
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                let _ = parse_attributes(reader, shared, parse_fn)?;
            }
            let _ = parse_fn(reader, shared)?;
            let _ = parse_fn(reader, shared)?;
            let _ = parse_fn(reader, shared)?;
            Ok(RObject::Null)
        }
        SEXPType::Prom => {
            // Promise: value, expression, environment
            if has_attributes(&header.raw) {
                use super::attributes::parse_attributes;
                let _ = parse_attributes(reader, shared, parse_fn)?;
            }
            let _ = parse_fn(reader, shared)?;
            let _ = parse_fn(reader, shared)?;
            let _ = parse_fn(reader, shared)?;
            Ok(RObject::Null)
        }
        SEXPType::Lang => {
            use super::language::parse_language_body;
            let _ = parse_language_body(reader, header, shared, parse_fn)?;
            Ok(RObject::Null)
        }
        SEXPType::S4 => {
            use super::s4::parse_s4_body;
            let _ = parse_s4_body(reader, header, shared, parse_fn)?;
            Ok(RObject::Null)
        }
        SEXPType::Builtin | SEXPType::Special => {
            use super::builtin::parse_builtin_body;
            let _ = parse_builtin_body(reader)?;
            Ok(RObject::Null)
        }
        SEXPType::Extptr => {
            use super::external_pointer::parse_external_pointer_body;
            let _ = parse_external_pointer_body(reader, header, shared, parse_fn)?;
            Ok(RObject::Null)
        }
        SEXPType::Ref => {
            let _ = shared.resolve_reference(&header.raw)?;
            Ok(RObject::Null)
        }
        _ => {
            // 其他类型 - 返回错误
            Err(crate::rds::error::RdsError::ParseError {
                context: "fallback".to_string(),
                message: format!("Unsupported type in fallback: {:?}", header.sexp_type),
            })
        }
    }
}
/// * `shared` - 共享解析信息，用于解析引用
///
/// # 返回
/// 解析后的字符串向量
///
/// # 简化版本
/// 解析字符串向量体（带引用支持）
///
/// 处理 CHAR、REF、SYM 等类型。
pub fn parse_string_body_with_ref<R: Read>(
    reader: &mut R,
    shared: &mut crate::rds::parse::shared_info::SharedParseInfo,
) -> Result<StringVector> {
    use crate::rds::sexp_type::SEXPType;
    use super::header::encoding_flags;
    
    let length = read_length(reader)?;
    
    let mut data = Vec::with_capacity(length);
    let mut encodings = Vec::with_capacity(length);
    let mut missing = Vec::with_capacity(length);
    
    for i in 0..length {
        // 读取头部
        let header = read_header(reader)?;
        let sexp_type = get_sexp_type(&header);
        
        if sexp_type == SEXPType::Ref as u8 {
            // 引用类型 - 从共享信息中获取字符串
            let (value, encoding, is_missing) = shared.get_string_or_symbol(&header)?;
            data.push(value);
            encodings.push(encoding);
            missing.push(is_missing);
        } else if sexp_type == SEXPType::Char as u8 {
            // 标准 CHARSXP - 请求新槽位并解析
            let str_index = shared.request_string();
            
            // 获取编码
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
                shared.store_string(str_index, String::new(), crate::rds::string_encoding::StringEncoding::None, true);
                data.push(String::new());
                encodings.push(crate::rds::string_encoding::StringEncoding::None);
                missing.push(true);
            } else {
                let str_len = length_i32 as usize;
                let bytes = quick_extract(reader, str_len)?;
                let value = String::from_utf8_lossy(&bytes).into_owned();
                
                shared.store_string(str_index, value.clone(), encoding, false);
                data.push(value);
                encodings.push(encoding);
                missing.push(false);
            }
        } else if sexp_type == SEXPType::Sym as u8 {
            // 符号类型 - 解析符号并使用其名称
            use super::symbol::parse_symbol_body;
            let sym_idx = parse_symbol_body(reader, shared)?;
            let symbol = shared.get_symbol(sym_idx.index)
                .ok_or_else(|| crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: format!("Symbol index {} not found", sym_idx.index),
                })?;
            
            data.push(symbol.name.clone());
            encodings.push(symbol.encoding);
            missing.push(false);
        } else if sexp_type == SEXPType::NilValue as u8 || sexp_type == SEXPType::Nil as u8 {
            // NULL 值 - 作为 NA 处理
            data.push(String::new());
            encodings.push(crate::rds::string_encoding::StringEncoding::None);
            missing.push(true);
        } else if sexp_type == SEXPType::Int as u8 {
            // 整数类型 - 这不应该出现在字符串向量中
            // 如果出现，可能是 ALTREP 或其他特殊情况
            // 读取整数向量并转换为单个字符串表示
            let int_len = read_length(reader)?;
            let mut int_values = Vec::with_capacity(int_len);
            for _ in 0..int_len {
                let mut val_buf = [0u8; 4];
                reader.read_exact(&mut val_buf)?;
                let val = i32::from_be_bytes(val_buf);
                int_values.push(val);
            }
            // 将整数向量转换为单个字符串
            if int_values.len() == 1 {
                let val = int_values[0];
                if val == i32::MIN {
                    data.push(String::new());
                    encodings.push(crate::rds::string_encoding::StringEncoding::None);
                    missing.push(true);
                } else {
                    data.push(val.to_string());
                    encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                    missing.push(false);
                }
            } else {
                // 多个整数 - 转换为字符串表示
                data.push(format!("<int[{}]>", int_values.len()));
                encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                missing.push(false);
            }
        } else if sexp_type == SEXPType::Str as u8 {
            // 嵌套的字符串向量 - 递归解析并取第一个元素
            // 这种情况可能出现在某些特殊的 R 对象中
            let nested_vec = parse_string_body_with_ref(reader, shared)?;
            if nested_vec.data.is_empty() {
                data.push(String::new());
                encodings.push(crate::rds::string_encoding::StringEncoding::None);
                missing.push(true);
            } else {
                // 取第一个元素
                data.push(nested_vec.data[0].clone());
                encodings.push(nested_vec.encodings[0]);
                missing.push(nested_vec.missing[0]);
            }
        } else if sexp_type == SEXPType::GlobalEnv as u8 {
            // 全局环境
            data.push("<GlobalEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::BaseEnv as u8 {
            // 基础环境
            data.push("<BaseEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::EmptyEnv as u8 {
            // 空环境
            data.push("<EmptyEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else {
            // 其他类型 - 使用通用解析器处理以保持流位置正确
            use super::header::ParsedHeader;
            
            if let Some(parsed_header) = ParsedHeader::from_raw(header) {
                // 创建递归解析闭包
                let mut parse_fn = |r: &mut R, s: &mut crate::rds::parse::shared_info::SharedParseInfo| -> crate::rds::error::Result<crate::rds::r_object::RObject> {
                    super::object::parse_object(r, s)
                };
                
                // 使用通用解析器处理
                match parse_object_with_header_fallback(reader, &parsed_header, shared, &mut parse_fn) {
                    Ok(_) => {
                        data.push(format!("<type:{}>", sexp_type));
                        encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                        missing.push(false);
                    }
                    Err(e) => {
                        return Err(crate::rds::error::RdsError::ParseError {
                            context: "string_body".to_string(),
                            message: format!(
                                "Failed to parse unexpected type {} at element {}/{}: {}",
                                sexp_type, i, length, e
                            ),
                        });
                    }
                }
            } else {
                return Err(crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: format!(
                        "Expected CHAR/REF/SYM at string element {}/{}, got {} (0x{:02x}), header=[{:02x},{:02x},{:02x},{:02x}]",
                        i, length, sexp_type, sexp_type,
                        header[0], header[1], header[2], header[3]
                    ),
                });
            }
        }
    }
    
    Ok(StringVector {
        data,
        encodings,
        missing,
        attributes: Default::default(),
    })
}

// 保留旧的复杂版本作为备用
#[allow(dead_code)]
fn parse_string_body_with_ref_complex<R: Read>(
    reader: &mut R,
    shared: &mut crate::rds::parse::shared_info::SharedParseInfo,
) -> Result<StringVector> {
    use crate::rds::sexp_type::SEXPType;
    use super::header::encoding_flags;
    
    let length = read_length(reader)?;
    let mut data = Vec::with_capacity(length);
    let mut encodings = Vec::with_capacity(length);
    let mut missing = Vec::with_capacity(length);
    
    for _ in 0..length {
        // 读取头部
        let header = read_header(reader)?;
        let sexp_type = get_sexp_type(&header);
        
        if sexp_type == SEXPType::Ref as u8 {
            // 引用类型 - 从共享信息中获取字符串
            let (value, encoding, is_missing) = shared.get_string_or_symbol(&header)?;
            data.push(value);
            encodings.push(encoding);
            missing.push(is_missing);
        } else if sexp_type == SEXPType::Char as u8 {
            // 标准 CHARSXP - 请求新槽位并解析
            let str_index = shared.request_string();
            
            // 获取编码
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
                shared.store_string(str_index, String::new(), crate::rds::string_encoding::StringEncoding::None, true);
                data.push(String::new());
                encodings.push(crate::rds::string_encoding::StringEncoding::None);
                missing.push(true);
            } else {
                let str_len = length_i32 as usize;
                let bytes = quick_extract(reader, str_len)?;
                let value = String::from_utf8_lossy(&bytes).into_owned();
                
                shared.store_string(str_index, value.clone(), encoding, false);
                data.push(value);
                encodings.push(encoding);
                missing.push(false);
            }
        } else if sexp_type == SEXPType::Sym as u8 {
            // 符号类型 - 解析符号并使用其名称
            // 使用 parse_symbol_body 来正确处理符号名（可能是 CHAR 或 Ref）
            use super::symbol::parse_symbol_body;
            let sym_idx = parse_symbol_body(reader, shared)?;
            let symbol = shared.get_symbol(sym_idx.index)
                .ok_or_else(|| crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: format!("Symbol index {} not found", sym_idx.index),
                })?;
            
            data.push(symbol.name.clone());
            encodings.push(symbol.encoding);
            missing.push(false);
        } else if sexp_type == SEXPType::NilValue as u8 || sexp_type == SEXPType::Nil as u8 {
            // NULL 值 - 作为 NA 处理
            data.push(String::new());
            encodings.push(crate::rds::string_encoding::StringEncoding::None);
            missing.push(true);
        } else if sexp_type == SEXPType::Int as u8 {
            // 整数类型 - 可能是 ALTREP 紧凑字符串表示
            // 读取整数向量并转换为字符串
            let int_len = read_length(reader)?;
            if int_len == 1 {
                let mut val_buf = [0u8; 4];
                reader.read_exact(&mut val_buf)?;
                let val = i32::from_be_bytes(val_buf);
                if val == i32::MIN {
                    // NA 整数
                    data.push(String::new());
                    encodings.push(crate::rds::string_encoding::StringEncoding::None);
                    missing.push(true);
                } else {
                    data.push(val.to_string());
                    encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                    missing.push(false);
                }
            } else {
                // 多个整数 - 跳过并返回占位符
                for _ in 0..int_len {
                    let mut val_buf = [0u8; 4];
                    reader.read_exact(&mut val_buf)?;
                }
                data.push(format!("<int[{}]>", int_len));
                encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                missing.push(false);
            }
        } else if sexp_type == SEXPType::Env as u8 {
            // 环境类型 - 在字符串向量中遇到环境（可能是 Seurat 对象中的特殊情况）
            // 需要完整解析环境以保持流位置正确
            // 注意：头部已经被读取，所以我们需要使用 parse_new_environment_body
            use super::environment::parse_new_environment_body;
            use super::header::ParsedHeader;
            
            let parsed_header = ParsedHeader::from_raw(header)
                .ok_or_else(|| crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: "Failed to parse Env header".to_string(),
                })?;
            
            // 创建递归解析闭包
            let mut parse_fn = |r: &mut R, s: &mut crate::rds::parse::shared_info::SharedParseInfo| -> crate::rds::error::Result<crate::rds::r_object::RObject> {
                super::object::parse_object(r, s)
            };
            
            let env_idx = parse_new_environment_body(reader, &parsed_header, shared, &mut parse_fn)?;
            data.push(format!("<env:{}>", env_idx.index));
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::Prom as u8 {
            // Promise 类型 - 解析并跳过
            // Promise 结构: value, expression, environment
            // 注意：头部已经被读取
            use super::header::{ParsedHeader, has_attributes};
            
            let parsed_header = ParsedHeader::from_raw(header)
                .ok_or_else(|| crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: "Failed to parse Promise header".to_string(),
                })?;
            
            // 创建递归解析闭包
            let mut parse_fn = |r: &mut R, s: &mut crate::rds::parse::shared_info::SharedParseInfo| -> crate::rds::error::Result<crate::rds::r_object::RObject> {
                super::object::parse_object(r, s)
            };
            
            // Promise 有属性时先解析属性
            if has_attributes(&parsed_header.raw) {
                use super::attributes::parse_attributes;
                let _ = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            // 解析 value
            let _ = parse_fn(reader, shared)?;
            // 解析 expression
            let _ = parse_fn(reader, shared)?;
            // 解析 environment
            let _ = parse_fn(reader, shared)?;
            
            data.push("<promise>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::Clo as u8 {
            // Closure 类型 - 解析并跳过
            // Closure 结构: formals, body, environment
            use super::header::{ParsedHeader, has_attributes};
            
            let parsed_header = ParsedHeader::from_raw(header)
                .ok_or_else(|| crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: "Failed to parse Closure header".to_string(),
                })?;
            
            // 创建递归解析闭包
            let mut parse_fn = |r: &mut R, s: &mut crate::rds::parse::shared_info::SharedParseInfo| -> crate::rds::error::Result<crate::rds::r_object::RObject> {
                super::object::parse_object(r, s)
            };
            
            // Closure 有属性时先解析属性
            if has_attributes(&parsed_header.raw) {
                use super::attributes::parse_attributes;
                let _ = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            // 解析 formals
            let _ = parse_fn(reader, shared)?;
            // 解析 body
            let _ = parse_fn(reader, shared)?;
            // 解析 environment
            let _ = parse_fn(reader, shared)?;
            
            data.push("<closure>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::GlobalEnv as u8 {
            // 全局环境
            data.push("<GlobalEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::BaseEnv as u8 {
            // 基础环境
            data.push("<BaseEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else if sexp_type == SEXPType::EmptyEnv as u8 {
            // 空环境
            data.push("<EmptyEnv>".to_string());
            encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
            missing.push(false);
        } else {
            // 未知类型 - 可能是流位置错误
            // 尝试使用 parse_object 来解析
            use super::header::ParsedHeader;
            
            if let Some(parsed_header) = ParsedHeader::from_raw(header) {
                // 创建递归解析闭包
                let mut parse_fn = |r: &mut R, s: &mut crate::rds::parse::shared_info::SharedParseInfo| -> crate::rds::error::Result<crate::rds::r_object::RObject> {
                    super::object::parse_object(r, s)
                };
                
                // 尝试使用 parse_object_with_header 来解析
                match parse_object_with_header_fallback(reader, &parsed_header, shared, &mut parse_fn) {
                    Ok(_) => {
                        data.push(format!("<type:{}>", sexp_type));
                        encodings.push(crate::rds::string_encoding::StringEncoding::Utf8);
                        missing.push(false);
                    }
                    Err(e) => {
                        return Err(crate::rds::error::RdsError::ParseError {
                            context: "string_body".to_string(),
                            message: format!("Failed to parse type {}: {}", sexp_type, e),
                        });
                    }
                }
            } else {
                return Err(crate::rds::error::RdsError::ParseError {
                    context: "string_body".to_string(),
                    message: format!("Expected CHAR (9), Ref (255), Sym (1), Int (13), Env (4), or Nil, got {}", sexp_type),
                });
            }
        }
    }
    
    Ok(StringVector {
        data,
        encodings,
        missing,
        attributes: Default::default(),
    })
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::string_encoding::StringEncoding;
    use crate::rds::parse::header::encoding_flags;

    /// 构造长度字节（大端序）
    fn make_length(len: i32) -> [u8; 4] {
        len.to_be_bytes()
    }

    /// 构造 i32 字节（大端序）
    fn make_i32(val: i32) -> [u8; 4] {
        val.to_be_bytes()
    }

    /// 构造 f64 字节（大端序）
    fn make_f64(val: f64) -> [u8; 8] {
        val.to_be_bytes()
    }

    /// 构造 CHAR 类型头部
    fn make_char_header(encoding: StringEncoding) -> [u8; 4] {
        let gp = match encoding {
            StringEncoding::Utf8 => encoding_flags::UTF8,
            StringEncoding::Latin1 => encoding_flags::LATIN1,
            StringEncoding::Ascii => encoding_flags::ASCII,
            StringEncoding::None => 0,
        };
        [0x00, gp, 0x00, SEXPType::Char as u8]
    }

    // ========================================
    // 整数向量测试
    // ========================================

    #[test]
    fn test_parse_integer_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_integer_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
    }

    #[test]
    fn test_parse_integer_single() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(&make_i32(42));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_integer_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![42]);
    }

    #[test]
    fn test_parse_integer_multiple() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&make_i32(1));
        data.extend_from_slice(&make_i32(2));
        data.extend_from_slice(&make_i32(3));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_integer_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![1, 2, 3]);
    }

    #[test]
    fn test_parse_integer_negative() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_i32(-100));
        data.extend_from_slice(&make_i32(i32::MIN));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_integer_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![-100, i32::MIN]);
    }

    #[test]
    fn test_parse_integer_max_values() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_i32(i32::MAX));
        data.extend_from_slice(&make_i32(i32::MIN));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_integer_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![i32::MAX, i32::MIN]);
    }


    // ========================================
    // 逻辑向量测试
    // ========================================

    #[test]
    fn test_parse_logical_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_logical_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
    }

    #[test]
    fn test_parse_logical_values() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&make_i32(0));  // FALSE
        data.extend_from_slice(&make_i32(1));  // TRUE
        data.extend_from_slice(&make_i32(i32::MIN));  // NA (NA_INTEGER)
        
        let mut cursor = Cursor::new(data);
        let vec = parse_logical_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![0, 1, i32::MIN]);
    }

    // ========================================
    // 浮点向量测试
    // ========================================

    #[test]
    fn test_parse_double_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_double_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
    }

    #[test]
    fn test_parse_double_single() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(&make_f64(3.14159));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_double_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data.len(), 1);
        assert!((vec.data[0] - 3.14159).abs() < 1e-10);
    }

    #[test]
    fn test_parse_double_multiple() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&make_f64(1.0));
        data.extend_from_slice(&make_f64(-2.5));
        data.extend_from_slice(&make_f64(0.0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_double_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![1.0, -2.5, 0.0]);
    }

    #[test]
    fn test_parse_double_special_values() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_f64(f64::INFINITY));
        data.extend_from_slice(&make_f64(f64::NEG_INFINITY));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_double_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data[0], f64::INFINITY);
        assert_eq!(vec.data[1], f64::NEG_INFINITY);
    }


    // ========================================
    // 原始向量测试
    // ========================================

    #[test]
    fn test_parse_raw_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_raw_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
    }

    #[test]
    fn test_parse_raw_bytes() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(5));
        data.extend_from_slice(&[0x01, 0x02, 0x03, 0x04, 0x05]);
        
        let mut cursor = Cursor::new(data);
        let vec = parse_raw_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![0x01, 0x02, 0x03, 0x04, 0x05]);
    }

    #[test]
    fn test_parse_raw_all_values() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&[0x00, 0xFF, 0x80]);
        
        let mut cursor = Cursor::new(data);
        let vec = parse_raw_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec![0x00, 0xFF, 0x80]);
    }

    // ========================================
    // 复数向量测试
    // ========================================

    #[test]
    fn test_parse_complex_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_complex_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
    }

    #[test]
    fn test_parse_complex_single() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(&make_f64(3.0));  // 实部
        data.extend_from_slice(&make_f64(4.0));  // 虚部
        
        let mut cursor = Cursor::new(data);
        let vec = parse_complex_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data.len(), 1);
        assert_eq!(vec.data[0].re, 3.0);
        assert_eq!(vec.data[0].im, 4.0);
    }

    #[test]
    fn test_parse_complex_multiple() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_f64(1.0));
        data.extend_from_slice(&make_f64(2.0));
        data.extend_from_slice(&make_f64(-1.0));
        data.extend_from_slice(&make_f64(-2.0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_complex_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data.len(), 2);
        assert_eq!(vec.data[0], Complex64::new(1.0, 2.0));
        assert_eq!(vec.data[1], Complex64::new(-1.0, -2.0));
    }


    // ========================================
    // 字符串向量测试
    // ========================================

    #[test]
    fn test_parse_string_empty() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let vec = parse_string_body(&mut cursor).unwrap();
        
        assert!(vec.data.is_empty());
        assert!(vec.encodings.is_empty());
        assert!(vec.missing.is_empty());
    }

    #[test]
    fn test_parse_string_single() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(1));
        // 添加一个 CHAR 对象
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(5));
        data.extend_from_slice(b"hello");
        
        let mut cursor = Cursor::new(data);
        let vec = parse_string_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec!["hello"]);
        assert_eq!(vec.encodings, vec![StringEncoding::Utf8]);
        assert_eq!(vec.missing, vec![false]);
    }

    #[test]
    fn test_parse_string_multiple() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        
        // 第一个字符串
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(b"foo");
        
        // 第二个字符串
        data.extend_from_slice(&make_char_header(StringEncoding::Ascii));
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(b"bar");
        
        let mut cursor = Cursor::new(data);
        let vec = parse_string_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec!["foo", "bar"]);
        assert_eq!(vec.encodings, vec![StringEncoding::Utf8, StringEncoding::Ascii]);
        assert_eq!(vec.missing, vec![false, false]);
    }

    #[test]
    fn test_parse_string_with_na() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        
        // 第一个字符串
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(b"a");
        
        // NA 字符串
        data.extend_from_slice(&make_char_header(StringEncoding::None));
        data.extend_from_slice(&make_length(-1));
        
        // 第三个字符串
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(b"b");
        
        let mut cursor = Cursor::new(data);
        let vec = parse_string_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec!["a", "", "b"]);
        assert_eq!(vec.missing, vec![false, true, false]);
    }

    #[test]
    fn test_parse_string_unicode() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(1));
        
        let unicode_str = "你好";
        let bytes = unicode_str.as_bytes();
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(bytes.len() as i32));
        data.extend_from_slice(bytes);
        
        let mut cursor = Cursor::new(data);
        let vec = parse_string_body(&mut cursor).unwrap();
        
        assert_eq!(vec.data, vec!["你好"]);
    }

    // ========================================
    // 错误处理测试
    // ========================================

    #[test]
    fn test_parse_integer_truncated() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&make_i32(1));
        // 缺少 2 个整数
        
        let mut cursor = Cursor::new(data);
        let result = parse_integer_body(&mut cursor);
        
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_double_truncated() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_f64(1.0));
        // 缺少 1 个浮点数
        
        let mut cursor = Cursor::new(data);
        let result = parse_double_body(&mut cursor);
        
        assert!(result.is_err());
    }
}
