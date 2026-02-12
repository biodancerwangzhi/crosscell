//! ALTREP 解析
//!
//! ALTREP (Alternative Representation) 是 R 3.5+ 引入的一种优化机制，
//! 允许向量使用替代的内部表示。
//!
//! 参考 rds2cpp 的 parse_altrep.hpp 实现。
//!
//! ## ALTREP 结构
//!
//! ALTREP 对象的结构：
//! 1. LIST 头部（pairlist）
//! 2. Pairlist 包含类型信息（第一个元素是符号，指示 ALTREP 类型）
//! 3. 根据类型不同，后续结构不同
//!
//! ## 支持的 ALTREP 类型
//!
//! - `compact_intseq`: 紧凑整数序列（如 1:1000）
//! - `wrap_integer`: 带属性的整数向量包装器
//! - `wrap_real`: 带属性的浮点向量包装器
//! - `deferred_string`: 延迟字符串转换

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::r_object::{RObject, IntegerVector, DoubleVector, StringVector};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use super::shared_info::SharedParseInfo;
use super::utils::read_header;
use super::header::ParsedHeader;
use super::atomic::{parse_integer_body, parse_double_body};

/// 解析 ALTREP 体
///
/// 参考 rds2cpp 的 parse_altrep_body 实现。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
/// * `shared` - 共享解析信息
/// * `parse_fn` - 递归解析函数
///
/// # 返回
/// 解析后的 RObject
pub fn parse_altrep_body<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<RObject>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 1. 读取 pairlist 头部
    let header = read_header(reader)?;
    let sexp_type = header[3];
    
    if sexp_type != SEXPType::List as u8 {
        return Err(RdsError::ParseError {
            context: "altrep".to_string(),
            message: format!("Expected ALTREP description to be a pairlist (LIST), got type {}", sexp_type),
        });
    }
    
    let parsed_header = ParsedHeader::from_raw(header)
        .ok_or_else(|| RdsError::ParseError {
            context: "altrep".to_string(),
            message: "Failed to parse pairlist header".to_string(),
        })?;
    
    // 2. 解析 pairlist 以获取类型信息
    use super::pairlist::parse_pairlist_body;
    let pairlist = parse_pairlist_body(reader, &parsed_header, shared, parse_fn)?;
    
    // 3. 从 pairlist 中提取类型符号
    if pairlist.data.is_empty() {
        return Err(RdsError::ParseError {
            context: "altrep".to_string(),
            message: "ALTREP pairlist is empty, expected type specification".to_string(),
        });
    }
    
    // 第一个元素应该是符号索引或字符串向量（引用到 CHARSXP）
    let altrep_type = match &pairlist.data[0] {
        RObject::SymbolIndex(sym_idx) => {
            let name = shared.get_symbol(sym_idx.index)
                .map(|s| s.name.clone())
                .unwrap_or_else(|| "<unknown>".to_string());
            name
        }
        RObject::StringVector(sv) => {
            // R 可能将符号序列化为 CHARSXP 引用
            // 在这种情况下，我们从字符串向量中提取类名
            if sv.data.is_empty() {
                return Err(RdsError::ParseError {
                    context: "altrep".to_string(),
                    message: "ALTREP type specification string is empty".to_string(),
                });
            }
            let name = sv.data[0].clone();
            name
        }
        _ => {
            return Err(RdsError::ParseError {
                context: "altrep".to_string(),
                message: format!("Expected type specification symbol or string in ALTREP description, got {:?}", pairlist.data[0]),
            });
        }
    };
    
    // 4. 根据类型分发到具体的解析器
    // 注意：R 的 ALTREP 序列化格式是 flags, info, state, attr
    // 但对于 wrap_* 类型，attr 是作为 state 的一部分被读取的
    // （在 parse_wrap_integer/parse_wrap_real 中读取）
    // 所以这里不需要单独读取 attr
    match altrep_type.as_str() {
        "wrap_integer" => parse_wrap_integer(reader, shared, parse_fn),
        "wrap_real" => parse_wrap_real(reader, shared, parse_fn),
        "compact_intseq" => parse_compact_intseq(reader),
        "compact_realseq" => parse_compact_realseq(reader),
        "deferred_string" => parse_deferred_string(reader, shared, parse_fn),
        _ => {
            // 未知的 ALTREP 类型 - 尝试通用解析
            // 某些 ALTREP 类型可能有不同的结构
            Err(RdsError::ParseError {
                context: "altrep".to_string(),
                message: format!("Unrecognized ALTREP type: '{}'", altrep_type),
            })
        }
    }
}

/// 解析 wrap_integer ALTREP
///
/// 结构：
/// 1. LIST 头部（pairlist 节点）
/// 2. 第一个元素：包装的整数向量（CAR）
/// 3. 第二个元素：元数据（长度为 2 的整数向量）
/// 4. 属性（pairlist 或 NULL）- 这些属性会被添加到 contents 的属性中
///
/// 注意：rds2cpp 不解析完整的 pairlist 结构，而是直接读取对象
fn parse_wrap_integer<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<RObject>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 读取 pairlist 头部
    let plist_header = read_header(reader)?;
    if plist_header[3] != SEXPType::List as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_integer".to_string(),
            message: format!("Expected pairlist in wrap_integer payload, got type {}", plist_header[3]),
        });
    }
    
    // 直接解析包装的整数向量（CAR）
    let contents = parse_fn(reader, shared)?;
    let mut int_vec = match &contents {
        RObject::IntegerVector(v) => {
            v.clone()
        }
        _ => {
            return Err(RdsError::ParseError {
                context: "wrap_integer".to_string(),
                message: format!("Expected IntegerVector in wrap_integer payload, got {:?}", contents),
            });
        }
    };
    
    // 解析元数据（长度为 2 的整数向量）
    let meta_header = read_header(reader)?;
    if meta_header[3] != SEXPType::Int as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_integer".to_string(),
            message: format!("Expected integer vector for metadata, got type {}", meta_header[3]),
        });
    }
    let metadata = parse_integer_body(reader)?;
    if metadata.data.len() != 2 {
        return Err(RdsError::ParseError {
            context: "wrap_integer".to_string(),
            message: format!("wrap_integer metadata should be length 2, got {}", metadata.data.len()),
        });
    }
    
    // 解析属性 - 这些属性会被添加到 contents 的属性中（参考 rds2cpp）
    let attr_header = read_header(reader)?;
    if attr_header[3] == SEXPType::List as u8 {
        let parsed_header = ParsedHeader::from_raw(attr_header)
            .ok_or_else(|| RdsError::ParseError {
                context: "wrap_integer".to_string(),
                message: "Failed to parse attribute header".to_string(),
            })?;
        use super::attributes::parse_attributes_body;
        let wrapper_attrs = parse_attributes_body(reader, &parsed_header, shared, parse_fn)?;
        
        // 将 wrapper 的属性添加到 contents 的属性中
        for i in 0..wrapper_attrs.names.len() {
            int_vec.attributes.names.push(wrapper_attrs.names[i].clone());
            int_vec.attributes.encodings.push(wrapper_attrs.encodings[i]);
            int_vec.attributes.values.push(wrapper_attrs.values[i].clone());
        }
    } else if attr_header[3] != SEXPType::NilValue as u8 && attr_header[3] != SEXPType::Nil as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_integer".to_string(),
            message: format!("Expected pairlist or NULL for attributes, got type {}", attr_header[3]),
        });
    }
    
    Ok(RObject::IntegerVector(int_vec))
}

/// 解析 wrap_real ALTREP
///
/// 结构与 wrap_integer 相同（参考 rds2cpp 的 parse_attribute_wrapper）：
/// 1. LIST 头部（CONS cell）
/// 2. 第一个元素：包装的浮点向量（CAR）
/// 3. 第二个元素：元数据（长度为 2 的整数向量）- 直接跟在 contents 后面
/// 4. 属性（pairlist 或 NULL）- 这些属性会被添加到 contents 的属性中
fn parse_wrap_real<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<RObject>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 读取 pairlist 头部
    let plist_header = read_header(reader)?;
    if plist_header[3] != SEXPType::List as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_real".to_string(),
            message: format!("Expected pairlist in wrap_real payload, got type {}", plist_header[3]),
        });
    }
    
    // 直接解析包装的浮点向量（CAR）
    let contents = parse_fn(reader, shared)?;
    let mut double_vec = match &contents {
        RObject::DoubleVector(v) => {
            v.clone()
        }
        _ => {
            return Err(RdsError::ParseError {
                context: "wrap_real".to_string(),
                message: format!("Expected DoubleVector in wrap_real payload, got {:?}", contents),
            });
        }
    };
    
    // 解析元数据（长度为 2 的整数向量）
    let meta_header = read_header(reader)?;
    if meta_header[3] != SEXPType::Int as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_real".to_string(),
            message: format!("Expected integer vector for metadata, got type {}", meta_header[3]),
        });
    }
    let metadata = parse_integer_body(reader)?;
    if metadata.data.len() != 2 {
        return Err(RdsError::ParseError {
            context: "wrap_real".to_string(),
            message: format!("wrap_real metadata should be length 2, got {}", metadata.data.len()),
        });
    }
    
    // 解析属性 - 这些属性会被添加到 contents 的属性中（参考 rds2cpp）
    let attr_header = read_header(reader)?;
    if attr_header[3] == SEXPType::List as u8 {
        let parsed_header = ParsedHeader::from_raw(attr_header)
            .ok_or_else(|| RdsError::ParseError {
                context: "wrap_real".to_string(),
                message: "Failed to parse attribute header".to_string(),
            })?;
        use super::attributes::parse_attributes_body;
        let wrapper_attrs = parse_attributes_body(reader, &parsed_header, shared, parse_fn)?;
        
        // 将 wrapper 的属性添加到 contents 的属性中
        for i in 0..wrapper_attrs.names.len() {
            double_vec.attributes.names.push(wrapper_attrs.names[i].clone());
            double_vec.attributes.encodings.push(wrapper_attrs.encodings[i]);
            double_vec.attributes.values.push(wrapper_attrs.values[i].clone());
        }
    } else if attr_header[3] != SEXPType::NilValue as u8 && attr_header[3] != SEXPType::Nil as u8 {
        return Err(RdsError::ParseError {
            context: "wrap_real".to_string(),
            message: format!("Expected pairlist or NULL for attributes, got type {}", attr_header[3]),
        });
    }
    
    Ok(RObject::DoubleVector(double_vec))
}

/// 解析 compact_intseq ALTREP
///
/// 紧凑整数序列，如 1:1000
///
/// 结构：
/// 1. REAL 头部
/// 2. 长度为 3 的浮点向量：[length, start, step]
/// 3. 终止符（254）
fn parse_compact_intseq<R: Read>(reader: &mut R) -> Result<RObject> {
    // 读取 REAL 头部
    let header = read_header(reader)?;
    if header[3] != SEXPType::Real as u8 {
        return Err(RdsError::ParseError {
            context: "compact_intseq".to_string(),
            message: format!("Expected REAL for sequence info, got type {}", header[3]),
        });
    }
    
    // 解析序列信息
    let info = parse_double_body(reader)?;
    if info.data.len() != 3 {
        return Err(RdsError::ParseError {
            context: "compact_intseq".to_string(),
            message: format!("Expected 3 elements for sequence info, got {}", info.data.len()),
        });
    }
    
    let len = info.data[0] as usize;
    let start = info.data[1];
    let step = info.data[2];
    
    // 生成整数序列
    let mut data = Vec::with_capacity(len);
    let mut current = start;
    for _ in 0..len {
        data.push(current as i32);
        current += step;
    }
    
    // 读取终止符
    let terminator = read_header(reader)?;
    if terminator[3] != 254 {
        return Err(RdsError::ParseError {
            context: "compact_intseq".to_string(),
            message: format!("Expected terminator (254), got {}", terminator[3]),
        });
    }
    
    Ok(RObject::IntegerVector(IntegerVector {
        data,
        attributes: Default::default(),
    }))
}

/// 解析 compact_realseq ALTREP
///
/// 紧凑浮点序列
fn parse_compact_realseq<R: Read>(reader: &mut R) -> Result<RObject> {
    // 读取 REAL 头部
    let header = read_header(reader)?;
    if header[3] != SEXPType::Real as u8 {
        return Err(RdsError::ParseError {
            context: "compact_realseq".to_string(),
            message: format!("Expected REAL for sequence info, got type {}", header[3]),
        });
    }
    
    // 解析序列信息
    let info = parse_double_body(reader)?;
    if info.data.len() != 3 {
        return Err(RdsError::ParseError {
            context: "compact_realseq".to_string(),
            message: format!("Expected 3 elements for sequence info, got {}", info.data.len()),
        });
    }
    
    let len = info.data[0] as usize;
    let start = info.data[1];
    let step = info.data[2];
    
    // 生成浮点序列
    let mut data = Vec::with_capacity(len);
    let mut current = start;
    for _ in 0..len {
        data.push(current);
        current += step;
    }
    
    // 读取终止符
    let terminator = read_header(reader)?;
    if terminator[3] != 254 {
        return Err(RdsError::ParseError {
            context: "compact_realseq".to_string(),
            message: format!("Expected terminator (254), got {}", terminator[3]),
        });
    }
    
    Ok(RObject::DoubleVector(DoubleVector {
        data,
        attributes: Default::default(),
    }))
}

/// 解析 deferred_string ALTREP
///
/// 延迟字符串转换，将整数或浮点向量延迟转换为字符串
///
/// 结构：
/// 1. LIST 头部
/// 2. 第一个元素：要转换的向量（INT 或 REAL）
/// 3. 第二个元素：元数据（长度为 1 的整数向量）
/// 4. 终止符（NILVALUE）
fn parse_deferred_string<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<RObject>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 读取 pairlist 头部
    let plist_header = read_header(reader)?;
    if plist_header[3] != SEXPType::List as u8 {
        return Err(RdsError::ParseError {
            context: "deferred_string".to_string(),
            message: format!("Expected pairlist in deferred_string payload, got type {}", plist_header[3]),
        });
    }
    
    // 解析要转换的向量
    let contents = parse_fn(reader, shared)?;
    
    let output = match contents {
        RObject::IntegerVector(iv) => {
            let n = iv.data.len();
            let mut data = Vec::with_capacity(n);
            let mut encodings = Vec::with_capacity(n);
            let mut missing = Vec::with_capacity(n);
            
            for val in iv.data {
                if val == i32::MIN {
                    // NA 值
                    data.push(String::new());
                    encodings.push(StringEncoding::None);
                    missing.push(true);
                } else {
                    data.push(val.to_string());
                    encodings.push(StringEncoding::Ascii);
                    missing.push(false);
                }
            }
            
            StringVector {
                data,
                encodings,
                missing,
                attributes: Default::default(),
            }
        }
        RObject::DoubleVector(dv) => {
            let n = dv.data.len();
            let mut data = Vec::with_capacity(n);
            let mut encodings = Vec::with_capacity(n);
            let mut missing = Vec::with_capacity(n);
            
            for val in dv.data {
                if val.is_nan() {
                    // 检查是否是 R 的 NA（特殊的 NaN）
                    let bits = val.to_bits();
                    let low_word = (bits & 0xFFFFFFFF) as u32;
                    if low_word == 1954 {
                        // R NA
                        data.push(String::new());
                        encodings.push(StringEncoding::None);
                        missing.push(true);
                    } else {
                        data.push("NaN".to_string());
                        encodings.push(StringEncoding::Ascii);
                        missing.push(false);
                    }
                } else if val.is_infinite() {
                    if val > 0.0 {
                        data.push("Inf".to_string());
                    } else {
                        data.push("-Inf".to_string());
                    }
                    encodings.push(StringEncoding::Ascii);
                    missing.push(false);
                } else {
                    data.push(format!("{}", val));
                    encodings.push(StringEncoding::Ascii);
                    missing.push(false);
                }
            }
            
            StringVector {
                data,
                encodings,
                missing,
                attributes: Default::default(),
            }
        }
        _ => {
            return Err(RdsError::ParseError {
                context: "deferred_string".to_string(),
                message: format!("Unsupported content type in deferred_string: {:?}", contents),
            });
        }
    };
    
    // 解析元数据
    let meta_header = read_header(reader)?;
    if meta_header[3] != SEXPType::Int as u8 {
        return Err(RdsError::ParseError {
            context: "deferred_string".to_string(),
            message: format!("Expected integer vector for metadata, got type {}", meta_header[3]),
        });
    }
    let _metadata = parse_integer_body(reader)?;
    
    // 读取终止符
    let terminator = read_header(reader)?;
    if terminator[3] != SEXPType::NilValue as u8 && terminator[3] != SEXPType::Nil as u8 {
        return Err(RdsError::ParseError {
            context: "deferred_string".to_string(),
            message: format!("Expected NULL terminator, got type {}", terminator[3]),
        });
    }
    
    Ok(RObject::StringVector(output))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compact_intseq_expansion() {
        // 测试紧凑整数序列的展开逻辑
        let len = 5usize;
        let start = 1.0f64;
        let step = 1.0f64;
        
        let mut data = Vec::with_capacity(len);
        let mut current = start;
        for _ in 0..len {
            data.push(current as i32);
            current += step;
        }
        
        assert_eq!(data, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_compact_intseq_step2() {
        let len = 3usize;
        let start = 0.0f64;
        let step = 2.0f64;
        
        let mut data = Vec::with_capacity(len);
        let mut current = start;
        for _ in 0..len {
            data.push(current as i32);
            current += step;
        }
        
        assert_eq!(data, vec![0, 2, 4]);
    }
    
    #[test]
    fn test_deferred_string_int_conversion() {
        // 测试整数到字符串的转换逻辑
        let int_data = vec![1, 2, i32::MIN, 42];
        let mut data = Vec::new();
        let mut missing = Vec::new();
        
        for val in int_data {
            if val == i32::MIN {
                data.push(String::new());
                missing.push(true);
            } else {
                data.push(val.to_string());
                missing.push(false);
            }
        }
        
        assert_eq!(data, vec!["1", "2", "", "42"]);
        assert_eq!(missing, vec![false, false, true, false]);
    }
}
