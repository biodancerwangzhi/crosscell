//! 对象解析分发
//!
//! 根据 SEXP 类型分发到具体的解析器。
//! 对应 rds2cpp 的 parse_object 函数。

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::r_object::RObject;
use super::utils::read_header;
use super::header::{ParsedHeader, has_attributes, has_tag};
use super::shared_info::SharedParseInfo;

// 导入各类型解析器
use super::atomic::{
    parse_integer_body, parse_logical_body, parse_double_body,
    parse_raw_body, parse_complex_body, parse_string_body_with_ref,
};
use super::list::parse_list_body;
use super::pairlist::parse_pairlist_body;
use super::symbol::parse_symbol_body;
use super::environment::{
    parse_new_environment_body, parse_global_environment_body,
    parse_base_environment_body, parse_empty_environment_body,
};
use super::external_pointer::parse_external_pointer_body;
use super::s4::parse_s4_body;
use super::language::parse_language_body;
use super::expression::parse_expression_body;
use super::builtin::parse_builtin_body;
use super::altrep::parse_altrep_body;
use super::attributes::parse_attributes;

/// 解析任意 R 对象
///
/// 从 reader 读取一个 R 对象，根据头部的 SEXP 类型分发到具体的解析器。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
/// * `shared` - 共享解析信息，用于处理引用和循环结构
///
/// # 返回
/// 解析后的 RObject
///
/// # 错误
/// 如果遇到不支持的类型或解析失败返回错误
pub fn parse_object<R: Read>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
) -> Result<RObject> {
    let header_bytes = read_header(reader)?;
    let header = ParsedHeader::from_raw(header_bytes)
        .ok_or_else(|| RdsError::UnsupportedType(header_bytes[3]))?;
    
    parse_object_with_header(reader, &header, shared)
}


/// 使用已解析的头部解析对象
///
/// 当头部已经被读取时使用此函数。
pub fn parse_object_with_header<R: Read>(
    reader: &mut R,
    header: &ParsedHeader,
    shared: &mut SharedParseInfo,
) -> Result<RObject> {
    // 创建递归解析闭包
    let mut parse_fn = |r: &mut R, s: &mut SharedParseInfo| -> Result<RObject> {
        parse_object(r, s)
    };
    
    let obj = match header.sexp_type {
        // ========================================
        // 特殊值
        // ========================================
        SEXPType::Nil | SEXPType::NilValue => RObject::Null,
        
        SEXPType::MissingArg | SEXPType::UnboundValue => RObject::Null,
        
        // ========================================
        // 引用类型
        // ========================================
        SEXPType::Ref => {
            return shared.resolve_reference(&header.raw);
        }
        
        // ========================================
        // 原子向量
        // ========================================
        SEXPType::Int => {
            let mut vec = parse_integer_body(reader)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::IntegerVector(vec)
        }
        
        SEXPType::Lgl => {
            let mut vec = parse_logical_body(reader)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::LogicalVector(vec)
        }
        
        SEXPType::Real => {
            let mut vec = parse_double_body(reader)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::DoubleVector(vec)
        }
        
        SEXPType::Raw => {
            let mut vec = parse_raw_body(reader)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::RawVector(vec)
        }
        
        SEXPType::Cplx => {
            let mut vec = parse_complex_body(reader)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::ComplexVector(vec)
        }
        
        SEXPType::Str => {
            let mut vec = parse_string_body_with_ref(reader, shared)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::StringVector(vec)
        }
        
        // ========================================
        // 复合类型
        // ========================================
        SEXPType::Vec => {
            let mut vec = parse_list_body(reader, shared, &mut parse_fn)?;
            if has_attributes(&header.raw) {
                vec.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::GenericVector(vec)
        }
        
        SEXPType::List | SEXPType::AttrList => {
            // 注意：pairlist 的属性在 parse_pairlist_body 内部的 recursive_parse 中处理
            // 不需要在这里再次解析属性
            let pairlist = parse_pairlist_body(reader, header, shared, &mut parse_fn)?;
            RObject::PairList(pairlist)
        }
        
        SEXPType::Expr => {
            let mut expr = parse_expression_body(reader, shared, &mut parse_fn)?;
            if has_attributes(&header.raw) {
                expr.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::ExpressionVector(expr)
        }
        
        // ========================================
        // 符号
        // ========================================
        SEXPType::Sym => {
            let sym_idx = parse_symbol_body(reader, shared)?;
            RObject::SymbolIndex(sym_idx)
        }
        
        // ========================================
        // 环境
        // ========================================
        SEXPType::Env => {
            let env_idx = parse_new_environment_body(reader, header, shared, &mut parse_fn)?;
            RObject::EnvironmentIndex(env_idx)
        }
        
        SEXPType::GlobalEnv => {
            RObject::EnvironmentIndex(parse_global_environment_body())
        }
        
        SEXPType::BaseEnv => {
            RObject::EnvironmentIndex(parse_base_environment_body())
        }
        
        SEXPType::EmptyEnv => {
            RObject::EnvironmentIndex(parse_empty_environment_body())
        }
        
        // ========================================
        // 外部指针
        // ========================================
        SEXPType::Extptr => {
            let ptr_idx = parse_external_pointer_body(reader, header, shared, &mut parse_fn)?;
            RObject::ExternalPointerIndex(ptr_idx)
        }
        
        // ========================================
        // S4 对象
        // ========================================
        SEXPType::S4 => {
            let s4 = parse_s4_body(reader, header, shared, &mut parse_fn)?;
            RObject::S4Object(s4)
        }
        
        // ========================================
        // 语言对象
        // ========================================
        SEXPType::Lang | SEXPType::AttrLang => {
            let mut lang = parse_language_body(reader, header, shared, &mut parse_fn)?;
            if has_attributes(&header.raw) && header.sexp_type != SEXPType::AttrLang {
                lang.attributes = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::LanguageObject(lang)
        }
        
        // ========================================
        // 内置函数
        // ========================================
        SEXPType::Builtin | SEXPType::Special => {
            let builtin = parse_builtin_body(reader)?;
            RObject::BuiltInFunction(builtin)
        }
        
        // ========================================
        // ALTREP
        // ========================================
        SEXPType::Altrep => {
            parse_altrep_body(reader, shared, &mut parse_fn)?
        }
        
        // ========================================
        // 闭包 (Closure)
        // R 序列化顺序: ATTRIB (if any), CLOENV, FORMALS, BODY
        // 注意: CLOSXP 的 HAS_TAG 标志总是设置，但不写入 TAG
        // ========================================
        SEXPType::Clo => {
            // R 序列化顺序: ATTRIB (if hasattr), CLOENV (always, hastag=true), FORMALS, BODY
            // 闭包有属性时先解析属性
            if has_attributes(&header.raw) {
                let _attrs = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            // 解析 environment (CLOENV = TAG, 对于 CLOSXP hastag 总是 true)
            let _env = parse_fn(reader, shared)?;
            // 解析 formals (参数列表 = CAR)
            let _formals = parse_fn(reader, shared)?;
            // 解析 body (函数体 = CDR, 通过 tailcall 写入)
            let _body = parse_fn(reader, shared)?;
            RObject::Null
        }
        
        // ========================================
        // Promise
        // 结构: value, expression, environment
        // ========================================
        SEXPType::Prom => {
            // R 序列化顺序: ATTRIB (if hasattr), TAG/env (if hastag), CAR/value, CDR/expression
            if has_attributes(&header.raw) {
                let _attrs = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            // TAG = environment (只在 hastag 时读取)
            if has_tag(&header.raw) {
                let _env = parse_fn(reader, shared)?;
            }
            // CAR = value
            let _value = parse_fn(reader, shared)?;
            // CDR = expression (通过 tailcall 写入)
            let _expr = parse_fn(reader, shared)?;
            RObject::Null
        }
        
        // ========================================
        // Dot-dot-dot (...)
        // 结构类似 pairlist
        // ========================================
        SEXPType::Dot => {
            // ... 参数是一个特殊的 pairlist
            let _pairlist = parse_pairlist_body(reader, header, shared, &mut parse_fn)?;
            RObject::Null
        }
        
        // ========================================
        // 字节码
        // 结构: reps_len, bc1 (code + consts)
        // ========================================
        SEXPType::Bcode => {
            // R 源码: ReadBC 读取 reps_len + ReadBC1, 然后在 ReadItem_Recursive
            // 的 switch 之后统一读取属性: SET_ATTRIB(s, hasattr ? ReadItem(...) : R_NilValue)
            let mut len_buf = [0u8; 4];
            reader.read_exact(&mut len_buf)?;
            let reps_len = i32::from_be_bytes(len_buf) as usize;
            
            // 解析 bc1: code + consts
            parse_bc1(reader, shared, &mut parse_fn, reps_len)?;
            
            // 属性在字节码内容之后读取（与 R 源码 ReadItem_Recursive 一致）
            if has_attributes(&header.raw) {
                let _attrs = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::Null
        }
        
        SEXPType::BcRepDef | SEXPType::BcRepRef => {
            // 字节码引用 - 这些不应该在顶层出现
            RObject::Null
        }
        
        // ========================================
        // 弱引用
        // 结构: key, value, finalizer, next
        // 注意：R 会将弱引用添加到引用表中
        // ========================================
        SEXPType::Weakref => {
            // R 源码: WEAKREFSXP 写入时只 HashAdd，不写入任何数据
            // 读取时只 R_MakeWeakRef + AddReadRef，不读取任何数据
            let ptr_idx = shared.request_external_pointer();
            // 属性在 switch 之后统一处理
            if has_attributes(&header.raw) {
                let _attrs = parse_attributes(reader, shared, &mut parse_fn)?;
            }
            RObject::ExternalPointerIndex(crate::rds::r_object::ExternalPointerIndex { index: ptr_idx })
        }
        
        // ========================================
        // 命名空间和包
        // R 序列化格式: 头部 + 字符串向量（包名/命名空间名）
        // 注意：R 会将这些添加到引用表中，所以我们也需要添加
        // ========================================
        SEXPType::Namespace | SEXPType::Package => {
            // R 源码: InStringVec 格式 = placeholder(0) + length + ReadItem*length
            // 注意：这不是标准的 ReadItem，而是特殊的 InStringVec 格式
            let env_idx = shared.request_environment();
            // 读取 InStringVec: placeholder integer (应为 0)
            let mut placeholder_buf = [0u8; 4];
            reader.read_exact(&mut placeholder_buf)?;
            // 读取字符串数量
            let mut len_buf = [0u8; 4];
            reader.read_exact(&mut len_buf)?;
            let n = i32::from_be_bytes(len_buf) as usize;
            // 读取每个字符串元素 (通过 ReadItem)
            for _ in 0..n {
                let _item = parse_fn(reader, shared)?;
            }
            RObject::EnvironmentIndex(crate::rds::r_object::EnvironmentIndex {
                index: env_idx,
                env_type: header.sexp_type,
            })
        }
        
        SEXPType::BaseNamespace => {
            // BaseNamespace 没有额外数据，但也需要添加到引用表
            // 实际上 BaseNamespace 在 R 中不添加到引用表
            RObject::Null
        }
        
        SEXPType::Persist => {
            // R 源码: PERSISTSXP 也使用 InStringVec 格式
            let env_idx = shared.request_environment();
            // 读取 InStringVec: placeholder integer (应为 0)
            let mut placeholder_buf = [0u8; 4];
            reader.read_exact(&mut placeholder_buf)?;
            // 读取字符串数量
            let mut len_buf = [0u8; 4];
            reader.read_exact(&mut len_buf)?;
            let n = i32::from_be_bytes(len_buf) as usize;
            // 读取每个字符串元素
            for _ in 0..n {
                let _item = parse_fn(reader, shared)?;
            }
            RObject::EnvironmentIndex(crate::rds::r_object::EnvironmentIndex {
                index: env_idx,
                env_type: header.sexp_type,
            })
        }
        
        SEXPType::ClassRef | SEXPType::GenericRef => {
            // 类引用 - 简化处理为 Null
            RObject::Null
        }
        
        // CHAR 类型 - 在某些情况下可能出现在顶层（如 ALTREP）
        // 将其解析为单元素字符串向量
        SEXPType::Char => {
            use super::single_string::parse_single_string_with_header;
            let string_info = parse_single_string_with_header(reader, header)?;
            let mut vec = crate::rds::r_object::StringVector::default();
            vec.data.push(string_info.value);
            vec.missing.push(string_info.missing);
            vec.encodings.push(string_info.encoding);
            RObject::StringVector(vec)
        }
        
        // 其他未知类型
        SEXPType::Any => {
            RObject::Null
        }
    };
    
    Ok(obj)
}

/// 解析字节码 bc1 结构
fn parse_bc1<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
    reps_len: usize,
) -> Result<()>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    // 解析 code
    let _code = parse_fn(reader, shared)?;
    // 解析 consts
    parse_bc_consts(reader, shared, parse_fn, reps_len)?;
    Ok(())
}

/// 解析字节码常量
/// 
/// R 的 ReadBCConsts 格式：
/// - 先读取一个 i32 作为类型指示符
/// - 对于 BCODESXP: 递归调用 ReadBC1
/// - 对于 LANGSXP/LISTSXP/BCREPDEF/BCREPREF/ATTRLANGSXP/ATTRLISTSXP: 调用 ReadBCLang
/// - 对于其他类型（default）: 调用 ReadItem（读取完整对象，包括头部）
/// 
/// 关键点：类型指示符只是一个提示，对于 default 分支，实际对象（包括完整头部）紧随其后
fn parse_bc_consts<R: Read, F>(
    reader: &mut R,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
    reps_len: usize,
) -> Result<()>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    // 读取常量数量
    let mut len_buf = [0u8; 4];
    reader.read_exact(&mut len_buf)?;
    let n = i32::from_be_bytes(len_buf) as usize;
    
    for _i in 0..n {
        // 读取类型指示符（作为 i32）
        // 注意：这只是一个提示，不是实际的对象头部
        let mut type_buf = [0u8; 4];
        reader.read_exact(&mut type_buf)?;
        let type_int = i32::from_be_bytes(type_buf);
        let type_byte = type_int as u8;
        let const_type = SEXPType::from_u8(type_byte);

        match const_type {
            Some(SEXPType::Bcode) => {
                // 递归解析嵌套字节码
                // BCODESXP 特殊处理：直接调用 ReadBC1，不读取额外头部
                parse_bc1(reader, shared, parse_fn, reps_len)?;
            }
            Some(SEXPType::Lang) | Some(SEXPType::List) |
            Some(SEXPType::BcRepDef) | Some(SEXPType::BcRepRef) |
            Some(SEXPType::AttrLang) | Some(SEXPType::AttrList) => {
                // 解析字节码语言对象
                // 这些类型有特殊的序列化格式，不是标准对象
                parse_bc_lang(reader, type_byte, shared, parse_fn, reps_len)?;
            }

            _ => {
                // default 分支：调用 ReadItem 读取完整对象
                // 关键：类型指示符只是提示，实际对象（包括完整 4 字节头部）紧随其后
                // R 源码: SET_VECTOR_ELT(ans, i, ReadItem(ref_table, stream));
                let _obj = parse_fn(reader, shared)?;
            }
        }
    }
    Ok(())
}

/// 解析字节码语言对象
fn parse_bc_lang<R: Read, F>(
    reader: &mut R,
    type_byte: u8,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
    _reps_len: usize,
) -> Result<RObject>
where F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject> {
    let sexp_type = SEXPType::from_u8(type_byte);
    match sexp_type {
        Some(SEXPType::BcRepRef) => {
            // 读取引用位置
            let mut pos_buf = [0u8; 4];
            reader.read_exact(&mut pos_buf)?;
            let _pos = i32::from_be_bytes(pos_buf) as usize;
            Ok(RObject::Null)
        }
        Some(SEXPType::BcRepDef) => {
            // 读取定义位置
            let mut pos_buf = [0u8; 4];
            reader.read_exact(&mut pos_buf)?;
            let _pos = i32::from_be_bytes(pos_buf) as usize;
            
            // 读取实际类型
            let mut type_buf = [0u8; 4];
            reader.read_exact(&mut type_buf)?;
            let actual_type = i32::from_be_bytes(type_buf) as u8;
            
            // 确定是否有属性
            let has_attr = actual_type == SEXPType::AttrLang as u8 || 
                          actual_type == SEXPType::AttrList as u8;
            
            // 如果有属性，读取属性
            if has_attr {
                let _attr = parse_fn(reader, shared)?;
            }
            
            // 解析 tag
            let _tag = parse_fn(reader, shared)?;
            
            // 解析 car (递归调用 bc_lang)
            let mut car_type_buf = [0u8; 4];
            reader.read_exact(&mut car_type_buf)?;
            let car_type = i32::from_be_bytes(car_type_buf) as u8;
            let _car = parse_bc_lang(reader, car_type, shared, parse_fn, _reps_len)?;
            
            // 解析 cdr (递归调用 bc_lang)
            let mut cdr_type_buf = [0u8; 4];
            reader.read_exact(&mut cdr_type_buf)?;
            let cdr_type = i32::from_be_bytes(cdr_type_buf) as u8;
            let _cdr = parse_bc_lang(reader, cdr_type, shared, parse_fn, _reps_len)?;
            
            Ok(RObject::Null)
        }
        Some(SEXPType::Lang) | Some(SEXPType::AttrLang) |
        Some(SEXPType::List) | Some(SEXPType::AttrList) => {
            let has_attr = sexp_type == Some(SEXPType::AttrLang) || sexp_type == Some(SEXPType::AttrList);
            
            // 如果有属性，读取属性
            if has_attr {
                let _attr = parse_fn(reader, shared)?;
            }
            
            // 解析 tag
            let _tag = parse_fn(reader, shared)?;
            
            // 解析 car (递归调用 bc_lang)
            let mut car_type_buf = [0u8; 4];
            reader.read_exact(&mut car_type_buf)?;
            let car_type = i32::from_be_bytes(car_type_buf) as u8;
            let _car = parse_bc_lang(reader, car_type, shared, parse_fn, _reps_len)?;
            
            // 解析 cdr (递归调用 bc_lang)
            let mut cdr_type_buf = [0u8; 4];
            reader.read_exact(&mut cdr_type_buf)?;
            let cdr_type = i32::from_be_bytes(cdr_type_buf) as u8;
            let _cdr = parse_bc_lang(reader, cdr_type, shared, parse_fn, _reps_len)?;
            
            Ok(RObject::Null)
        }
        Some(SEXPType::Nil) => {
            // NILSXP (0) 作为类型指示符时，是 WriteBCLang 的 pad 字节
            // WriteBCLang: OutInteger(stream, 0); /* pad */ WriteItem(s, ...);
            // 所以后面跟着一个完整对象，需要调用 ReadItem
            parse_fn(reader, shared)
        }
        Some(SEXPType::NilValue) => {
            // NILVALUE_SXP (254) 不应该出现在 ReadBCLang 的类型指示符中
            // 但如果出现，按 default 分支处理：调用 ReadItem
            parse_fn(reader, shared)
        }
        _ => {
            // 其他类型作为普通对象解析
            // 注意：这里的 type_byte 只是类型指示符，实际对象跟在后面
            parse_fn(reader, shared)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::string_encoding::StringEncoding;
    use crate::rds::parse::header::encoding_flags;

    fn make_header(sexp_type: SEXPType, flags: u8) -> [u8; 4] {
        [0x00, 0x00, flags, sexp_type as u8]
    }

    fn make_length(len: i32) -> [u8; 4] {
        len.to_be_bytes()
    }

    fn make_i32(val: i32) -> [u8; 4] {
        val.to_be_bytes()
    }

    fn make_f64(val: f64) -> [u8; 8] {
        val.to_be_bytes()
    }

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
    // 基础类型测试
    // ========================================

    #[test]
    fn test_parse_null() {
        let data = make_header(SEXPType::NilValue, 0).to_vec();
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        assert!(matches!(obj, RObject::Null));
    }

    #[test]
    fn test_parse_integer_vector() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Int, 0));
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&make_i32(1));
        data.extend_from_slice(&make_i32(2));
        data.extend_from_slice(&make_i32(3));
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::IntegerVector(vec) = obj {
            assert_eq!(vec.data, vec![1, 2, 3]);
        } else {
            panic!("Expected IntegerVector");
        }
    }

    #[test]
    fn test_parse_double_vector() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Real, 0));
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_f64(1.5));
        data.extend_from_slice(&make_f64(2.5));
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::DoubleVector(vec) = obj {
            assert_eq!(vec.data, vec![1.5, 2.5]);
        } else {
            panic!("Expected DoubleVector");
        }
    }

    #[test]
    fn test_parse_logical_vector() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Lgl, 0));
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_i32(1));  // TRUE
        data.extend_from_slice(&make_i32(0));  // FALSE
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::LogicalVector(vec) = obj {
            assert_eq!(vec.data, vec![1, 0]);
        } else {
            panic!("Expected LogicalVector");
        }
    }

    #[test]
    fn test_parse_raw_vector() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Raw, 0));
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(&[0x01, 0x02, 0x03]);
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::RawVector(vec) = obj {
            assert_eq!(vec.data, vec![0x01, 0x02, 0x03]);
        } else {
            panic!("Expected RawVector");
        }
    }

    #[test]
    fn test_parse_string_vector() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Str, 0));
        data.extend_from_slice(&make_length(1));
        data.extend_from_slice(&make_char_header(StringEncoding::Utf8));
        data.extend_from_slice(&make_length(5));
        data.extend_from_slice(b"hello");
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::StringVector(vec) = obj {
            assert_eq!(vec.data, vec!["hello"]);
        } else {
            panic!("Expected StringVector");
        }
    }

    #[test]
    fn test_parse_empty_list() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Vec, 0));
        data.extend_from_slice(&make_length(0));
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::GenericVector(vec) = obj {
            assert!(vec.data.is_empty());
        } else {
            panic!("Expected GenericVector");
        }
    }

    #[test]
    fn test_parse_list_with_nulls() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Vec, 0));
        data.extend_from_slice(&make_length(2));
        data.extend_from_slice(&make_header(SEXPType::NilValue, 0));
        data.extend_from_slice(&make_header(SEXPType::NilValue, 0));
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::GenericVector(vec) = obj {
            assert_eq!(vec.data.len(), 2);
            assert!(matches!(vec.data[0], RObject::Null));
            assert!(matches!(vec.data[1], RObject::Null));
        } else {
            panic!("Expected GenericVector");
        }
    }

    #[test]
    fn test_parse_builtin() {
        let mut data = Vec::new();
        data.extend_from_slice(&make_header(SEXPType::Builtin, 0));
        data.extend_from_slice(&make_length(3));
        data.extend_from_slice(b"sum");
        
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::BuiltInFunction(f) = obj {
            assert_eq!(f.name, "sum");
        } else {
            panic!("Expected BuiltInFunction");
        }
    }

    #[test]
    fn test_parse_global_env() {
        let data = make_header(SEXPType::GlobalEnv, 0).to_vec();
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::EnvironmentIndex(idx) = obj {
            assert_eq!(idx.env_type, SEXPType::GlobalEnv);
        } else {
            panic!("Expected EnvironmentIndex");
        }
    }

    #[test]
    fn test_parse_empty_env() {
        let data = make_header(SEXPType::EmptyEnv, 0).to_vec();
        let mut cursor = Cursor::new(data);
        let mut shared = SharedParseInfo::new();
        
        let obj = parse_object(&mut cursor, &mut shared).unwrap();
        if let RObject::EnvironmentIndex(idx) = obj {
            assert_eq!(idx.env_type, SEXPType::EmptyEnv);
        } else {
            panic!("Expected EnvironmentIndex");
        }
    }
}