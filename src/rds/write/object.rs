//! 对象写入分发
//!
//! 根据 RObject 类型分发到具体的写入器。
//! 对应 rds2cpp 的 write_object 函数。

use super::atomic::{
    write_complex_body, write_double_body, write_integer_body, write_logical_body, write_raw_body,
    write_string_body,
};
use super::attributes::write_attributes;
use super::builtin::write_builtin_body;
use super::expression::write_expression_body;
use super::language::write_language;
use super::list::write_list_body;
use super::pairlist::write_pairlist;
use super::s4::write_s4;
use super::shared_info::SharedWriteInfo;
use super::utils::{write_header, write_header_with_attributes};
use crate::rds::error::Result;
use crate::rds::r_object::RObject;
use crate::rds::sexp_type::SEXPType;
use std::io::Write;

/// 写入任意 R 对象
///
/// 根据 RObject 类型分发到具体的写入器。
///
/// # 参数
/// * `object` - 要写入的 R 对象
/// * `writer` - 实现 Write trait 的目标
/// * `shared` - 共享写入信息，用于处理引用和去重
///
/// # 错误
/// 如果遇到不支持的类型或写入失败返回错误
pub fn write_object<W: Write>(
    object: &RObject,
    writer: &mut W,
    shared: &mut SharedWriteInfo,
) -> Result<()> {
    // 创建递归写入闭包
    let mut write_fn = |obj: &RObject, w: &mut W, s: &mut SharedWriteInfo| -> Result<()> {
        write_object(obj, w, s)
    };

    match object {
        // ========================================
        // 特殊值
        // ========================================
        RObject::Null => {
            write_header(SEXPType::NilValue, writer)?;
        }

        // ========================================
        // 原子向量
        // ========================================
        RObject::IntegerVector(vec) => {
            write_header_with_attributes(SEXPType::Int, &vec.attributes, writer)?;
            write_integer_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::LogicalVector(vec) => {
            write_header_with_attributes(SEXPType::Lgl, &vec.attributes, writer)?;
            write_logical_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::DoubleVector(vec) => {
            write_header_with_attributes(SEXPType::Real, &vec.attributes, writer)?;
            write_double_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::RawVector(vec) => {
            write_header_with_attributes(SEXPType::Raw, &vec.attributes, writer)?;
            write_raw_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::ComplexVector(vec) => {
            write_header_with_attributes(SEXPType::Cplx, &vec.attributes, writer)?;
            write_complex_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::StringVector(vec) => {
            write_header_with_attributes(SEXPType::Str, &vec.attributes, writer)?;
            write_string_body(vec, writer)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        // ========================================
        // 复合类型
        // ========================================
        RObject::GenericVector(vec) => {
            write_header_with_attributes(SEXPType::Vec, &vec.attributes, writer)?;
            write_list_body(vec, writer, shared, &mut write_fn)?;
            if !vec.attributes.is_empty() {
                write_attributes(&vec.attributes, writer, shared, &mut write_fn)?;
            }
        }

        RObject::PairList(pairlist) => {
            // PairList 头部在 write_pairlist 内部处理
            write_pairlist(pairlist, writer, shared, &mut write_fn)?;
        }

        RObject::ExpressionVector(expr) => {
            write_header_with_attributes(SEXPType::Expr, &expr.attributes, writer)?;
            write_expression_body(expr, writer, shared, &mut write_fn)?;
            if !expr.attributes.is_empty() {
                write_attributes(&expr.attributes, writer, shared, &mut write_fn)?;
            }
        }

        // ========================================
        // 符号
        // ========================================
        RObject::SymbolIndex(sym_idx) => {
            shared.write_symbol_ref(sym_idx, writer)?;
        }

        // ========================================
        // 环境
        // ========================================
        RObject::EnvironmentIndex(env_idx) => {
            shared.write_environment(env_idx, writer)?;
        }

        // ========================================
        // 外部指针
        // ========================================
        RObject::ExternalPointerIndex(ptr_idx) => {
            shared.write_external_pointer(ptr_idx, writer)?;
        }

        // ========================================
        // S4 对象
        // ========================================
        RObject::S4Object(s4) => {
            // S4 头部 - 参考 rds2cpp write_s4.hpp
            // byte 0: 0
            // byte 1: 0x1 (see logic in parse_s4)
            // byte 2: 0x1 | 0x2 (is_object | has_attributes)
            // byte 3: S4 type
            let header = [0, 0x01, 0x01 | 0x02, SEXPType::S4 as u8];
            writer.write_all(&header)?;
            write_s4(s4, writer, shared, &mut write_fn)?;
        }

        // ========================================
        // 语言对象
        // ========================================
        RObject::LanguageObject(lang) => {
            // LANG 头部
            let has_attrs = !lang.attributes.is_empty();
            let flags = if has_attrs { 0x02 } else { 0x00 };
            let header = [0, 0, flags, SEXPType::Lang as u8];
            writer.write_all(&header)?;
            write_language(lang, writer, shared, &mut write_fn)?;
        }

        // ========================================
        // 内置函数
        // ========================================
        RObject::BuiltInFunction(func) => {
            write_header(SEXPType::Builtin, writer)?;
            write_builtin_body(func, writer)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::environment::Environment;
    use crate::rds::external_pointer::ExternalPointer;
    use crate::rds::r_object::{
        Attributes, BuiltInFunction, ComplexVector, DoubleVector, EnvironmentIndex, GenericVector,
        IntegerVector, LanguageObject, LogicalVector, PairList, RawVector, S4Object, StringVector,
        SymbolIndex,
    };
    use crate::rds::string_encoding::StringEncoding;
    use crate::rds::symbol::Symbol;
    use num_complex::Complex64;
    use std::io::Cursor;

    fn make_shared<'a>(
        symbols: &'a [Symbol],
        environments: &'a [Environment],
        external_pointers: &'a [ExternalPointer],
    ) -> SharedWriteInfo<'a> {
        SharedWriteInfo::new(symbols, environments, external_pointers)
    }

    #[test]
    fn test_write_null() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        write_object(&RObject::Null, &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data.len(), 4);
        assert_eq!(data[3], SEXPType::NilValue as u8);
    }

    #[test]
    fn test_write_integer_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = IntegerVector {
            data: vec![1, 2, 3],
            attributes: Default::default(),
        };
        write_object(&RObject::IntegerVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Int as u8);
    }

    #[test]
    fn test_write_logical_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = LogicalVector {
            data: vec![1, 0, i32::MIN],
            attributes: Default::default(),
        };
        write_object(&RObject::LogicalVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Lgl as u8);
    }

    #[test]
    fn test_write_double_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = DoubleVector {
            data: vec![1.5, 2.5, 3.5],
            attributes: Default::default(),
        };
        write_object(&RObject::DoubleVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Real as u8);
    }

    #[test]
    fn test_write_raw_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = RawVector {
            data: vec![0x01, 0x02, 0x03],
            attributes: Default::default(),
        };
        write_object(&RObject::RawVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Raw as u8);
    }

    #[test]
    fn test_write_complex_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = ComplexVector {
            data: vec![Complex64::new(1.0, 2.0)],
            attributes: Default::default(),
        };
        write_object(&RObject::ComplexVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Cplx as u8);
    }

    #[test]
    fn test_write_string_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = StringVector {
            data: vec!["hello".to_string(), "world".to_string()],
            encodings: vec![StringEncoding::Utf8, StringEncoding::Utf8],
            missing: vec![false, false],
            attributes: Default::default(),
        };
        write_object(&RObject::StringVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Str as u8);
    }

    #[test]
    fn test_write_generic_vector() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let vec = GenericVector {
            data: vec![RObject::Null, RObject::Null],
            attributes: Default::default(),
        };
        write_object(&RObject::GenericVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Vec as u8);
    }

    #[test]
    fn test_write_pairlist() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let pairlist = PairList {
            data: vec![RObject::Null],
            has_tag: vec![true],
            tag_names: vec!["x".to_string()],
            tag_encodings: vec![StringEncoding::Utf8],
            attributes: Default::default(),
        };
        write_object(&RObject::PairList(pairlist), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        // PairList 以 LIST 头部开始
        assert_eq!(data[3], SEXPType::List as u8);
    }

    #[test]
    fn test_write_builtin() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let func = BuiltInFunction {
            name: "sum".to_string(),
        };
        write_object(&RObject::BuiltInFunction(func), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Builtin as u8);
    }

    #[test]
    fn test_write_language_object() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let lang = LanguageObject {
            function_name: "print".to_string(),
            function_encoding: StringEncoding::Utf8,
            argument_values: vec![],
            argument_names: vec![],
            argument_has_name: vec![],
            argument_encodings: vec![],
            attributes: Default::default(),
        };
        write_object(&RObject::LanguageObject(lang), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::Lang as u8);
    }

    #[test]
    fn test_write_s4_object() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let s4 = S4Object {
            class_name: "MyClass".to_string(),
            class_encoding: StringEncoding::Utf8,
            package_name: "mypackage".to_string(),
            package_encoding: StringEncoding::Utf8,
            attributes: Default::default(),
        };
        write_object(&RObject::S4Object(s4), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::S4 as u8);
    }

    #[test]
    fn test_write_symbol_ref() {
        let symbols = vec![Symbol::new("x".to_string(), StringEncoding::Utf8)];
        let (e, p) = (vec![], vec![]);
        let mut shared = make_shared(&symbols, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let sym_idx = SymbolIndex { index: 0 };
        write_object(&RObject::SymbolIndex(sym_idx), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        // 第一次写入符号，应该是 SYM 头部
        assert_eq!(data[3], SEXPType::Sym as u8);
    }

    #[test]
    fn test_write_global_env() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let env_idx = EnvironmentIndex {
            index: usize::MAX,
            env_type: SEXPType::GlobalEnv,
        };
        write_object(
            &RObject::EnvironmentIndex(env_idx),
            &mut buffer,
            &mut shared,
        )
        .unwrap();

        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::GlobalEnv as u8);
    }

    #[test]
    fn test_write_vector_with_attributes() {
        let (s, e, p) = (vec![], vec![], vec![]);
        let mut shared = make_shared(&s, &e, &p);
        let mut buffer = Cursor::new(Vec::new());

        let mut attrs = Attributes::new();
        attrs.add(
            "names".to_string(),
            RObject::StringVector(StringVector {
                data: vec!["a".to_string(), "b".to_string()],
                encodings: vec![StringEncoding::Utf8, StringEncoding::Utf8],
                missing: vec![false, false],
                attributes: Default::default(),
            }),
            StringEncoding::Utf8,
        );

        let vec = IntegerVector {
            data: vec![1, 2],
            attributes: attrs,
        };
        write_object(&RObject::IntegerVector(vec), &mut buffer, &mut shared).unwrap();

        let data = buffer.into_inner();
        // 头部应该有属性标志
        assert_eq!(data[2] & 0x02, 0x02);
    }
}
