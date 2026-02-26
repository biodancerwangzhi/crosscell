//! 写入工具函数
//!
//! 提供底层写入工具，包括字节序转换、长度写入等。
//! 对应 rds2cpp 的写入工具函数。

use crate::rds::error::Result;
use crate::rds::r_object::{Attributes, RObject};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use std::io::Write;

/// 4 字节头部类型
pub type Header = [u8; 4];

/// 将 i32 转换为大端序字节
///
/// # 参数
/// * `value` - 要转换的 i32 值
///
/// # 返回
/// 4 字节大端序数组
#[inline]
pub fn i32_to_be_bytes(value: i32) -> [u8; 4] {
    value.to_be_bytes()
}

/// 将 f64 转换为大端序字节
///
/// # 参数
/// * `value` - 要转换的 f64 值
///
/// # 返回
/// 8 字节大端序数组
#[inline]
pub fn f64_to_be_bytes(value: f64) -> [u8; 8] {
    value.to_be_bytes()
}

/// 注入大端序 i32 到 buffer
///
/// # 参数
/// * `value` - 要写入的 i32 值
/// * `buffer` - 目标缓冲区
#[inline]
pub fn inject_i32_be(value: i32, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&value.to_be_bytes());
}

/// 注入大端序 f64 到 buffer
///
/// # 参数
/// * `value` - 要写入的 f64 值
/// * `buffer` - 目标缓冲区
#[inline]
pub fn inject_f64_be(value: f64, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&value.to_be_bytes());
}

/// 写入大端序 i32 到 writer
///
/// # 参数
/// * `value` - 要写入的 i32 值
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_i32<W: Write>(value: i32, writer: &mut W) -> Result<()> {
    writer.write_all(&value.to_be_bytes())?;
    Ok(())
}

/// 写入大端序 f64 到 writer
///
/// # 参数
/// * `value` - 要写入的 f64 值
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_f64<W: Write>(value: f64, writer: &mut W) -> Result<()> {
    writer.write_all(&value.to_be_bytes())?;
    Ok(())
}

/// 注入长度到 buffer
///
/// RDS 格式支持两种长度格式：
/// - 标准格式：4 字节大端序无符号整数（最大 2^31-1）
/// - 大长度格式：当长度 >= 2^31 时，先写 -1，再写 8 字节实际长度
///
/// # 参数
/// * `value` - 要写入的长度值
/// * `buffer` - 目标缓冲区
pub fn inject_length(value: usize, buffer: &mut Vec<u8>) {
    if value < 0x80000000 {
        // 标准格式
        inject_i32_be(value as i32, buffer);
    } else {
        // 大长度格式
        inject_i32_be(-1, buffer);
        let high = (value >> 32) as i32;
        let low = (value & 0xFFFFFFFF) as i32;
        inject_i32_be(high, buffer);
        inject_i32_be(low, buffer);
    }
}

/// 写入长度到 writer
///
/// # 参数
/// * `value` - 要写入的长度值
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_length<W: Write>(value: usize, writer: &mut W) -> Result<()> {
    if value < 0x80000000 {
        // 标准格式
        write_i32(value as i32, writer)?;
    } else {
        // 大长度格式
        write_i32(-1, writer)?;
        let high = (value >> 32) as i32;
        let low = (value & 0xFFFFFFFF) as i32;
        write_i32(high, writer)?;
        write_i32(low, writer)?;
    }
    Ok(())
}

/// 注入字符串到 buffer（带长度前缀）
///
/// # 参数
/// * `s` - 要写入的字符串
/// * `buffer` - 目标缓冲区
pub fn inject_string(s: &str, buffer: &mut Vec<u8>) {
    let bytes = s.as_bytes();
    inject_length(bytes.len(), buffer);
    buffer.extend_from_slice(bytes);
}

/// 写入字符串到 writer（带长度前缀）
///
/// # 参数
/// * `s` - 要写入的字符串
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_string<W: Write>(s: &str, writer: &mut W) -> Result<()> {
    let bytes = s.as_bytes();
    write_length(bytes.len(), writer)?;
    writer.write_all(bytes)?;
    Ok(())
}

/// 构建对象头部
///
/// RDS 头部格式：
/// - 字节 0-1: 保留/gp 字段
/// - 字节 2: 标志位 (bit 1: has_attributes, bit 2: has_tag, bit 0: is_object)
/// - 字节 3: SEXP 类型
///
/// # 参数
/// * `sexp_type` - SEXP 类型
///
/// # 返回
/// 4 字节头部数组
#[inline]
pub fn make_header(sexp_type: SEXPType) -> Header {
    [0, 0, 0, sexp_type as u8]
}

/// 构建带属性标志的对象头部
///
/// # 参数
/// * `sexp_type` - SEXP 类型
/// * `has_attributes` - 是否有属性
/// * `is_object` - 是否有 class 属性（R 的 is.object() 标志）
/// * `is_s4` - 是否是 S4 对象（扩展基本类型的 S4 类，如 LogMap）
///
/// # 返回
/// 4 字节头部数组
#[inline]
pub fn make_header_with_flags(
    sexp_type: SEXPType,
    has_attributes: bool,
    is_object: bool,
    is_s4: bool,
) -> Header {
    let mut flags: u8 = 0x00;
    if is_object {
        flags |= 0x01;
    } // bit 0 = is_object
    if has_attributes {
        flags |= 0x02;
    } // bit 1 = has_attributes
      // gp high byte: bit 0 = S4 flag (corresponds to gp bit 4 in R)
    let gp_high: u8 = if is_s4 { 0x01 } else { 0x00 };
    [0, gp_high, flags, sexp_type as u8]
}

/// 注入对象头部到 buffer
///
/// # 参数
/// * `sexp_type` - SEXP 类型
/// * `buffer` - 目标缓冲区
pub fn inject_header(sexp_type: SEXPType, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&make_header(sexp_type));
}

/// 注入带属性的对象头部到 buffer
///
/// # 参数
/// * `sexp_type` - SEXP 类型
/// * `attributes` - 属性（用于判断是否有属性和 class）
/// * `buffer` - 目标缓冲区
pub fn inject_header_with_attributes(
    sexp_type: SEXPType,
    attributes: &Attributes,
    buffer: &mut Vec<u8>,
) {
    let has_attrs = !attributes.is_empty();
    let is_object = attributes.get("class").is_some();
    let is_s4 = is_s4_class(attributes);
    buffer.extend_from_slice(&make_header_with_flags(
        sexp_type, has_attrs, is_object, is_s4,
    ));
}

/// 写入对象头部到 writer
///
/// # 参数
/// * `sexp_type` - SEXP 类型
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_header<W: Write>(sexp_type: SEXPType, writer: &mut W) -> Result<()> {
    writer.write_all(&make_header(sexp_type))?;
    Ok(())
}

/// 写入带属性的对象头部到 writer
///
/// # 参数
/// * `sexp_type` - SEXP 类型
/// * `attributes` - 属性（用于判断是否有属性和 class）
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_header_with_attributes<W: Write>(
    sexp_type: SEXPType,
    attributes: &Attributes,
    writer: &mut W,
) -> Result<()> {
    let has_attrs = !attributes.is_empty();
    let is_object = attributes.get("class").is_some();
    let is_s4 = is_s4_class(attributes);
    writer.write_all(&make_header_with_flags(
        sexp_type, has_attrs, is_object, is_s4,
    ))?;
    Ok(())
}

/// Check if the class attribute indicates an S4 class.
/// S4 classes have a "package" sub-attribute on the class string vector.
fn is_s4_class(attributes: &Attributes) -> bool {
    if let Some(class_obj) = attributes.get("class") {
        if let RObject::StringVector(sv) = class_obj {
            // If the class StringVector has a "package" attribute, it's S4
            return sv.attributes.get("package").is_some();
        }
    }
    false
}

/// 构建配对列表节点头部
///
/// 配对列表头部格式：
/// - 字节 2 bit 2: has_tag 标志
/// - 字节 3: List 类型 (2)
///
/// # 参数
/// * `tagged` - 是否有标签
///
/// # 返回
/// 4 字节头部数组
#[inline]
pub fn make_pairlist_header(tagged: bool) -> Header {
    let flags = if tagged { 0x04 } else { 0x00 }; // bit 2 = has_tag
    [0, 0, flags, SEXPType::List as u8]
}

/// 注入配对列表节点头部到 buffer
///
/// # 参数
/// * `tagged` - 是否有标签
/// * `buffer` - 目标缓冲区
pub fn inject_pairlist_header(tagged: bool, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&make_pairlist_header(tagged));
}

/// 写入配对列表节点头部到 writer
///
/// # 参数
/// * `tagged` - 是否有标签
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_pairlist_header<W: Write>(tagged: bool, writer: &mut W) -> Result<()> {
    writer.write_all(&make_pairlist_header(tagged))?;
    Ok(())
}

/// 构建单字符串头部
///
/// 单字符串头部格式：
/// - 字节 1: ASCII 编码标志 (0x04)
/// - 字节 2: 其他编码标志 (0x20 = NONE, 0x40 = Latin1, 0x80 = UTF-8)
/// - 字节 3: Char 类型 (9)
///
/// 参考 rds2cpp write_single_string.hpp:
/// - NONE: buffer[2] = 1 << (12 - 8 + 1) = 0x20
/// - LATIN1: buffer[2] = 1 << (12 - 8 + 2) = 0x40
/// - UTF8: buffer[2] = 1 << (12 - 8 + 3) = 0x80
/// - ASCII: buffer[1] = 1 << (12 - 16 + 6) = 0x04
///
/// # 参数
/// * `encoding` - 字符串编码
/// * `missing` - 是否为 NA 值
///
/// # 返回
/// 4 字节头部数组
pub fn make_string_header(encoding: StringEncoding, missing: bool) -> Header {
    if missing {
        // NA 字符串使用特殊标志（长度为 -1）
        [0, 0, 0, SEXPType::Char as u8]
    } else {
        match encoding {
            StringEncoding::Utf8 => [0, 0, 0x80, SEXPType::Char as u8],
            StringEncoding::Latin1 => [0, 0, 0x40, SEXPType::Char as u8],
            StringEncoding::Ascii => [0, 0x04, 0, SEXPType::Char as u8],
            StringEncoding::None => [0, 0, 0x20, SEXPType::Char as u8],
        }
    }
}

/// 注入单字符串头部到 buffer
///
/// # 参数
/// * `encoding` - 字符串编码
/// * `missing` - 是否为 NA 值
/// * `buffer` - 目标缓冲区
pub fn inject_string_header(encoding: StringEncoding, missing: bool, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&make_string_header(encoding, missing));
}

/// 写入单字符串头部到 writer
///
/// # 参数
/// * `encoding` - 字符串编码
/// * `missing` - 是否为 NA 值
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_string_header<W: Write>(
    encoding: StringEncoding,
    missing: bool,
    writer: &mut W,
) -> Result<()> {
    writer.write_all(&make_string_header(encoding, missing))?;
    Ok(())
}

/// 构建引用头部
///
/// 引用头部格式：
/// - 字节 0-2: 引用索引（大端序）
/// - 字节 3: Ref 类型 (255)
///
/// # 参数
/// * `index` - 引用索引
///
/// # 返回
/// 4 字节头部数组
pub fn make_ref_header(index: usize) -> Header {
    let idx = index as u32;
    [
        ((idx >> 16) & 0xFF) as u8,
        ((idx >> 8) & 0xFF) as u8,
        (idx & 0xFF) as u8,
        SEXPType::Ref as u8,
    ]
}

/// 注入引用头部到 buffer
///
/// # 参数
/// * `index` - 引用索引
/// * `buffer` - 目标缓冲区
pub fn inject_ref_header(index: usize, buffer: &mut Vec<u8>) {
    buffer.extend_from_slice(&make_ref_header(index));
}

/// 写入引用头部到 writer
///
/// # 参数
/// * `index` - 引用索引
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_ref_header<W: Write>(index: usize, writer: &mut W) -> Result<()> {
    writer.write_all(&make_ref_header(index))?;
    Ok(())
}

/// 写入字节数组到 writer
///
/// # 参数
/// * `bytes` - 要写入的字节数组
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_bytes<W: Write>(bytes: &[u8], writer: &mut W) -> Result<()> {
    writer.write_all(bytes)?;
    Ok(())
}

/// 刷新 buffer 到 writer
///
/// # 参数
/// * `buffer` - 源缓冲区
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn flush_buffer<W: Write>(buffer: &mut Vec<u8>, writer: &mut W) -> Result<()> {
    writer.write_all(buffer)?;
    buffer.clear();
    Ok(())
}
