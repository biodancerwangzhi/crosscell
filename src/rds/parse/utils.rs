//! 解析工具函数
//!
//! 提供底层解析工具，包括字节序转换、长度读取等。
//! 对应 rds2cpp 的解析工具函数。

use crate::rds::error::Result;
use std::io::Read;

/// 4 字节头部类型
pub type Header = [u8; 4];

/// 检查系统是否为小端序
///
/// RDS 使用大端序（XDR 格式），在小端系统上需要进行字节序转换。
///
/// # 返回
/// 如果系统为小端序返回 `true`，否则返回 `false`
#[inline]
pub fn is_little_endian() -> bool {
    cfg!(target_endian = "little")
}

/// 从大端序字节读取 i32
///
/// # 参数
/// * `bytes` - 4 字节大端序数据
///
/// # 返回
/// 转换后的 i32 值
#[inline]
pub fn read_i32_be(bytes: &[u8; 4]) -> i32 {
    i32::from_be_bytes(*bytes)
}

/// 从大端序字节读取 f64
///
/// # 参数
/// * `bytes` - 8 字节大端序数据
///
/// # 返回
/// 转换后的 f64 值
#[inline]
pub fn read_f64_be(bytes: &[u8; 8]) -> f64 {
    f64::from_be_bytes(*bytes)
}

/// 从 reader 读取大端序 i32
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 读取并转换后的 i32 值
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn read_i32<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(read_i32_be(&buf))
}

/// 从 reader 读取大端序 f64
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 读取并转换后的 f64 值
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn read_f64<R: Read>(reader: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(read_f64_be(&buf))
}

/// 读取向量长度
///
/// RDS 格式支持两种长度格式：
/// - 标准格式：4 字节大端序无符号整数（最大 2^31-1）
/// - 大长度格式：当标准长度为 -1 时，后续 8 字节为实际长度
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 读取的长度值
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn read_length<R: Read>(reader: &mut R) -> Result<usize> {
    let len = read_i32(reader)?;

    if len == -1 {
        // 大长度格式：读取 8 字节（两个 4 字节整数）
        let high = read_i32(reader)? as u32 as u64;
        let low = read_i32(reader)? as u32 as u64;
        Ok(((high << 32) | low) as usize)
    } else {
        Ok(len as usize)
    }
}

/// 快速提取指定长度的字节
///
/// 从 reader 读取指定长度的字节到新分配的 Vec 中。
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
/// * `len` - 要读取的字节数
///
/// # 返回
/// 包含读取字节的 Vec
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn quick_extract<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u8>> {
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

/// 将 u8 转换为 char
///
/// # 参数
/// * `val` - 要转换的字节值
///
/// # 返回
/// 对应的 ASCII 字符
#[inline]
pub fn as_char(val: u8) -> char {
    val as char
}

/// 读取 4 字节头部
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
///
/// # 返回
/// 4 字节头部数组
///
/// # 错误
/// 如果读取失败返回 IO 错误
pub fn read_header<R: Read>(reader: &mut R) -> Result<Header> {
    let mut header = [0u8; 4];
    reader.read_exact(&mut header)?;
    Ok(header)
}

/// 从头部提取 SEXP 类型
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// SEXP 类型字节
#[inline]
pub fn get_sexp_type(header: &Header) -> u8 {
    header[3]
}

/// 从头部提取标志字节
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// 标志字节
#[inline]
pub fn get_flags(header: &Header) -> u8 {
    header[2]
}

/// 从头部提取 gp 字段（用于编码信息）
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// gp 字段的高 8 位
#[inline]
pub fn get_gp_high(header: &Header) -> u8 {
    header[1]
}

/// 从头部提取 gp 字段的低 8 位
///
/// # 参数
/// * `header` - 4 字节头部
///
/// # 返回
/// gp 字段的低 8 位
#[inline]
pub fn get_gp_low(header: &Header) -> u8 {
    header[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_is_little_endian() {
        // 在大多数现代系统上应该返回 true
        let result = is_little_endian();
        // 只验证函数可以调用，不验证具体值（取决于平台）
        assert!(result || !result);
    }

    #[test]
    fn test_read_i32_be() {
        // 测试大端序转换
        assert_eq!(read_i32_be(&[0x00, 0x00, 0x00, 0x01]), 1);
        assert_eq!(read_i32_be(&[0x00, 0x00, 0x01, 0x00]), 256);
        assert_eq!(read_i32_be(&[0x7F, 0xFF, 0xFF, 0xFF]), i32::MAX);
        assert_eq!(read_i32_be(&[0x80, 0x00, 0x00, 0x00]), i32::MIN);
        assert_eq!(read_i32_be(&[0xFF, 0xFF, 0xFF, 0xFF]), -1);
    }

    #[test]
    fn test_read_f64_be() {
        // 测试 0.0
        assert_eq!(
            read_f64_be(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            0.0
        );
        // 测试 1.0 (IEEE 754: 0x3FF0000000000000)
        assert_eq!(
            read_f64_be(&[0x3F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            1.0
        );
        // 测试 -1.0 (IEEE 754: 0xBFF0000000000000)
        assert_eq!(
            read_f64_be(&[0xBF, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            -1.0
        );
    }

    #[test]
    fn test_read_i32_from_reader() {
        let data = vec![0x00, 0x00, 0x00, 0x2A]; // 42 in big-endian
        let mut cursor = Cursor::new(data);
        assert_eq!(read_i32(&mut cursor).unwrap(), 42);
    }

    #[test]
    fn test_read_f64_from_reader() {
        // 1.0 in big-endian IEEE 754
        let data = vec![0x3F, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut cursor = Cursor::new(data);
        assert_eq!(read_f64(&mut cursor).unwrap(), 1.0);
    }

    #[test]
    fn test_read_length_standard() {
        // 标准长度格式
        let data = vec![0x00, 0x00, 0x00, 0x64]; // 100
        let mut cursor = Cursor::new(data);
        assert_eq!(read_length(&mut cursor).unwrap(), 100);
    }

    #[test]
    fn test_read_length_large() {
        // 大长度格式：-1 后跟 8 字节
        let mut data = vec![0xFF, 0xFF, 0xFF, 0xFF]; // -1
        data.extend_from_slice(&[0x00, 0x00, 0x00, 0x01]); // high = 1
        data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]); // low = 0
        let mut cursor = Cursor::new(data);
        // 结果应该是 (1 << 32) = 4294967296
        assert_eq!(read_length(&mut cursor).unwrap(), 4294967296);
    }

    #[test]
    fn test_quick_extract() {
        let data = vec![0x01, 0x02, 0x03, 0x04, 0x05];
        let mut cursor = Cursor::new(data);
        let result = quick_extract(&mut cursor, 3).unwrap();
        assert_eq!(result, vec![0x01, 0x02, 0x03]);
    }

    #[test]
    fn test_quick_extract_empty() {
        let data = vec![0x01, 0x02];
        let mut cursor = Cursor::new(data);
        let result = quick_extract(&mut cursor, 0).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_as_char() {
        assert_eq!(as_char(65), 'A');
        assert_eq!(as_char(97), 'a');
        assert_eq!(as_char(48), '0');
    }

    #[test]
    fn test_read_header() {
        let data = vec![0x01, 0x02, 0x03, 0x04];
        let mut cursor = Cursor::new(data);
        let header = read_header(&mut cursor).unwrap();
        assert_eq!(header, [0x01, 0x02, 0x03, 0x04]);
    }

    #[test]
    fn test_get_sexp_type() {
        let header: Header = [0x00, 0x00, 0x00, 0x0D]; // Int type = 13
        assert_eq!(get_sexp_type(&header), 13);
    }

    #[test]
    fn test_get_flags() {
        let header: Header = [0x00, 0x00, 0x02, 0x00]; // has_attributes flag
        assert_eq!(get_flags(&header), 0x02);
    }

    #[test]
    fn test_get_gp_high() {
        let header: Header = [0x00, 0x08, 0x00, 0x00]; // UTF-8 encoding flag
        assert_eq!(get_gp_high(&header), 0x08);
    }

    #[test]
    fn test_get_gp_low() {
        let header: Header = [0x42, 0x00, 0x00, 0x00];
        assert_eq!(get_gp_low(&header), 0x42);
    }

    #[test]
    fn test_read_i32_eof() {
        let data = vec![0x00, 0x00]; // 只有 2 字节
        let mut cursor = Cursor::new(data);
        let result = read_i32(&mut cursor);
        assert!(result.is_err());
    }

    #[test]
    fn test_quick_extract_eof() {
        let data = vec![0x01, 0x02];
        let mut cursor = Cursor::new(data);
        let result = quick_extract(&mut cursor, 5); // 请求 5 字节但只有 2 字节
        assert!(result.is_err());
    }
}
