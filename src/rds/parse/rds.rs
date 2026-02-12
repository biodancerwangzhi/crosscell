//! RDS 解析入口
//!
//! 提供 RDS 文件的主要解析入口函数。
//! 对应 rds2cpp 的 parse_rds 函数。

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use flate2::read::GzDecoder;
use crate::rds::error::{Result, RdsError};
use crate::rds::rds_file::RdsFile;
use super::utils::{read_i32_be, quick_extract};
use super::shared_info::SharedParseInfo;
use super::object::parse_object;

/// RDS 文件魔数
const RDS_MAGIC_XDR: &[u8] = b"X\n";
const RDS_MAGIC_ASCII: &[u8] = b"A\n";
const RDS_MAGIC_BINARY: &[u8] = b"B\n";

/// 解析选项
#[derive(Debug, Clone, Default)]
pub struct ParseRdsOptions {
    /// 是否启用并行解析（暂未实现）
    pub parallel: bool,
}

/// 从 Read 实现解析 RDS
///
/// # 参数
/// * `reader` - 实现 Read trait 的数据源
/// * `options` - 解析选项
///
/// # 返回
/// 解析后的 RdsFile 结构
///
/// # 错误
/// 如果格式无效或解析失败返回错误
pub fn parse_rds<R: Read>(mut reader: R, _options: &ParseRdsOptions) -> Result<RdsFile> {
    // 读取魔数（2 字节）
    let mut magic = [0u8; 2];
    reader.read_exact(&mut magic)?;
    
    // 验证魔数
    if magic != *RDS_MAGIC_XDR && magic != *RDS_MAGIC_ASCII && magic != *RDS_MAGIC_BINARY {
        return Err(RdsError::InvalidFormat(format!(
            "Invalid RDS magic number: {:?}", magic
        )));
    }
    
    // 目前只支持 XDR 格式
    if magic != *RDS_MAGIC_XDR {
        return Err(RdsError::InvalidFormat(format!(
            "Only XDR format is supported, got: {:?}", magic
        )));
    }
    
    // 读取格式版本（4 字节大端序）
    let mut version_bytes = [0u8; 4];
    reader.read_exact(&mut version_bytes)?;
    let format_version = read_i32_be(&version_bytes) as u32;
    
    // 验证版本
    if format_version < 2 || format_version > 3 {
        return Err(RdsError::InvalidFormat(format!(
            "Unsupported RDS format version: {}", format_version
        )));
    }
    
    // 读取写入器版本（4 字节）
    let mut writer_version_bytes = [0u8; 4];
    reader.read_exact(&mut writer_version_bytes)?;
    let writer_version = parse_r_version(read_i32_be(&writer_version_bytes) as u32);
    
    // 读取读取器版本（4 字节）
    let mut reader_version_bytes = [0u8; 4];
    reader.read_exact(&mut reader_version_bytes)?;
    let reader_version = parse_r_version(read_i32_be(&reader_version_bytes) as u32);
    
    // 版本 3 有编码字符串
    let encoding = if format_version >= 3 {
        // 读取编码字符串长度
        let mut len_bytes = [0u8; 4];
        reader.read_exact(&mut len_bytes)?;
        let len = read_i32_be(&len_bytes) as usize;
        
        // 读取编码字符串
        let bytes = quick_extract(&mut reader, len)?;
        String::from_utf8_lossy(&bytes).into_owned()
    } else {
        "native".to_string()
    };
    
    // 创建共享解析信息
    let mut shared = SharedParseInfo::new();
    
    // 解析根对象
    let object = parse_object(&mut reader, &mut shared)?;
    
    // 构建 RdsFile
    Ok(RdsFile {
        format_version,
        writer_version,
        reader_version,
        encoding,
        object,
        environments: shared.environments,
        symbols: shared.symbols,
        external_pointers: shared.external_pointers,
    })
}

/// 解析 R 版本号
///
/// R 版本号编码为单个整数：major * 65536 + minor * 256 + patch
fn parse_r_version(encoded: u32) -> [u8; 3] {
    let major = (encoded >> 16) as u8;
    let minor = ((encoded >> 8) & 0xFF) as u8;
    let patch = (encoded & 0xFF) as u8;
    [major, minor, patch]
}


/// 从文件路径解析 RDS
///
/// 自动检测并处理 gzip 压缩文件。
///
/// # 参数
/// * `path` - 文件路径
/// * `options` - 解析选项
///
/// # 返回
/// 解析后的 RdsFile 结构
///
/// # 错误
/// 如果文件不存在、格式无效或解析失败返回错误
pub fn parse_rds_file<P: AsRef<Path>>(path: P, options: &ParseRdsOptions) -> Result<RdsFile> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut buf_reader = BufReader::new(file);
    
    // 检查是否为 gzip 压缩
    let mut magic = [0u8; 2];
    buf_reader.read_exact(&mut magic)?;
    
    // gzip 魔数: 0x1f 0x8b
    let is_gzip = magic[0] == 0x1f && magic[1] == 0x8b;
    
    // 重新打开文件
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);
    
    if is_gzip {
        let decoder = GzDecoder::new(buf_reader);
        parse_rds(decoder, options)
    } else {
        parse_rds(buf_reader, options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::sexp_type::SEXPType;
    use crate::rds::r_object::RObject;

    fn make_rds_header_v3(encoding: &str) -> Vec<u8> {
        let mut data = Vec::new();
        // 魔数
        data.extend_from_slice(b"X\n");
        // 版本 3
        data.extend_from_slice(&3i32.to_be_bytes());
        // 写入器版本 4.2.0 = 4*65536 + 2*256 + 0 = 262656
        data.extend_from_slice(&262656i32.to_be_bytes());
        // 读取器版本 3.5.0 = 3*65536 + 5*256 + 0 = 197888
        data.extend_from_slice(&197888i32.to_be_bytes());
        // 编码字符串长度
        data.extend_from_slice(&(encoding.len() as i32).to_be_bytes());
        // 编码字符串
        data.extend_from_slice(encoding.as_bytes());
        data
    }

    fn make_rds_header_v2() -> Vec<u8> {
        let mut data = Vec::new();
        // 魔数
        data.extend_from_slice(b"X\n");
        // 版本 2
        data.extend_from_slice(&2i32.to_be_bytes());
        // 写入器版本
        data.extend_from_slice(&262656i32.to_be_bytes());
        // 读取器版本 3.5.0 = 3*65536 + 5*256 + 0 = 197888
        data.extend_from_slice(&197888i32.to_be_bytes());
        data
    }

    #[test]
    fn test_parse_r_version() {
        // 4.2.0 = 4*65536 + 2*256 + 0 = 262656
        assert_eq!(parse_r_version(262656), [4, 2, 0]);
        // 3.5.0 = 3*65536 + 5*256 + 0 = 197888
        assert_eq!(parse_r_version(197888), [3, 5, 0]);
        // 4.1.3 = 4*65536 + 1*256 + 3 = 262403
        assert_eq!(parse_r_version(262403), [4, 1, 3]);
    }

    #[test]
    fn test_parse_rds_null_v3() {
        let mut data = make_rds_header_v3("UTF-8");
        // NULL 对象
        data.extend_from_slice(&[0x00, 0x00, 0x00, SEXPType::NilValue as u8]);
        
        let cursor = Cursor::new(data);
        let options = ParseRdsOptions::default();
        let rds = parse_rds(cursor, &options).unwrap();
        
        assert_eq!(rds.format_version, 3);
        assert_eq!(rds.writer_version, [4, 2, 0]);
        assert_eq!(rds.reader_version, [3, 5, 0]);
        assert_eq!(rds.encoding, "UTF-8");
        assert!(matches!(rds.object, RObject::Null));
    }

    #[test]
    fn test_parse_rds_null_v2() {
        let mut data = make_rds_header_v2();
        // NULL 对象
        data.extend_from_slice(&[0x00, 0x00, 0x00, SEXPType::NilValue as u8]);
        
        let cursor = Cursor::new(data);
        let options = ParseRdsOptions::default();
        let rds = parse_rds(cursor, &options).unwrap();
        
        assert_eq!(rds.format_version, 2);
        assert_eq!(rds.encoding, "native");
        assert!(matches!(rds.object, RObject::Null));
    }

    #[test]
    fn test_parse_rds_integer_vector() {
        let mut data = make_rds_header_v3("UTF-8");
        // 整数向量头部
        data.extend_from_slice(&[0x00, 0x00, 0x00, SEXPType::Int as u8]);
        // 长度 3
        data.extend_from_slice(&3i32.to_be_bytes());
        // 数据
        data.extend_from_slice(&1i32.to_be_bytes());
        data.extend_from_slice(&2i32.to_be_bytes());
        data.extend_from_slice(&3i32.to_be_bytes());
        
        let cursor = Cursor::new(data);
        let options = ParseRdsOptions::default();
        let rds = parse_rds(cursor, &options).unwrap();
        
        if let RObject::IntegerVector(vec) = rds.object {
            assert_eq!(vec.data, vec![1, 2, 3]);
        } else {
            panic!("Expected IntegerVector");
        }
    }

    #[test]
    fn test_invalid_magic() {
        let data = vec![0x00, 0x00, 0x00, 0x00];
        let cursor = Cursor::new(data);
        let options = ParseRdsOptions::default();
        let result = parse_rds(cursor, &options);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_version() {
        let mut data = Vec::new();
        data.extend_from_slice(b"X\n");
        data.extend_from_slice(&99i32.to_be_bytes()); // 无效版本
        
        let cursor = Cursor::new(data);
        let options = ParseRdsOptions::default();
        let result = parse_rds(cursor, &options);
        assert!(result.is_err());
    }
}