//! RDS 写入入口
//!
//! 提供 RDS 文件的主要写入入口函数。
//! 对应 rds2cpp 的 write_rds 函数。

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;
use crate::rds::error::Result;
use crate::rds::rds_file::RdsFile;
use super::utils::{write_i32, write_length, write_bytes};
use super::shared_info::SharedWriteInfo;
use super::object::write_object;

/// RDS 文件魔数（XDR 格式）
const RDS_MAGIC_XDR: &[u8] = b"X\n";

/// 写入选项
#[derive(Debug, Clone)]
pub struct WriteRdsOptions {
    /// 是否启用 gzip 压缩
    pub compress: bool,
    /// 压缩级别 (0-9)
    pub compression_level: u32,
}

impl Default for WriteRdsOptions {
    fn default() -> Self {
        Self {
            compress: true,
            compression_level: 6,
        }
    }
}

/// 将 RdsFile 写入到 Write 实现
///
/// # 参数
/// * `rds` - 要写入的 RdsFile 结构
/// * `writer` - 实现 Write trait 的目标
///
/// # 错误
/// 如果写入失败返回 IO 错误
pub fn write_rds<W: Write>(rds: &RdsFile, writer: &mut W) -> Result<()> {
    // 写入魔数
    write_bytes(RDS_MAGIC_XDR, writer)?;

    // 写入格式版本
    write_i32(rds.format_version as i32, writer)?;

    // 写入写入器版本（编码为单个整数）
    let writer_version = encode_r_version(&rds.writer_version);
    write_i32(writer_version, writer)?;

    // 写入读取器版本
    let reader_version = encode_r_version(&rds.reader_version);
    write_i32(reader_version, writer)?;

    // 版本 3 需要写入编码字符串
    if rds.format_version >= 3 {
        write_length(rds.encoding.len(), writer)?;
        write_bytes(rds.encoding.as_bytes(), writer)?;
    }

    // 创建共享写入信息
    let mut shared = SharedWriteInfo::new(
        &rds.symbols,
        &rds.environments,
        &rds.external_pointers,
    );

    // 写入根对象
    write_object(&rds.object, writer, &mut shared)?;

    Ok(())
}

/// 编码 R 版本号
///
/// R 版本号编码为单个整数：major * 65536 + minor * 256 + patch
fn encode_r_version(version: &[u8; 3]) -> i32 {
    let major = version[0] as i32;
    let minor = version[1] as i32;
    let patch = version[2] as i32;
    (major << 16) | (minor << 8) | patch
}

/// 将 RdsFile 写入到文件
///
/// 根据选项决定是否使用 gzip 压缩。
///
/// # 参数
/// * `rds` - 要写入的 RdsFile 结构
/// * `path` - 输出文件路径
/// * `options` - 写入选项
///
/// # 错误
/// 如果文件创建失败或写入失败返回错误
pub fn write_rds_file<P: AsRef<Path>>(
    rds: &RdsFile,
    path: P,
    options: &WriteRdsOptions,
) -> Result<()> {
    let file = File::create(path)?;
    let buf_writer = BufWriter::new(file);

    if options.compress {
        let level = Compression::new(options.compression_level);
        let mut encoder = GzEncoder::new(buf_writer, level);
        write_rds(rds, &mut encoder)?;
        encoder.finish()?;
    } else {
        let mut writer = buf_writer;
        write_rds(rds, &mut writer)?;
    }

    Ok(())
}

/// 将 RdsFile 写入到文件（使用默认选项）
///
/// 默认启用 gzip 压缩，压缩级别为 6。
///
/// # 参数
/// * `rds` - 要写入的 RdsFile 结构
/// * `path` - 输出文件路径
///
/// # 错误
/// 如果文件创建失败或写入失败返回错误
pub fn write_rds_file_default<P: AsRef<Path>>(rds: &RdsFile, path: P) -> Result<()> {
    write_rds_file(rds, path, &WriteRdsOptions::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::rds::r_object::{RObject, IntegerVector, StringVector};
    use crate::rds::string_encoding::StringEncoding;
    use crate::rds::parse::rds::{parse_rds, ParseRdsOptions};

    fn make_rds_v3(object: RObject) -> RdsFile {
        RdsFile {
            format_version: 3,
            writer_version: [4, 2, 0],
            reader_version: [3, 5, 0],
            encoding: "UTF-8".to_string(),
            object,
            symbols: vec![],
            environments: vec![],
            external_pointers: vec![],
        }
    }

    fn make_rds_v2(object: RObject) -> RdsFile {
        RdsFile {
            format_version: 2,
            writer_version: [4, 2, 0],
            reader_version: [3, 5, 0],
            encoding: "native".to_string(),
            object,
            symbols: vec![],
            environments: vec![],
            external_pointers: vec![],
        }
    }

    #[test]
    fn test_encode_r_version() {
        // 4.2.0 = 4*65536 + 2*256 + 0 = 262656
        assert_eq!(encode_r_version(&[4, 2, 0]), 262656);
        // 3.5.0 = 3*65536 + 5*256 + 0 = 197888
        assert_eq!(encode_r_version(&[3, 5, 0]), 197888);
        // 4.1.3 = 4*65536 + 1*256 + 3 = 262403
        assert_eq!(encode_r_version(&[4, 1, 3]), 262403);
    }

    #[test]
    fn test_write_rds_null_v3() {
        let rds = make_rds_v3(RObject::Null);
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        // 魔数 (2) + 版本 (4) + 写入器版本 (4) + 读取器版本 (4) + 编码长度 (4) + 编码 (5) + NULL (4)
        assert!(!data.is_empty());
        assert_eq!(&data[0..2], b"X\n");
    }

    #[test]
    fn test_write_rds_null_v2() {
        let rds = make_rds_v2(RObject::Null);
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        assert!(!data.is_empty());
        assert_eq!(&data[0..2], b"X\n");
    }

    #[test]
    fn test_roundtrip_null() {
        let rds = make_rds_v3(RObject::Null);
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        let parsed = parse_rds(Cursor::new(data), &ParseRdsOptions::default()).unwrap();

        assert_eq!(parsed.format_version, 3);
        assert_eq!(parsed.writer_version, [4, 2, 0]);
        assert_eq!(parsed.reader_version, [3, 5, 0]);
        assert_eq!(parsed.encoding, "UTF-8");
        assert!(matches!(parsed.object, RObject::Null));
    }

    #[test]
    fn test_roundtrip_integer_vector() {
        let vec = IntegerVector {
            data: vec![1, 2, 3, -1, 0, i32::MAX, i32::MIN],
            attributes: Default::default(),
        };
        let rds = make_rds_v3(RObject::IntegerVector(vec));
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        let parsed = parse_rds(Cursor::new(data), &ParseRdsOptions::default()).unwrap();

        if let RObject::IntegerVector(v) = parsed.object {
            assert_eq!(v.data, vec![1, 2, 3, -1, 0, i32::MAX, i32::MIN]);
        } else {
            panic!("Expected IntegerVector");
        }
    }

    #[test]
    fn test_roundtrip_string_vector() {
        let vec = StringVector {
            data: vec!["hello".to_string(), "world".to_string(), "".to_string()],
            encodings: vec![StringEncoding::Utf8, StringEncoding::Utf8, StringEncoding::Utf8],
            missing: vec![false, false, true],
            attributes: Default::default(),
        };
        let rds = make_rds_v3(RObject::StringVector(vec));
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        let parsed = parse_rds(Cursor::new(data), &ParseRdsOptions::default()).unwrap();

        if let RObject::StringVector(v) = parsed.object {
            assert_eq!(v.data[0], "hello");
            assert_eq!(v.data[1], "world");
            assert!(v.missing[2]);
        } else {
            panic!("Expected StringVector");
        }
    }

    #[test]
    fn test_roundtrip_v2() {
        let vec = IntegerVector {
            data: vec![42],
            attributes: Default::default(),
        };
        let rds = make_rds_v2(RObject::IntegerVector(vec));
        let mut buffer = Cursor::new(Vec::new());

        write_rds(&rds, &mut buffer).unwrap();

        let data = buffer.into_inner();
        let parsed = parse_rds(Cursor::new(data), &ParseRdsOptions::default()).unwrap();

        assert_eq!(parsed.format_version, 2);
        if let RObject::IntegerVector(v) = parsed.object {
            assert_eq!(v.data, vec![42]);
        } else {
            panic!("Expected IntegerVector");
        }
    }
}
