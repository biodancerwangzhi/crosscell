//! RDS 错误类型定义
//!
//! 提供 RDS 解析和写入过程中的错误类型。

use thiserror::Error;

/// RDS 操作错误类型
///
/// 包含解析和写入 RDS 文件时可能发生的所有错误。
#[derive(Error, Debug)]
pub enum RdsError {
    /// IO 错误
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// 无效的 RDS 格式
    #[error("Invalid RDS format: {0}")]
    InvalidFormat(String),

    /// 不支持的 SEXP 类型
    #[error("Unsupported SEXP type: {0}")]
    UnsupportedType(u8),

    /// 意外的数据结束
    #[error("Unexpected end of data")]
    UnexpectedEof,

    /// 无效的引用索引
    #[error("Invalid reference index: {0}")]
    InvalidReference(usize),

    /// 解析错误（带上下文）
    #[error("Parse error in {context}: {message}")]
    ParseError {
        /// 错误发生的上下文
        context: String,
        /// 错误消息
        message: String,
    },

    /// 写入错误（带上下文）
    #[error("Write error in {context}: {message}")]
    WriteError {
        /// 错误发生的上下文
        context: String,
        /// 错误消息
        message: String,
    },

    /// 字符串过长
    #[error("String too long: {0} bytes (max 2^31-1)")]
    StringTooLong(usize),
}

/// RDS 操作结果类型别名
pub type Result<T> = std::result::Result<T, RdsError>;

/// 使用上下文信息包装解析错误
///
/// # 参数
/// * `context` - 错误发生的上下文描述
/// * `f` - 可能产生错误的闭包
///
/// # 返回
/// 如果闭包成功，返回其结果；否则返回带上下文的 ParseError
pub fn parse_with_context<T, F>(context: &str, f: F) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    f().map_err(|e| RdsError::ParseError {
        context: context.to_string(),
        message: e.to_string(),
    })
}

/// 使用上下文信息包装写入错误
///
/// # 参数
/// * `context` - 错误发生的上下文描述
/// * `f` - 可能产生错误的闭包
///
/// # 返回
/// 如果闭包成功，返回其结果；否则返回带上下文的 WriteError
pub fn write_with_context<T, F>(context: &str, f: F) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    f().map_err(|e| RdsError::WriteError {
        context: context.to_string(),
        message: e.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    #[test]
    fn test_io_error_conversion() {
        let io_err = io::Error::new(io::ErrorKind::NotFound, "file not found");
        let rds_err: RdsError = io_err.into();
        assert!(matches!(rds_err, RdsError::Io(_)));
        assert!(rds_err.to_string().contains("IO error"));
    }

    #[test]
    fn test_invalid_format_error() {
        let err = RdsError::InvalidFormat("not XDR format".to_string());
        assert_eq!(err.to_string(), "Invalid RDS format: not XDR format");
    }

    #[test]
    fn test_unsupported_type_error() {
        let err = RdsError::UnsupportedType(99);
        assert_eq!(err.to_string(), "Unsupported SEXP type: 99");
    }

    #[test]
    fn test_unexpected_eof_error() {
        let err = RdsError::UnexpectedEof;
        assert_eq!(err.to_string(), "Unexpected end of data");
    }

    #[test]
    fn test_invalid_reference_error() {
        let err = RdsError::InvalidReference(42);
        assert_eq!(err.to_string(), "Invalid reference index: 42");
    }

    #[test]
    fn test_parse_error() {
        let err = RdsError::ParseError {
            context: "header".to_string(),
            message: "invalid magic bytes".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "Parse error in header: invalid magic bytes"
        );
    }

    #[test]
    fn test_write_error() {
        let err = RdsError::WriteError {
            context: "integer vector".to_string(),
            message: "buffer overflow".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "Write error in integer vector: buffer overflow"
        );
    }

    #[test]
    fn test_string_too_long_error() {
        let err = RdsError::StringTooLong(3_000_000_000);
        assert_eq!(
            err.to_string(),
            "String too long: 3000000000 bytes (max 2^31-1)"
        );
    }

    #[test]
    fn test_parse_with_context() {
        let result: Result<i32> = parse_with_context("test context", || {
            Err(RdsError::UnexpectedEof)
        });
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, RdsError::ParseError { .. }));
        assert!(err.to_string().contains("test context"));
    }

    #[test]
    fn test_write_with_context() {
        let result: Result<i32> = write_with_context("test context", || {
            Err(RdsError::UnexpectedEof)
        });
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, RdsError::WriteError { .. }));
        assert!(err.to_string().contains("test context"));
    }

    #[test]
    fn test_result_type_alias() {
        fn returns_ok() -> Result<i32> {
            Ok(42)
        }

        fn returns_err() -> Result<i32> {
            Err(RdsError::UnexpectedEof)
        }

        assert_eq!(returns_ok().unwrap(), 42);
        assert!(returns_err().is_err());
    }
}
