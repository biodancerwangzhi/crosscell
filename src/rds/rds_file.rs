//! RDS 文件结构
//!
//! 定义 RDS 文件的完整结构。
//! 对应 rds2cpp 的 RdsFile.hpp

use super::environment::Environment;
use super::external_pointer::ExternalPointer;
use super::r_object::RObject;
use super::symbol::Symbol;

/// RDS 文件内容
///
/// 表示一个完整的 RDS 文件，包含：
/// - 格式版本信息
/// - 根对象
/// - 符号、环境、外部指针的引用表
#[derive(Debug, Clone)]
pub struct RdsFile {
    /// 格式版本号（通常为 2 或 3）
    pub format_version: u32,
    /// 写入器版本（R 版本，如 [4, 2, 0]）
    pub writer_version: [u8; 3],
    /// 读取器最低版本要求
    pub reader_version: [u8; 3],
    /// 编码字符串（如 "UTF-8"）
    pub encoding: String,
    /// 根对象
    pub object: RObject,
    /// 环境向量（用于引用解析）
    pub environments: Vec<Environment>,
    /// 符号向量（用于引用解析）
    pub symbols: Vec<Symbol>,
    /// 外部指针向量（用于引用解析）
    pub external_pointers: Vec<ExternalPointer>,
}

impl Default for RdsFile {
    fn default() -> Self {
        Self {
            format_version: 3,
            writer_version: [4, 2, 0],
            reader_version: [3, 5, 0],
            encoding: "UTF-8".to_string(),
            object: RObject::Null,
            environments: Vec::new(),
            symbols: Vec::new(),
            external_pointers: Vec::new(),
        }
    }
}

impl RdsFile {
    /// 创建新的 RDS 文件结构
    pub fn new() -> Self {
        Self::default()
    }

    /// 使用指定的根对象创建 RDS 文件
    pub fn with_object(object: RObject) -> Self {
        Self {
            object,
            ..Self::default()
        }
    }

    /// 获取写入器版本字符串
    pub fn writer_version_string(&self) -> String {
        format!(
            "{}.{}.{}",
            self.writer_version[0], self.writer_version[1], self.writer_version[2]
        )
    }

    /// 获取读取器版本字符串
    pub fn reader_version_string(&self) -> String {
        format!(
            "{}.{}.{}",
            self.reader_version[0], self.reader_version[1], self.reader_version[2]
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rds_file_default() {
        let file = RdsFile::default();
        assert_eq!(file.format_version, 3);
        assert_eq!(file.writer_version, [4, 2, 0]);
        assert_eq!(file.reader_version, [3, 5, 0]);
        assert_eq!(file.encoding, "UTF-8");
        assert_eq!(file.object, RObject::Null);
        assert!(file.environments.is_empty());
        assert!(file.symbols.is_empty());
        assert!(file.external_pointers.is_empty());
    }

    #[test]
    fn test_rds_file_with_object() {
        use super::super::r_object::IntegerVector;
        let vec = IntegerVector { data: vec![1, 2, 3], ..Default::default() };
        let file = RdsFile::with_object(RObject::IntegerVector(vec));
        assert!(matches!(file.object, RObject::IntegerVector(_)));
    }

    #[test]
    fn test_version_strings() {
        let file = RdsFile::default();
        assert_eq!(file.writer_version_string(), "4.2.0");
        assert_eq!(file.reader_version_string(), "3.5.0");
    }
}
