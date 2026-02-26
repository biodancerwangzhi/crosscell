//! R 符号结构
//!
//! 定义 R 符号（用于变量名和函数名）。
//! 对应 rds2cpp 的 Symbol.hpp

use super::string_encoding::StringEncoding;

/// R 符号
///
/// 符号用于表示 R 中的变量名、函数名等标识符。
/// 每个符号包含名称字符串和其编码信息。
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Symbol {
    /// 符号名称
    pub name: String,
    /// 名称的编码
    pub encoding: StringEncoding,
}

impl Symbol {
    /// 创建新符号
    pub fn new(name: String, encoding: StringEncoding) -> Self {
        Self { name, encoding }
    }

    /// 使用 UTF-8 编码创建符号
    pub fn utf8(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            encoding: StringEncoding::Utf8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_symbol_default() {
        let sym = Symbol::default();
        assert_eq!(sym.name, "");
        assert_eq!(sym.encoding, StringEncoding::Utf8);
    }

    #[test]
    fn test_symbol_new() {
        let sym = Symbol::new("test".to_string(), StringEncoding::Latin1);
        assert_eq!(sym.name, "test");
        assert_eq!(sym.encoding, StringEncoding::Latin1);
    }

    #[test]
    fn test_symbol_utf8() {
        let sym = Symbol::utf8("myvar");
        assert_eq!(sym.name, "myvar");
        assert_eq!(sym.encoding, StringEncoding::Utf8);
    }
}
