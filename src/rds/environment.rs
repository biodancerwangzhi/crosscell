//! R 环境结构
//!
//! 定义 R 环境（存储变量绑定）。
//! 对应 rds2cpp 的 Environment.hpp

use super::r_object::{Attributes, RObject};
use super::sexp_type::SEXPType;
use super::string_encoding::StringEncoding;

/// R 环境
///
/// 环境是 R 中存储变量绑定的数据结构。
/// 每个环境包含变量名-值对，以及指向父环境的引用。
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Environment {
    /// 环境是否被锁定（不可修改）
    pub locked: bool,
    /// 环境是否使用哈希表存储
    pub hashed: bool,
    /// 父环境类型
    pub parent_type: SEXPType,
    /// 父环境索引（对于特殊环境为 usize::MAX）
    pub parent: usize,
    /// 变量名称
    pub variable_names: Vec<String>,
    /// 变量名称编码
    pub variable_encodings: Vec<StringEncoding>,
    /// 变量值
    pub variable_values: Vec<RObject>,
    /// 环境属性
    pub attributes: Attributes,
}

impl Environment {
    /// 创建新的空环境
    pub fn new() -> Self {
        Self::default()
    }

    /// 添加变量
    pub fn add_variable(&mut self, name: String, value: RObject, encoding: StringEncoding) {
        self.variable_names.push(name);
        self.variable_values.push(value);
        self.variable_encodings.push(encoding);
    }

    /// 获取变量数量
    pub fn len(&self) -> usize {
        self.variable_names.len()
    }

    /// 检查是否为空
    pub fn is_empty(&self) -> bool {
        self.variable_names.is_empty()
    }

    /// 根据名称获取变量值
    pub fn get(&self, name: &str) -> Option<&RObject> {
        self.variable_names
            .iter()
            .position(|n| n == name)
            .map(|i| &self.variable_values[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_environment_default() {
        let env = Environment::default();
        assert!(!env.locked);
        assert!(!env.hashed);
        assert!(env.is_empty());
        assert_eq!(env.len(), 0);
    }

    #[test]
    fn test_environment_add_variable() {
        let mut env = Environment::new();
        env.add_variable("x".to_string(), RObject::Null, StringEncoding::Utf8);
        assert_eq!(env.len(), 1);
        assert!(!env.is_empty());
        assert!(env.get("x").is_some());
        assert!(env.get("y").is_none());
    }
}
