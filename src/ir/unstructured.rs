//! 非结构化数据 IR 类型
//!
//! 用于存储 AnnData.uns 和 Seurat@misc 中的任意数据

use std::collections::HashMap;
use std::fmt;

/// 非结构化数据值（用于 AnnData .uns 和 Seurat @misc）
///
/// 支持嵌套的字典和数组结构
#[derive(Debug, Clone, PartialEq)]
pub enum UnstructuredValue {
    String(String),
    Integer(i64),
    Float(f64),
    Boolean(bool),
    Array(Vec<UnstructuredValue>),
    Dict(HashMap<String, UnstructuredValue>),
    Null,
}

impl UnstructuredValue {
    /// 判断是否为 null
    pub fn is_null(&self) -> bool {
        matches!(self, UnstructuredValue::Null)
    }

    /// 尝试转换为字符串
    pub fn as_str(&self) -> Option<&str> {
        match self {
            UnstructuredValue::String(s) => Some(s),
            _ => None,
        }
    }

    /// 尝试转换为整数
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            UnstructuredValue::Integer(i) => Some(*i),
            _ => None,
        }
    }

    /// 尝试转换为浮点数
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            UnstructuredValue::Float(f) => Some(*f),
            UnstructuredValue::Integer(i) => Some(*i as f64),
            _ => None,
        }
    }

    /// 尝试转换为布尔值
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            UnstructuredValue::Boolean(b) => Some(*b),
            _ => None,
        }
    }

    /// 尝试转换为数组
    pub fn as_array(&self) -> Option<&Vec<UnstructuredValue>> {
        match self {
            UnstructuredValue::Array(arr) => Some(arr),
            _ => None,
        }
    }

    /// 尝试转换为字典
    pub fn as_dict(&self) -> Option<&HashMap<String, UnstructuredValue>> {
        match self {
            UnstructuredValue::Dict(dict) => Some(dict),
            _ => None,
        }
    }
}

impl fmt::Display for UnstructuredValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UnstructuredValue::String(s) => write!(f, "\"{}\"", s),
            UnstructuredValue::Integer(i) => write!(f, "{}", i),
            UnstructuredValue::Float(fl) => write!(f, "{}", fl),
            UnstructuredValue::Boolean(b) => write!(f, "{}", b),
            UnstructuredValue::Array(arr) => write!(f, "[{} items]", arr.len()),
            UnstructuredValue::Dict(dict) => write!(f, "{{{} keys}}", dict.len()),
            UnstructuredValue::Null => write!(f, "null"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unstructured_primitives() {
        let s = UnstructuredValue::String("hello".to_string());
        assert_eq!(s.as_str(), Some("hello"));

        let i = UnstructuredValue::Integer(42);
        assert_eq!(i.as_i64(), Some(42));
        assert_eq!(i.as_f64(), Some(42.0));

        let f = UnstructuredValue::Float(3.14);
        assert_eq!(f.as_f64(), Some(3.14));

        let b = UnstructuredValue::Boolean(true);
        assert_eq!(b.as_bool(), Some(true));

        let n = UnstructuredValue::Null;
        assert!(n.is_null());
    }

    #[test]
    fn test_unstructured_array() {
        let arr = UnstructuredValue::Array(vec![
            UnstructuredValue::Integer(1),
            UnstructuredValue::Integer(2),
            UnstructuredValue::Integer(3),
        ]);

        assert!(arr.as_array().is_some());
        assert_eq!(arr.as_array().unwrap().len(), 3);
    }

    #[test]
    fn test_unstructured_dict() {
        let mut dict = HashMap::new();
        dict.insert(
            "key1".to_string(),
            UnstructuredValue::String("value1".to_string()),
        );
        dict.insert("key2".to_string(), UnstructuredValue::Integer(42));

        let val = UnstructuredValue::Dict(dict);
        assert!(val.as_dict().is_some());
        assert_eq!(val.as_dict().unwrap().len(), 2);
    }

    #[test]
    fn test_nested_structure() {
        let mut inner_dict = HashMap::new();
        inner_dict.insert("inner_key".to_string(), UnstructuredValue::Float(3.14));

        let mut outer_dict = HashMap::new();
        outer_dict.insert("nested".to_string(), UnstructuredValue::Dict(inner_dict));
        outer_dict.insert(
            "array".to_string(),
            UnstructuredValue::Array(vec![
                UnstructuredValue::Integer(1),
                UnstructuredValue::Integer(2),
            ]),
        );

        let val = UnstructuredValue::Dict(outer_dict);
        let dict = val.as_dict().unwrap();

        assert!(dict.contains_key("nested"));
        assert!(dict.contains_key("array"));
        assert!(dict.get("nested").unwrap().as_dict().is_some());
        assert!(dict.get("array").unwrap().as_array().is_some());
    }
}
