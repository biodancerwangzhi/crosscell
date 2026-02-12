//! 外部指针结构
//!
//! 定义 R 外部指针（指向非 R 对象的指针）。
//! 对应 rds2cpp 的 ExternalPointer.hpp

use super::r_object::{Attributes, RObject};

/// 外部指针
///
/// 外部指针用于在 R 中引用非 R 对象（如 C/C++ 数据结构）。
/// 包含保护值、标签和属性。
#[derive(Debug, Clone)]
pub struct ExternalPointer {
    /// 保护值 - 防止被垃圾回收的 R 对象
    pub protection: Box<RObject>,
    /// 标签 - 用于标识指针类型的 R 对象
    pub tag: Box<RObject>,
    /// 属性
    pub attributes: Attributes,
}

impl Default for ExternalPointer {
    fn default() -> Self {
        Self {
            protection: Box::new(RObject::Null),
            tag: Box::new(RObject::Null),
            attributes: Attributes::default(),
        }
    }
}

impl PartialEq for ExternalPointer {
    fn eq(&self, other: &Self) -> bool {
        self.protection == other.protection
            && self.tag == other.tag
            && self.attributes == other.attributes
    }
}

impl ExternalPointer {
    /// 创建新的外部指针
    pub fn new() -> Self {
        Self::default()
    }

    /// 使用指定的保护值和标签创建外部指针
    pub fn with_values(protection: RObject, tag: RObject) -> Self {
        Self {
            protection: Box::new(protection),
            tag: Box::new(tag),
            attributes: Attributes::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_external_pointer_default() {
        let ptr = ExternalPointer::default();
        assert_eq!(*ptr.protection, RObject::Null);
        assert_eq!(*ptr.tag, RObject::Null);
        assert!(ptr.attributes.is_empty());
    }

    #[test]
    fn test_external_pointer_with_values() {
        let ptr = ExternalPointer::with_values(RObject::Null, RObject::Null);
        assert_eq!(*ptr.protection, RObject::Null);
        assert_eq!(*ptr.tag, RObject::Null);
    }

    #[test]
    fn test_external_pointer_equality() {
        let ptr1 = ExternalPointer::default();
        let ptr2 = ExternalPointer::default();
        assert_eq!(ptr1, ptr2);
    }
}
