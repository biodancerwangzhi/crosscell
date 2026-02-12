//! SEXP 类型枚举
//!
//! 定义 R 语言的所有数据类型（S-Expression Type）。
//! 对应 rds2cpp 的 SEXPType.hpp
//!
//! ## 类型分类
//!
//! ### 标准类型 (0-25)
//! 这些是 R 语言的核心数据类型，用于表示各种 R 对象。
//!
//! ### 序列化特殊类型 (238-255)
//! 这些类型仅在 RDS 序列化格式中使用，用于表示引用、特殊环境等。

/// R 数据类型枚举
///
/// 表示 R 语言中所有可能的 SEXP（S-Expression）类型。
/// 使用 `#[repr(u8)]` 确保与 RDS 二进制格式兼容。
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SEXPType {
    // ========================================
    // 标准类型 (0-25)
    // ========================================
    
    /// NULL 值 (NILSXP)
    Nil = 0,
    /// 符号 (SYMSXP) - 变量名、函数名等
    Sym = 1,
    /// 配对列表 (LISTSXP) - R 的链表结构
    List = 2,
    /// 闭包 (CLOSXP) - 用户定义的函数
    Clo = 3,
    /// 环境 (ENVSXP) - 变量绑定的作用域
    Env = 4,
    /// Promise (PROMSXP) - 延迟求值对象
    Prom = 5,
    /// 语言对象 (LANGSXP) - 函数调用表达式
    Lang = 6,
    /// 特殊函数 (SPECIALSXP) - 内置特殊函数
    Special = 7,
    /// 内置函数 (BUILTINSXP) - 内置普通函数
    Builtin = 8,
    /// 字符 (CHARSXP) - 单个字符串（内部使用）
    Char = 9,
    /// 逻辑向量 (LGLSXP)
    Lgl = 10,
    // 注意：11 和 12 在 R 中未使用
    /// 整数向量 (INTSXP)
    Int = 13,
    /// 实数向量 (REALSXP) - 双精度浮点数
    Real = 14,
    /// 复数向量 (CPLXSXP)
    Cplx = 15,
    /// 字符串向量 (STRSXP)
    Str = 16,
    /// 点参数 (DOTSXP) - ... 参数
    Dot = 17,
    /// ANY 类型 (ANYSXP) - 用于泛型函数
    Any = 18,
    /// 通用向量 (VECSXP) - R 的 list 类型
    Vec = 19,
    /// 表达式向量 (EXPRSXP)
    Expr = 20,
    /// 字节码 (BCODESXP) - 编译后的 R 代码
    Bcode = 21,
    /// 外部指针 (EXTPTRSXP)
    Extptr = 22,
    /// 弱引用 (WEAKREFSXP)
    Weakref = 23,
    /// 原始向量 (RAWSXP) - 字节向量
    Raw = 24,
    /// S4 对象 (S4SXP)
    S4 = 25,

    // ========================================
    // 序列化特殊类型 (238-255)
    // ========================================
    
    /// ALTREP 对象 - 替代表示（紧凑数据格式）
    Altrep = 238,
    /// 带属性的配对列表
    AttrList = 239,
    /// 带属性的语言对象
    AttrLang = 240,
    /// 基础环境 (base environment)
    BaseEnv = 241,
    /// 空环境 (empty environment)
    EmptyEnv = 242,
    /// 字节码引用定义
    BcRepRef = 243,
    /// 字节码重复定义
    BcRepDef = 244,
    /// 泛型引用
    GenericRef = 245,
    /// 类引用
    ClassRef = 246,
    /// 持久化对象
    Persist = 247,
    /// 包
    Package = 248,
    /// 命名空间
    Namespace = 249,
    /// 基础命名空间
    BaseNamespace = 250,
    /// 缺失参数
    MissingArg = 251,
    /// 未绑定值
    UnboundValue = 252,
    /// 全局环境 (global environment)
    GlobalEnv = 253,
    /// NIL 值（序列化格式）
    NilValue = 254,
    /// 引用 - 指向之前序列化的对象
    Ref = 255,
}

impl SEXPType {
    /// 从 u8 值创建 SEXPType
    ///
    /// # 参数
    /// * `value` - RDS 格式中的类型字节
    ///
    /// # 返回
    /// 如果值对应有效的 SEXP 类型，返回 `Some(SEXPType)`；
    /// 否则返回 `None`。
    #[inline]
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            // 标准类型 (0-25)
            0 => Some(Self::Nil),
            1 => Some(Self::Sym),
            2 => Some(Self::List),
            3 => Some(Self::Clo),
            4 => Some(Self::Env),
            5 => Some(Self::Prom),
            6 => Some(Self::Lang),
            7 => Some(Self::Special),
            8 => Some(Self::Builtin),
            9 => Some(Self::Char),
            10 => Some(Self::Lgl),
            // 11 和 12 未使用
            13 => Some(Self::Int),
            14 => Some(Self::Real),
            15 => Some(Self::Cplx),
            16 => Some(Self::Str),
            17 => Some(Self::Dot),
            18 => Some(Self::Any),
            19 => Some(Self::Vec),
            20 => Some(Self::Expr),
            21 => Some(Self::Bcode),
            22 => Some(Self::Extptr),
            23 => Some(Self::Weakref),
            24 => Some(Self::Raw),
            25 => Some(Self::S4),
            
            // 序列化特殊类型 (238-255)
            238 => Some(Self::Altrep),
            239 => Some(Self::AttrList),
            240 => Some(Self::AttrLang),
            241 => Some(Self::BaseEnv),
            242 => Some(Self::EmptyEnv),
            243 => Some(Self::BcRepRef),
            244 => Some(Self::BcRepDef),
            245 => Some(Self::GenericRef),
            246 => Some(Self::ClassRef),
            247 => Some(Self::Persist),
            248 => Some(Self::Package),
            249 => Some(Self::Namespace),
            250 => Some(Self::BaseNamespace),
            251 => Some(Self::MissingArg),
            252 => Some(Self::UnboundValue),
            253 => Some(Self::GlobalEnv),
            254 => Some(Self::NilValue),
            255 => Some(Self::Ref),
            
            // 未知类型
            _ => None,
        }
    }

    /// 检查是否为原子向量类型
    #[inline]
    pub fn is_atomic_vector(&self) -> bool {
        matches!(
            self,
            Self::Lgl | Self::Int | Self::Real | Self::Cplx | Self::Str | Self::Raw
        )
    }

    /// 检查是否为特殊环境类型
    #[inline]
    pub fn is_special_environment(&self) -> bool {
        matches!(
            self,
            Self::GlobalEnv | Self::BaseEnv | Self::EmptyEnv
        )
    }

    /// 检查是否为序列化特殊类型
    #[inline]
    pub fn is_serialization_type(&self) -> bool {
        (*self as u8) >= 238
    }
}

impl Default for SEXPType {
    fn default() -> Self {
        Self::Nil
    }
}

impl std::fmt::Display for SEXPType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_types_values() {
        // 验证标准类型的数值
        assert_eq!(SEXPType::Nil as u8, 0);
        assert_eq!(SEXPType::Sym as u8, 1);
        assert_eq!(SEXPType::List as u8, 2);
        assert_eq!(SEXPType::Clo as u8, 3);
        assert_eq!(SEXPType::Env as u8, 4);
        assert_eq!(SEXPType::Prom as u8, 5);
        assert_eq!(SEXPType::Lang as u8, 6);
        assert_eq!(SEXPType::Special as u8, 7);
        assert_eq!(SEXPType::Builtin as u8, 8);
        assert_eq!(SEXPType::Char as u8, 9);
        assert_eq!(SEXPType::Lgl as u8, 10);
        assert_eq!(SEXPType::Int as u8, 13);
        assert_eq!(SEXPType::Real as u8, 14);
        assert_eq!(SEXPType::Cplx as u8, 15);
        assert_eq!(SEXPType::Str as u8, 16);
        assert_eq!(SEXPType::Dot as u8, 17);
        assert_eq!(SEXPType::Any as u8, 18);
        assert_eq!(SEXPType::Vec as u8, 19);
        assert_eq!(SEXPType::Expr as u8, 20);
        assert_eq!(SEXPType::Bcode as u8, 21);
        assert_eq!(SEXPType::Extptr as u8, 22);
        assert_eq!(SEXPType::Weakref as u8, 23);
        assert_eq!(SEXPType::Raw as u8, 24);
        assert_eq!(SEXPType::S4 as u8, 25);
    }

    #[test]
    fn test_serialization_types_values() {
        // 验证序列化特殊类型的数值
        assert_eq!(SEXPType::Altrep as u8, 238);
        assert_eq!(SEXPType::AttrList as u8, 239);
        assert_eq!(SEXPType::AttrLang as u8, 240);
        assert_eq!(SEXPType::BaseEnv as u8, 241);
        assert_eq!(SEXPType::EmptyEnv as u8, 242);
        assert_eq!(SEXPType::BcRepRef as u8, 243);
        assert_eq!(SEXPType::BcRepDef as u8, 244);
        assert_eq!(SEXPType::GenericRef as u8, 245);
        assert_eq!(SEXPType::ClassRef as u8, 246);
        assert_eq!(SEXPType::Persist as u8, 247);
        assert_eq!(SEXPType::Package as u8, 248);
        assert_eq!(SEXPType::Namespace as u8, 249);
        assert_eq!(SEXPType::BaseNamespace as u8, 250);
        assert_eq!(SEXPType::MissingArg as u8, 251);
        assert_eq!(SEXPType::UnboundValue as u8, 252);
        assert_eq!(SEXPType::GlobalEnv as u8, 253);
        assert_eq!(SEXPType::NilValue as u8, 254);
        assert_eq!(SEXPType::Ref as u8, 255);
    }

    #[test]
    fn test_from_u8_standard_types() {
        // 测试标准类型的 from_u8
        assert_eq!(SEXPType::from_u8(0), Some(SEXPType::Nil));
        assert_eq!(SEXPType::from_u8(1), Some(SEXPType::Sym));
        assert_eq!(SEXPType::from_u8(10), Some(SEXPType::Lgl));
        assert_eq!(SEXPType::from_u8(13), Some(SEXPType::Int));
        assert_eq!(SEXPType::from_u8(14), Some(SEXPType::Real));
        assert_eq!(SEXPType::from_u8(19), Some(SEXPType::Vec));
        assert_eq!(SEXPType::from_u8(25), Some(SEXPType::S4));
    }

    #[test]
    fn test_from_u8_serialization_types() {
        // 测试序列化特殊类型的 from_u8
        assert_eq!(SEXPType::from_u8(238), Some(SEXPType::Altrep));
        assert_eq!(SEXPType::from_u8(253), Some(SEXPType::GlobalEnv));
        assert_eq!(SEXPType::from_u8(254), Some(SEXPType::NilValue));
        assert_eq!(SEXPType::from_u8(255), Some(SEXPType::Ref));
    }

    #[test]
    fn test_from_u8_invalid_values() {
        // 测试无效值返回 None
        assert_eq!(SEXPType::from_u8(11), None);  // 未使用
        assert_eq!(SEXPType::from_u8(12), None);  // 未使用
        assert_eq!(SEXPType::from_u8(26), None);  // 超出标准范围
        assert_eq!(SEXPType::from_u8(100), None); // 中间未定义区域
        assert_eq!(SEXPType::from_u8(237), None); // 序列化类型之前
    }

    #[test]
    fn test_is_atomic_vector() {
        assert!(SEXPType::Lgl.is_atomic_vector());
        assert!(SEXPType::Int.is_atomic_vector());
        assert!(SEXPType::Real.is_atomic_vector());
        assert!(SEXPType::Cplx.is_atomic_vector());
        assert!(SEXPType::Str.is_atomic_vector());
        assert!(SEXPType::Raw.is_atomic_vector());
        
        assert!(!SEXPType::Nil.is_atomic_vector());
        assert!(!SEXPType::Vec.is_atomic_vector());
        assert!(!SEXPType::List.is_atomic_vector());
    }

    #[test]
    fn test_is_special_environment() {
        assert!(SEXPType::GlobalEnv.is_special_environment());
        assert!(SEXPType::BaseEnv.is_special_environment());
        assert!(SEXPType::EmptyEnv.is_special_environment());
        
        assert!(!SEXPType::Env.is_special_environment());
        assert!(!SEXPType::Nil.is_special_environment());
    }

    #[test]
    fn test_is_serialization_type() {
        assert!(SEXPType::Altrep.is_serialization_type());
        assert!(SEXPType::Ref.is_serialization_type());
        assert!(SEXPType::NilValue.is_serialization_type());
        assert!(SEXPType::GlobalEnv.is_serialization_type());
        
        assert!(!SEXPType::Nil.is_serialization_type());
        assert!(!SEXPType::Int.is_serialization_type());
        assert!(!SEXPType::S4.is_serialization_type());
    }

    #[test]
    fn test_default() {
        assert_eq!(SEXPType::default(), SEXPType::Nil);
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", SEXPType::Nil), "Nil");
        assert_eq!(format!("{}", SEXPType::Int), "Int");
        assert_eq!(format!("{}", SEXPType::Ref), "Ref");
    }

    #[test]
    fn test_clone_and_copy() {
        let t1 = SEXPType::Int;
        let t2 = t1;  // Copy
        let t3 = t1.clone();  // Clone
        assert_eq!(t1, t2);
        assert_eq!(t1, t3);
    }

    #[test]
    fn test_eq_and_hash() {
        use std::collections::HashSet;
        
        let mut set = HashSet::new();
        set.insert(SEXPType::Int);
        set.insert(SEXPType::Real);
        set.insert(SEXPType::Int);  // 重复
        
        assert_eq!(set.len(), 2);
        assert!(set.contains(&SEXPType::Int));
        assert!(set.contains(&SEXPType::Real));
    }

    #[test]
    fn test_roundtrip_all_valid_types() {
        // 测试所有有效类型的往返转换
        let valid_values: Vec<u8> = (0..=25)
            .chain(238..=255)
            .filter(|&v| v != 11 && v != 12)  // 排除未使用的值
            .collect();
        
        for value in valid_values {
            let sexp_type = SEXPType::from_u8(value);
            assert!(sexp_type.is_some(), "Value {} should be valid", value);
            assert_eq!(sexp_type.unwrap() as u8, value, "Roundtrip failed for {}", value);
        }
    }
}
