//! 共享解析信息
//!
//! 维护解析过程中的共享状态，处理引用和循环结构。
//! 对应 rds2cpp 的 SharedParseInfo 结构。

use super::header::get_reference_index;
use super::utils::Header;
use crate::rds::environment::Environment;
use crate::rds::error::{RdsError, Result};
use crate::rds::external_pointer::ExternalPointer;
use crate::rds::r_object::{EnvironmentIndex, ExternalPointerIndex, RObject, SymbolIndex};
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use crate::rds::symbol::Symbol;

/// 存储的字符串信息
#[derive(Debug, Clone, Default)]
pub struct StoredString {
    pub value: String,
    pub encoding: StringEncoding,
    pub missing: bool,
}

/// 引用类型
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceType {
    /// 符号引用
    Symbol,
    /// 环境引用
    Environment,
    /// 外部指针引用
    ExternalPointer,
    /// 字符串引用 (CHARSXP)
    Char,
}

/// 共享解析信息
///
/// 在解析过程中维护符号、环境和外部指针的共享状态，
/// 支持引用解析和循环结构处理。
#[derive(Debug, Default)]
pub struct SharedParseInfo {
    /// 已解析的环境列表
    pub environments: Vec<Environment>,
    /// 已解析的符号列表
    pub symbols: Vec<Symbol>,
    /// 已解析的外部指针列表
    pub external_pointers: Vec<ExternalPointer>,
    /// 已解析的字符串列表 (CHARSXP)
    pub strings: Vec<StoredString>,
    /// 引用映射表：(类型, 索引)
    mappings: Vec<(ReferenceType, usize)>,
}

impl SharedParseInfo {
    /// 创建新的共享解析信息
    pub fn new() -> Self {
        Self::default()
    }

    /// 请求新符号索引
    ///
    /// 预分配一个符号槽位并返回其索引。
    /// 实际的符号数据将在后续解析中填充。
    ///
    /// # 返回
    /// 新分配的符号索引
    pub fn request_symbol(&mut self) -> usize {
        let index = self.symbols.len();
        self.symbols.push(Symbol::default());
        self.mappings.push((ReferenceType::Symbol, index));
        index
    }

    /// 获取符号索引（从引用头部）
    ///
    /// 从 REF 类型的头部提取引用索引，并验证它指向一个符号。
    /// 如果引用指向 CHARSXP，则创建一个新的 Symbol。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// 符号索引
    ///
    /// # 错误
    /// 如果引用索引无效返回错误
    pub fn get_symbol_index(&mut self, header: &Header) -> Result<usize> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            return Err(RdsError::InvalidReference(ref_index));
        }

        let (ref_type, index) = self.mappings[ref_index - 1];

        match ref_type {
            ReferenceType::Symbol => Ok(index),
            ReferenceType::Char => {
                // 引用指向 CHARSXP，创建一个新的 Symbol
                let stored = self
                    .strings
                    .get(index)
                    .ok_or_else(|| RdsError::InvalidReference(ref_index))?;
                let sym_index = self.symbols.len();
                self.symbols
                    .push(Symbol::new(stored.value.clone(), stored.encoding));
                // 注意：不添加到 mappings，因为这是一个派生的 Symbol
                Ok(sym_index)
            }
            _ => Err(RdsError::ParseError {
                context: "symbol reference".to_string(),
                message: format!("Expected symbol or char reference, got {:?}", ref_type),
            }),
        }
    }

    /// 请求新环境索引
    ///
    /// 预分配一个环境槽位并返回其索引。
    /// 实际的环境数据将在后续解析中填充。
    ///
    /// # 返回
    /// 新分配的环境索引
    pub fn request_environment(&mut self) -> usize {
        let index = self.environments.len();
        self.environments.push(Environment::default());
        self.mappings.push((ReferenceType::Environment, index));
        index
    }

    /// 获取环境索引（从引用头部）
    ///
    /// 从 REF 类型的头部提取引用索引，并验证它指向一个环境。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// 环境索引
    ///
    /// # 错误
    /// 如果引用索引无效或不指向环境返回错误
    pub fn get_environment_index(&self, header: &Header) -> Result<usize> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            return Err(RdsError::InvalidReference(ref_index));
        }

        let (ref_type, index) = self.mappings[ref_index - 1];
        if ref_type != ReferenceType::Environment {
            return Err(RdsError::ParseError {
                context: "environment reference".to_string(),
                message: format!("Expected environment reference, got {:?}", ref_type),
            });
        }

        Ok(index)
    }

    /// 请求新外部指针索引
    ///
    /// 预分配一个外部指针槽位并返回其索引。
    /// 实际的外部指针数据将在后续解析中填充。
    ///
    /// # 返回
    /// 新分配的外部指针索引
    pub fn request_external_pointer(&mut self) -> usize {
        let index = self.external_pointers.len();
        self.external_pointers.push(ExternalPointer::default());
        self.mappings.push((ReferenceType::ExternalPointer, index));
        index
    }

    /// 获取外部指针索引（从引用头部）
    ///
    /// 从 REF 类型的头部提取引用索引，并验证它指向一个外部指针。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// 外部指针索引
    ///
    /// # 错误
    /// 如果引用索引无效或不指向外部指针返回错误
    pub fn get_external_pointer_index(&self, header: &Header) -> Result<usize> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            return Err(RdsError::InvalidReference(ref_index));
        }

        let (ref_type, index) = self.mappings[ref_index - 1];
        if ref_type != ReferenceType::ExternalPointer {
            return Err(RdsError::ParseError {
                context: "external pointer reference".to_string(),
                message: format!("Expected external pointer reference, got {:?}", ref_type),
            });
        }

        Ok(index)
    }

    /// 请求新字符串索引
    ///
    /// 预分配一个字符串槽位并返回其索引。
    /// 实际的字符串数据将在后续解析中填充。
    ///
    /// 注意：CHARSXP 不添加到 mappings 中，因为 R 的引用机制
    /// 只用于符号、环境和外部指针，不用于 CHARSXP。
    ///
    /// # 返回
    /// 新分配的字符串索引
    pub fn request_string(&mut self) -> usize {
        let index = self.strings.len();
        self.strings.push(StoredString::default());
        // 注意：不添加到 mappings 中！
        // rds2cpp 也不将 CHARSXP 添加到 mappings 中
        index
    }

    /// 存储字符串
    ///
    /// 用实际解析的字符串数据更新预分配的槽位。
    ///
    /// # 参数
    /// * `index` - 字符串索引
    /// * `value` - 字符串值
    /// * `encoding` - 字符串编码
    /// * `missing` - 是否为 NA
    pub fn store_string(
        &mut self,
        index: usize,
        value: String,
        encoding: StringEncoding,
        missing: bool,
    ) {
        if index < self.strings.len() {
            self.strings[index] = StoredString {
                value,
                encoding,
                missing,
            };
        }
    }

    /// 获取字符串（从引用头部）
    ///
    /// 从 REF 类型的头部提取引用索引，并返回对应的字符串。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// 存储的字符串信息
    ///
    /// # 错误
    /// 如果引用索引无效返回错误
    pub fn get_string(&self, header: &Header) -> Result<&StoredString> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            return Err(RdsError::InvalidReference(ref_index));
        }

        let (ref_type, index) = self.mappings[ref_index - 1];

        match ref_type {
            ReferenceType::Char => self
                .strings
                .get(index)
                .ok_or_else(|| RdsError::InvalidReference(ref_index)),
            ReferenceType::Symbol => {
                // 符号也可以作为字符串使用
                // 返回一个临时的 StoredString 不太好，所以返回错误让调用者处理
                Err(RdsError::ParseError {
                    context: "string reference".to_string(),
                    message: format!("Reference points to symbol, not char (index {})", index),
                })
            }
            _ => Err(RdsError::ParseError {
                context: "string reference".to_string(),
                message: format!("Expected char reference, got {:?}", ref_type),
            }),
        }
    }

    /// 获取字符串或符号名（从引用头部）
    ///
    /// 从 REF 类型的头部提取引用索引，返回字符串值。
    /// 如果引用指向符号，返回符号名。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// (字符串值, 编码, 是否为NA)
    pub fn get_string_or_symbol(&self, header: &Header) -> Result<(String, StringEncoding, bool)> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            return Err(RdsError::InvalidReference(ref_index));
        }

        let (ref_type, index) = self.mappings[ref_index - 1];

        match ref_type {
            ReferenceType::Char => {
                let stored = self
                    .strings
                    .get(index)
                    .ok_or_else(|| RdsError::InvalidReference(ref_index))?;
                Ok((stored.value.clone(), stored.encoding, stored.missing))
            }
            ReferenceType::Symbol => {
                let sym = self
                    .symbols
                    .get(index)
                    .ok_or_else(|| RdsError::InvalidReference(ref_index))?;
                Ok((sym.name.clone(), sym.encoding, false))
            }
            ReferenceType::Environment => {
                // 环境引用 - 返回占位符
                Ok(("<environment>".to_string(), StringEncoding::Utf8, false))
            }
            ReferenceType::ExternalPointer => {
                // 外部指针引用 - 返回占位符
                Ok(("<externalptr>".to_string(), StringEncoding::Utf8, false))
            }
        }
    }

    /// 解析引用
    ///
    /// 从 REF 类型的头部解析引用，返回对应的 RObject。
    /// 如果引用索引无效，返回 Null 而不是错误，以支持字节码解析中的容错。
    ///
    /// # 参数
    /// * `header` - REF 类型的 4 字节头部
    ///
    /// # 返回
    /// 引用对应的 RObject，如果引用无效则返回 Null
    pub fn resolve_reference(&self, header: &Header) -> Result<RObject> {
        let ref_index = get_reference_index(header);

        if ref_index == 0 || ref_index > self.mappings.len() {
            // 返回 Null 而不是错误，以支持字节码解析中的容错
            return Ok(RObject::Null);
        }

        let (ref_type, index) = self.mappings[ref_index - 1];

        match ref_type {
            ReferenceType::Symbol => Ok(RObject::SymbolIndex(SymbolIndex { index })),
            ReferenceType::Environment => {
                // 对于环境引用，我们需要确定环境类型
                // 默认使用 Env 类型，实际类型在解析时确定
                Ok(RObject::EnvironmentIndex(EnvironmentIndex {
                    index,
                    env_type: SEXPType::Env,
                }))
            }
            ReferenceType::ExternalPointer => {
                Ok(RObject::ExternalPointerIndex(ExternalPointerIndex {
                    index,
                }))
            }
            ReferenceType::Char => {
                // 字符串引用 - 返回字符串向量
                if let Some(stored) = self.strings.get(index) {
                    use crate::rds::r_object::StringVector;
                    Ok(RObject::StringVector(StringVector {
                        data: vec![stored.value.clone()],
                        encodings: vec![stored.encoding],
                        missing: vec![stored.missing],
                        attributes: Default::default(),
                    }))
                } else {
                    // 返回 Null 而不是错误
                    Ok(RObject::Null)
                }
            }
        }
    }

    /// 更新符号
    ///
    /// 用实际解析的符号数据更新预分配的槽位。
    ///
    /// # 参数
    /// * `index` - 符号索引
    /// * `symbol` - 符号数据
    pub fn update_symbol(&mut self, index: usize, symbol: Symbol) {
        if index < self.symbols.len() {
            self.symbols[index] = symbol;
        }
    }

    /// 更新环境
    ///
    /// 用实际解析的环境数据更新预分配的槽位。
    ///
    /// # 参数
    /// * `index` - 环境索引
    /// * `environment` - 环境数据
    pub fn update_environment(&mut self, index: usize, environment: Environment) {
        if index < self.environments.len() {
            self.environments[index] = environment;
        }
    }

    /// 更新外部指针
    ///
    /// 用实际解析的外部指针数据更新预分配的槽位。
    ///
    /// # 参数
    /// * `index` - 外部指针索引
    /// * `external_pointer` - 外部指针数据
    pub fn update_external_pointer(&mut self, index: usize, external_pointer: ExternalPointer) {
        if index < self.external_pointers.len() {
            self.external_pointers[index] = external_pointer;
        }
    }

    /// 获取当前引用计数
    ///
    /// 返回已注册的引用总数。
    pub fn reference_count(&self) -> usize {
        self.mappings.len()
    }

    /// 获取符号
    ///
    /// # 参数
    /// * `index` - 符号索引
    ///
    /// # 返回
    /// 符号引用，如果索引无效返回 None
    pub fn get_symbol(&self, index: usize) -> Option<&Symbol> {
        self.symbols.get(index)
    }

    /// 获取环境
    ///
    /// # 参数
    /// * `index` - 环境索引
    ///
    /// # 返回
    /// 环境引用，如果索引无效返回 None
    pub fn get_environment(&self, index: usize) -> Option<&Environment> {
        self.environments.get(index)
    }

    /// 获取外部指针
    ///
    /// # 参数
    /// * `index` - 外部指针索引
    ///
    /// # 返回
    /// 外部指针引用，如果索引无效返回 None
    pub fn get_external_pointer(&self, index: usize) -> Option<&ExternalPointer> {
        self.external_pointers.get(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::string_encoding::StringEncoding;

    #[test]
    fn test_shared_parse_info_new() {
        let info = SharedParseInfo::new();
        assert!(info.symbols.is_empty());
        assert!(info.environments.is_empty());
        assert!(info.external_pointers.is_empty());
        assert_eq!(info.reference_count(), 0);
    }

    #[test]
    fn test_request_symbol() {
        let mut info = SharedParseInfo::new();

        let idx1 = info.request_symbol();
        assert_eq!(idx1, 0);
        assert_eq!(info.symbols.len(), 1);
        assert_eq!(info.reference_count(), 1);

        let idx2 = info.request_symbol();
        assert_eq!(idx2, 1);
        assert_eq!(info.symbols.len(), 2);
        assert_eq!(info.reference_count(), 2);
    }

    #[test]
    fn test_request_environment() {
        let mut info = SharedParseInfo::new();

        let idx1 = info.request_environment();
        assert_eq!(idx1, 0);
        assert_eq!(info.environments.len(), 1);
        assert_eq!(info.reference_count(), 1);

        let idx2 = info.request_environment();
        assert_eq!(idx2, 1);
        assert_eq!(info.environments.len(), 2);
        assert_eq!(info.reference_count(), 2);
    }

    #[test]
    fn test_request_external_pointer() {
        let mut info = SharedParseInfo::new();

        let idx1 = info.request_external_pointer();
        assert_eq!(idx1, 0);
        assert_eq!(info.external_pointers.len(), 1);
        assert_eq!(info.reference_count(), 1);

        let idx2 = info.request_external_pointer();
        assert_eq!(idx2, 1);
        assert_eq!(info.external_pointers.len(), 2);
        assert_eq!(info.reference_count(), 2);
    }

    #[test]
    fn test_get_symbol_index() {
        let mut info = SharedParseInfo::new();
        info.request_symbol();

        // 引用索引从 1 开始，存储在头部的前 3 字节（大端序）
        let header: Header = [0x00, 0x00, 0x01, 0xFF]; // ref_index = 1
        let result = info.get_symbol_index(&header);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 0);
    }

    #[test]
    fn test_get_symbol_index_invalid() {
        let mut info = SharedParseInfo::new();

        // 没有注册任何引用
        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.get_symbol_index(&header);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_symbol_index_wrong_type() {
        let mut info = SharedParseInfo::new();
        info.request_environment(); // 注册环境而不是符号

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.get_symbol_index(&header);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_symbol_index_from_char() {
        let mut info = SharedParseInfo::new();
        // request_string does NOT add to mappings, so we need to use request_symbol
        let sym_idx = info.request_symbol();
        let symbol = Symbol::new("test_symbol".to_string(), StringEncoding::Utf8);
        info.update_symbol(sym_idx, symbol);

        // ref_index 1 -> first mapping (symbol at index 0)
        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.get_symbol_index(&header);
        assert!(result.is_ok());
        assert_eq!(
            info.get_symbol(result.unwrap()).unwrap().name,
            "test_symbol"
        );
    }

    #[test]
    fn test_get_environment_index() {
        let mut info = SharedParseInfo::new();
        info.request_environment();

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.get_environment_index(&header);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 0);
    }

    #[test]
    fn test_get_external_pointer_index() {
        let mut info = SharedParseInfo::new();
        info.request_external_pointer();

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.get_external_pointer_index(&header);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 0);
    }

    #[test]
    fn test_resolve_reference_symbol() {
        let mut info = SharedParseInfo::new();
        info.request_symbol();

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.resolve_reference(&header);
        assert!(result.is_ok());

        match result.unwrap() {
            RObject::SymbolIndex(idx) => assert_eq!(idx.index, 0),
            _ => panic!("Expected SymbolIndex"),
        }
    }

    #[test]
    fn test_resolve_reference_environment() {
        let mut info = SharedParseInfo::new();
        info.request_environment();

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.resolve_reference(&header);
        assert!(result.is_ok());

        match result.unwrap() {
            RObject::EnvironmentIndex(idx) => {
                assert_eq!(idx.index, 0);
                assert_eq!(idx.env_type, SEXPType::Env);
            }
            _ => panic!("Expected EnvironmentIndex"),
        }
    }

    #[test]
    fn test_resolve_reference_external_pointer() {
        let mut info = SharedParseInfo::new();
        info.request_external_pointer();

        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.resolve_reference(&header);
        assert!(result.is_ok());

        match result.unwrap() {
            RObject::ExternalPointerIndex(idx) => assert_eq!(idx.index, 0),
            _ => panic!("Expected ExternalPointerIndex"),
        }
    }

    #[test]
    fn test_resolve_reference_invalid() {
        let info = SharedParseInfo::new();

        // resolve_reference returns Ok(Null) for invalid refs (bytecode tolerance)
        let header: Header = [0x00, 0x00, 0x01, 0xFF];
        let result = info.resolve_reference(&header);
        assert!(result.is_ok());
        assert!(matches!(result.unwrap(), RObject::Null));
    }

    #[test]
    fn test_update_symbol() {
        let mut info = SharedParseInfo::new();
        let idx = info.request_symbol();

        let symbol = Symbol::new("test".to_string(), StringEncoding::Utf8);
        info.update_symbol(idx, symbol.clone());

        assert_eq!(info.get_symbol(idx).unwrap().name, "test");
    }

    #[test]
    fn test_update_environment() {
        let mut info = SharedParseInfo::new();
        let idx = info.request_environment();

        let mut env = Environment::new();
        env.locked = true;
        info.update_environment(idx, env);

        assert!(info.get_environment(idx).unwrap().locked);
    }

    #[test]
    fn test_update_external_pointer() {
        let mut info = SharedParseInfo::new();
        let idx = info.request_external_pointer();

        let ptr = ExternalPointer::with_values(RObject::Null, RObject::Null);
        info.update_external_pointer(idx, ptr);

        assert!(info.get_external_pointer(idx).is_some());
    }

    #[test]
    fn test_mixed_references() {
        let mut info = SharedParseInfo::new();

        // 注册不同类型的引用
        let sym_idx = info.request_symbol();
        let env_idx = info.request_environment();
        let ptr_idx = info.request_external_pointer();

        assert_eq!(sym_idx, 0);
        assert_eq!(env_idx, 0);
        assert_eq!(ptr_idx, 0);
        assert_eq!(info.reference_count(), 3);

        // 验证引用解析（大端序）
        let header1: Header = [0x00, 0x00, 0x01, 0xFF]; // ref 1 -> symbol
        let header2: Header = [0x00, 0x00, 0x02, 0xFF]; // ref 2 -> environment
        let header3: Header = [0x00, 0x00, 0x03, 0xFF]; // ref 3 -> external pointer

        assert!(matches!(
            info.resolve_reference(&header1).unwrap(),
            RObject::SymbolIndex(_)
        ));
        assert!(matches!(
            info.resolve_reference(&header2).unwrap(),
            RObject::EnvironmentIndex(_)
        ));
        assert!(matches!(
            info.resolve_reference(&header3).unwrap(),
            RObject::ExternalPointerIndex(_)
        ));
    }

    #[test]
    fn test_reference_index_zero() {
        let mut info = SharedParseInfo::new();
        info.request_symbol();

        // resolve_reference returns Ok(Null) for ref_index 0 (bytecode tolerance)
        let header: Header = [0x00, 0x00, 0x00, 0xFF];
        let result = info.resolve_reference(&header);
        assert!(result.is_ok());
        assert!(matches!(result.unwrap(), RObject::Null));
    }
}
