//! 共享写入信息
//!
//! 维护写入过程中的共享状态，处理引用和去重。
//! 对应 rds2cpp 的写入共享信息结构。

use std::collections::HashMap;
use std::io::Write;

use crate::rds::error::Result;
use crate::rds::sexp_type::SEXPType;
use crate::rds::string_encoding::StringEncoding;
use crate::rds::symbol::Symbol;
use crate::rds::environment::Environment;
use crate::rds::external_pointer::ExternalPointer;
use crate::rds::r_object::{SymbolIndex, EnvironmentIndex, ExternalPointerIndex};

use super::utils::{
    write_header, write_ref_header, write_string_header, write_length, write_bytes,
};

/// 共享写入信息
///
/// 在写入过程中维护符号、环境和外部指针的共享状态，
/// 支持引用去重和循环结构处理。
pub struct SharedWriteInfo<'a> {
    /// 当前引用计数（从 1 开始）
    reference_count: usize,
    /// 符号映射表：(名称, 编码) -> 引用索引
    symbol_mappings: HashMap<(String, StringEncoding), usize>,
    /// 已知符号列表（来自 RdsFile）
    known_symbols: &'a [Symbol],
    /// 已知环境列表（来自 RdsFile）
    known_environments: &'a [Environment],
    /// 已知外部指针列表（来自 RdsFile）
    known_external_pointers: &'a [ExternalPointer],
    /// 已知符号的引用索引映射
    known_symbol_mappings: Vec<Option<usize>>,
    /// 已知环境的引用索引映射
    known_environment_mappings: Vec<Option<usize>>,
    /// 已知外部指针的引用索引映射
    known_external_pointer_mappings: Vec<Option<usize>>,
}

impl<'a> SharedWriteInfo<'a> {
    /// 创建新的共享写入信息
    ///
    /// # 参数
    /// * `symbols` - 已知符号列表
    /// * `environments` - 已知环境列表
    /// * `external_pointers` - 已知外部指针列表
    pub fn new(
        symbols: &'a [Symbol],
        environments: &'a [Environment],
        external_pointers: &'a [ExternalPointer],
    ) -> Self {
        Self {
            // 从 1 开始，与 rds2cpp 保持一致
            // 参考 rds2cpp SharedWriteInfo.hpp: reference_count(1)
            reference_count: 1,
            symbol_mappings: HashMap::new(),
            known_symbols: symbols,
            known_environments: environments,
            known_external_pointers: external_pointers,
            known_symbol_mappings: vec![None; symbols.len()],
            known_environment_mappings: vec![None; environments.len()],
            known_external_pointer_mappings: vec![None; external_pointers.len()],
        }
    }


    /// 分配新的引用索引
    /// 
    /// 与 rds2cpp 保持一致：先返回当前值，再增加计数。
    /// 参考：refs/rds2cpp/include/rds2cpp/SharedWriteInfo.hpp
    fn allocate_reference(&mut self) -> usize {
        let index = self.reference_count;
        self.reference_count += 1;
        index
    }

    /// 获取当前引用计数
    pub fn reference_count(&self) -> usize {
        self.reference_count
    }

    /// 写入符号（带去重）
    ///
    /// 如果符号已存在，写入引用；否则写入完整符号并注册。
    ///
    /// # 参数
    /// * `name` - 符号名称
    /// * `encoding` - 符号编码
    /// * `writer` - 写入目标
    ///
    /// # 返回
    /// 符号的引用索引
    pub fn write_symbol<W: Write>(
        &mut self,
        name: &str,
        encoding: StringEncoding,
        writer: &mut W,
    ) -> Result<usize> {
        let key = (name.to_string(), encoding);
        
        // 检查是否已存在
        if let Some(&ref_index) = self.symbol_mappings.get(&key) {
            // 写入引用
            write_ref_header(ref_index, writer)?;
            return Ok(ref_index);
        }
        
        // 分配新引用索引
        let ref_index = self.allocate_reference();
        self.symbol_mappings.insert(key, ref_index);
        
        // 写入符号头部
        write_header(SEXPType::Sym, writer)?;
        
        // 写入符号名称（作为单字符串）
        write_string_header(encoding, false, writer)?;
        write_length(name.len(), writer)?;
        write_bytes(name.as_bytes(), writer)?;
        
        Ok(ref_index)
    }

    /// 写入符号引用
    ///
    /// 根据 SymbolIndex 写入对应的符号。
    ///
    /// # 参数
    /// * `obj` - 符号索引
    /// * `writer` - 写入目标
    pub fn write_symbol_ref<W: Write>(
        &mut self,
        obj: &SymbolIndex,
        writer: &mut W,
    ) -> Result<()> {
        let index = obj.index;
        
        // 检查是否已写入过
        if let Some(ref_index) = self.known_symbol_mappings.get(index).copied().flatten() {
            write_ref_header(ref_index, writer)?;
            return Ok(());
        }
        
        // 获取符号数据
        if index >= self.known_symbols.len() {
            return Err(crate::rds::error::RdsError::InvalidReference(index));
        }
        
        let symbol = &self.known_symbols[index];
        let ref_index = self.write_symbol(&symbol.name, symbol.encoding, writer)?;
        
        // 记录映射
        if index < self.known_symbol_mappings.len() {
            self.known_symbol_mappings[index] = Some(ref_index);
        }
        
        Ok(())
    }


    /// 写入环境
    ///
    /// 根据 EnvironmentIndex 写入对应的环境。
    /// 支持特殊环境（GLOBALENV_, BASEENV_, EMPTYENV_）的直接写入。
    ///
    /// # 参数
    /// * `obj` - 环境索引
    /// * `writer` - 写入目标
    pub fn write_environment<W: Write>(
        &mut self,
        obj: &EnvironmentIndex,
        writer: &mut W,
    ) -> Result<()> {
        // 处理特殊环境
        match obj.env_type {
            SEXPType::GlobalEnv => {
                write_header(SEXPType::GlobalEnv, writer)?;
                return Ok(());
            }
            SEXPType::BaseEnv => {
                write_header(SEXPType::BaseEnv, writer)?;
                return Ok(());
            }
            SEXPType::EmptyEnv => {
                write_header(SEXPType::EmptyEnv, writer)?;
                return Ok(());
            }
            _ => {}
        }
        
        let index = obj.index;
        
        // 检查是否已写入过
        if let Some(ref_index) = self.known_environment_mappings.get(index).copied().flatten() {
            write_ref_header(ref_index, writer)?;
            return Ok(());
        }
        
        // 获取环境数据
        if index >= self.known_environments.len() {
            return Err(crate::rds::error::RdsError::InvalidReference(index));
        }
        
        // 分配引用索引（在写入环境内容之前，以支持循环引用）
        let ref_index = self.allocate_reference();
        if index < self.known_environment_mappings.len() {
            self.known_environment_mappings[index] = Some(ref_index);
        }
        
        let env = &self.known_environments[index];
        
        // 写入环境头部
        write_header(SEXPType::Env, writer)?;
        
        // 写入锁定标志
        let locked_flag: i32 = if env.locked { 1 } else { 0 };
        writer.write_all(&locked_flag.to_be_bytes())?;
        
        // 写入父环境
        let parent_env = EnvironmentIndex {
            index: env.parent,
            env_type: env.parent_type,
        };
        self.write_environment(&parent_env, writer)?;
        
        // 写入变量（作为配对列表）
        self.write_environment_variables(env, writer)?;
        
        // 写入哈希表（如果有）
        if env.hashed {
            // 哈希环境：写入 NILVALUE 表示无哈希表
            write_header(SEXPType::NilValue, writer)?;
        } else {
            // 非哈希环境：写入 NILVALUE
            write_header(SEXPType::NilValue, writer)?;
        }
        
        // 写入属性（如果有）
        if !env.attributes.is_empty() {
            // 属性作为配对列表写入
            // 这里简化处理，实际需要调用 write_attributes
            write_header(SEXPType::NilValue, writer)?;
        } else {
            write_header(SEXPType::NilValue, writer)?;
        }
        
        Ok(())
    }

    /// 写入环境变量（作为配对列表）
    fn write_environment_variables<W: Write>(
        &mut self,
        env: &Environment,
        writer: &mut W,
    ) -> Result<()> {
        if env.variable_names.is_empty() {
            write_header(SEXPType::NilValue, writer)?;
            return Ok(());
        }
        
        // 环境变量作为配对列表写入
        // 每个变量是一个 (tag, value) 对
        for i in 0..env.variable_names.len() {
            // 写入 LIST 头部（带标签）
            let header = [0, 0, 0x04, SEXPType::List as u8]; // has_tag = true
            writer.write_all(&header)?;
            
            // 写入标签（符号）
            self.write_symbol(
                &env.variable_names[i],
                env.variable_encodings[i],
                writer,
            )?;
            
            // 写入值
            // 注意：这里需要调用 write_object，但为避免循环依赖，
            // 我们在这里只写入 Null 作为占位符
            // 实际实现需要在 object.rs 中处理
            write_header(SEXPType::NilValue, writer)?;
        }
        
        // 终止配对列表
        write_header(SEXPType::NilValue, writer)?;
        
        Ok(())
    }


    /// 写入外部指针
    ///
    /// 根据 ExternalPointerIndex 写入对应的外部指针。
    ///
    /// # 参数
    /// * `obj` - 外部指针索引
    /// * `writer` - 写入目标
    pub fn write_external_pointer<W: Write>(
        &mut self,
        obj: &ExternalPointerIndex,
        writer: &mut W,
    ) -> Result<()> {
        let index = obj.index;
        
        // 检查是否已写入过
        if let Some(ref_index) = self.known_external_pointer_mappings.get(index).copied().flatten() {
            write_ref_header(ref_index, writer)?;
            return Ok(());
        }
        
        // 获取外部指针数据
        if index >= self.known_external_pointers.len() {
            return Err(crate::rds::error::RdsError::InvalidReference(index));
        }
        
        // 分配引用索引
        let ref_index = self.allocate_reference();
        if index < self.known_external_pointer_mappings.len() {
            self.known_external_pointer_mappings[index] = Some(ref_index);
        }
        
        let ptr = &self.known_external_pointers[index];
        
        // 写入外部指针头部
        let has_attrs = !ptr.attributes.is_empty();
        let flags = if has_attrs { 0x02 } else { 0x00 };
        let header = [0, 0, flags, SEXPType::Extptr as u8];
        writer.write_all(&header)?;
        
        // 写入保护值
        // 注意：这里需要调用 write_object，但为避免循环依赖，
        // 我们在这里只写入 Null 作为占位符
        write_header(SEXPType::NilValue, writer)?;
        
        // 写入标签
        write_header(SEXPType::NilValue, writer)?;
        
        // 写入属性（如果有）
        if has_attrs {
            write_header(SEXPType::NilValue, writer)?;
        }
        
        Ok(())
    }

    /// 检查符号是否已注册
    pub fn has_symbol(&self, name: &str, encoding: StringEncoding) -> bool {
        self.symbol_mappings.contains_key(&(name.to_string(), encoding))
    }

    /// 获取符号的引用索引（如果已注册）
    pub fn get_symbol_ref(&self, name: &str, encoding: StringEncoding) -> Option<usize> {
        self.symbol_mappings.get(&(name.to_string(), encoding)).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_shared_write_info_new() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        // reference_count 从 1 开始，与 rds2cpp 保持一致
        assert_eq!(info.reference_count(), 1);
    }

    #[test]
    fn test_write_symbol_new() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        // reference_count 从 1 开始，第一个符号分配 1，然后 reference_count 变为 2
        let ref_index = info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        assert_eq!(ref_index, 1);
        assert_eq!(info.reference_count(), 2);
        
        // 验证写入的数据
        let data = buffer.into_inner();
        assert!(!data.is_empty());
        // 第一个字节应该是 SYM 类型头部
        assert_eq!(data[3], SEXPType::Sym as u8);
    }

    #[test]
    fn test_write_symbol_dedup() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        // 第一次写入（reference_count 从 1 开始，第一个符号分配 1）
        let ref1 = info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        assert_eq!(ref1, 1);
        
        // 第二次写入相同符号
        let ref2 = info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        assert_eq!(ref2, 1); // 应该返回相同的引用索引
        assert_eq!(info.reference_count(), 2); // 引用计数不应增加
    }

    #[test]
    fn test_write_symbol_different_encoding() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        // UTF-8 编码（reference_count 从 1 开始，第一个符号分配 1）
        let ref1 = info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        assert_eq!(ref1, 1);
        
        // Latin1 编码（不同编码应该是不同的符号，分配 2）
        let ref2 = info.write_symbol("test", StringEncoding::Latin1, &mut buffer).unwrap();
        assert_eq!(ref2, 2);
        assert_eq!(info.reference_count(), 3);
    }

    #[test]
    fn test_write_symbol_ref() {
        let symbols = vec![
            Symbol::new("x".to_string(), StringEncoding::Utf8),
            Symbol::new("y".to_string(), StringEncoding::Utf8),
        ];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        // 写入第一个符号引用（reference_count 从 1 开始，第一个符号分配 1）
        let idx = SymbolIndex { index: 0 };
        info.write_symbol_ref(&idx, &mut buffer).unwrap();
        assert_eq!(info.reference_count(), 2);
        
        // 再次写入相同符号引用
        info.write_symbol_ref(&idx, &mut buffer).unwrap();
        assert_eq!(info.reference_count(), 2); // 应该使用引用
    }

    #[test]
    fn test_write_global_environment() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        let env_idx = EnvironmentIndex {
            index: usize::MAX,
            env_type: SEXPType::GlobalEnv,
        };
        
        info.write_environment(&env_idx, &mut buffer).unwrap();
        
        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::GlobalEnv as u8);
    }

    #[test]
    fn test_write_base_environment() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        let env_idx = EnvironmentIndex {
            index: usize::MAX,
            env_type: SEXPType::BaseEnv,
        };
        
        info.write_environment(&env_idx, &mut buffer).unwrap();
        
        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::BaseEnv as u8);
    }

    #[test]
    fn test_write_empty_environment() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        let env_idx = EnvironmentIndex {
            index: usize::MAX,
            env_type: SEXPType::EmptyEnv,
        };
        
        info.write_environment(&env_idx, &mut buffer).unwrap();
        
        let data = buffer.into_inner();
        assert_eq!(data[3], SEXPType::EmptyEnv as u8);
    }

    #[test]
    fn test_has_symbol() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        assert!(!info.has_symbol("test", StringEncoding::Utf8));
        
        info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        
        assert!(info.has_symbol("test", StringEncoding::Utf8));
        assert!(!info.has_symbol("test", StringEncoding::Latin1));
    }

    #[test]
    fn test_get_symbol_ref() {
        let symbols: Vec<Symbol> = vec![];
        let environments: Vec<Environment> = vec![];
        let external_pointers: Vec<ExternalPointer> = vec![];
        
        let mut info = SharedWriteInfo::new(&symbols, &environments, &external_pointers);
        let mut buffer = Cursor::new(Vec::new());
        
        assert!(info.get_symbol_ref("test", StringEncoding::Utf8).is_none());
        
        info.write_symbol("test", StringEncoding::Utf8, &mut buffer).unwrap();
        
        // reference_count 从 1 开始，第一个符号分配 1
        assert_eq!(info.get_symbol_ref("test", StringEncoding::Utf8), Some(1));
    }
}
