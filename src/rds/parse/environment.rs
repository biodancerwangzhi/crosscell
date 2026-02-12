//! 环境解析
//!
//! 参考 rds2cpp 的 parse_environment.hpp 实现。
//!
//! ## 环境结构
//!
//! 环境在 RDS 中的结构：
//! 1. locked (4 bytes) - 是否锁定
//! 2. parent (header) - 父环境（可以是 REF, ENV, GLOBALENV_, BASEENV_, EMPTYENV_, NILVALUE_）
//! 3. frame/unhashed (header) - 未哈希的变量（LIST 或 NILVALUE_）
//! 4. hashtab (header) - 哈希表（如果 frame 是 NILVALUE_，则这里是 VEC）
//! 5. attributes (header) - 属性（LIST 或 NILVALUE_）

use std::io::Read;
use crate::rds::error::{Result, RdsError};
use crate::rds::sexp_type::SEXPType;
use crate::rds::r_object::{EnvironmentIndex, RObject};
use crate::rds::environment::Environment;
use super::utils::{read_header, read_i32};
use super::header::ParsedHeader;
use super::shared_info::SharedParseInfo;

fn parse_header_from_bytes(h: &[u8; 4]) -> Result<ParsedHeader> {
    ParsedHeader::from_raw(*h).ok_or_else(|| RdsError::UnsupportedType(h[3]))
}

/// 解析新环境体
/// 
/// 参考 rds2cpp parse_environment.hpp
pub fn parse_new_environment_body<R: Read, F>(
    reader: &mut R,
    _header: &ParsedHeader,
    shared: &mut SharedParseInfo,
    parse_fn: &mut F,
) -> Result<EnvironmentIndex>
where
    F: FnMut(&mut R, &mut SharedParseInfo) -> Result<RObject>,
{
    // 预先分配环境索引，以便内部引用有效
    let index = shared.request_environment();
    let mut env = Environment::new();
    
    // 1. 读取 locked 标志 (4 bytes)
    env.locked = read_i32(reader)? != 0;
    
    // 2. 读取父环境头部
    let parent_header = read_header(reader)?;
    let parent_type = parent_header[3];
    
    match parent_type {
        t if t == SEXPType::Ref as u8 => {
            // 引用到之前的环境
            env.parent = shared.get_environment_index(&parent_header)?;
            env.parent_type = SEXPType::Env;
        }
        t if t == SEXPType::Env as u8 => {
            // 嵌套的新环境
            let parsed = parse_header_from_bytes(&parent_header)?;
            let parent_env = parse_new_environment_body(reader, &parsed, shared, parse_fn)?;
            env.parent = parent_env.index;
            env.parent_type = SEXPType::Env;
        }
        t if t == SEXPType::GlobalEnv as u8 => {
            env.parent_type = SEXPType::GlobalEnv;
        }
        t if t == SEXPType::BaseEnv as u8 => {
            env.parent_type = SEXPType::BaseEnv;
        }
        t if t == SEXPType::EmptyEnv as u8 => {
            env.parent_type = SEXPType::EmptyEnv;
        }
        t if t == SEXPType::NilValue as u8 => {
            // 处理 igraph 等特殊情况
            env.parent_type = SEXPType::BaseEnv;
        }
        _ => {
            return Err(RdsError::ParseError {
                context: "environment".to_string(),
                message: format!("Could not resolve parent environment type: {}", parent_type),
            });
        }
    }
    
    // 3. 读取 frame/unhashed 头部
    let unhashed_header = read_header(reader)?;
    let unhashed_type = unhashed_header[3];
    
    if unhashed_type == SEXPType::List as u8 {
        // 未哈希的变量存储为 pairlist
        let parsed = parse_header_from_bytes(&unhashed_header)?;
        use super::pairlist::parse_pairlist_body;
        let plist = parse_pairlist_body(reader, &parsed, shared, parse_fn)?;
        
        // 提取变量
        for i in 0..plist.data.len() {
            if plist.has_tag[i] {
                env.variable_names.push(plist.tag_names[i].clone());
                env.variable_encodings.push(plist.tag_encodings[i]);
                env.variable_values.push(plist.data[i].clone());
            }
        }
        
        // 读取哈希表头部（应该是 NILVALUE_）
        let hashed_header = read_header(reader)?;
        if hashed_header[3] != SEXPType::NilValue as u8 {
            return Err(RdsError::ParseError {
                context: "environment".to_string(),
                message: format!("Unhashed environment should not contain a non-NULL hash table, got type {}", hashed_header[3]),
            });
        }
    } else if unhashed_type == SEXPType::NilValue as u8 {
        // 哈希环境 - 读取哈希表
        env.hashed = true;
        
        let hash_header = read_header(reader)?;
        if hash_header[3] != SEXPType::Vec as u8 {
            return Err(RdsError::ParseError {
                context: "environment".to_string(),
                message: format!("Environment's hash table should be a list (VEC), got type {}", hash_header[3]),
            });
        }
        
        // 解析哈希表（VEC）
        use super::list::parse_list_body;
        let vec = parse_list_body(reader, shared, parse_fn)?;
        
        // 从哈希表中提取变量
        for item in vec.data {
            match item {
                RObject::Null => continue,
                RObject::PairList(plist) => {
                    for j in 0..plist.data.len() {
                        if plist.has_tag[j] {
                            env.variable_names.push(plist.tag_names[j].clone());
                            env.variable_encodings.push(plist.tag_encodings[j]);
                            env.variable_values.push(plist.data[j].clone());
                        }
                    }
                }
                _ => {
                    // 忽略其他类型
                }
            }
        }
    } else {
        return Err(RdsError::ParseError {
            context: "environment".to_string(),
            message: format!("Environment's frame should be LIST or NILVALUE_, got type {}", unhashed_type),
        });
    }
    
    // 5. 读取属性头部
    let attr_header = read_header(reader)?;
    let attr_type = attr_header[3];
    
    if attr_type == SEXPType::List as u8 {
        // 有属性
        let parsed = parse_header_from_bytes(&attr_header)?;
        use super::attributes::parse_attributes_body;
        env.attributes = parse_attributes_body(reader, &parsed, shared, parse_fn)?;
    } else if attr_type != SEXPType::NilValue as u8 {
        return Err(RdsError::ParseError {
            context: "environment".to_string(),
            message: format!("Environment should be terminated by NULL, got type {}", attr_type),
        });
    }
    
    shared.update_environment(index, env);
    Ok(EnvironmentIndex { index, env_type: SEXPType::Env })
}

pub fn parse_global_environment_body() -> EnvironmentIndex {
    EnvironmentIndex { index: usize::MAX, env_type: SEXPType::GlobalEnv }
}

pub fn parse_base_environment_body() -> EnvironmentIndex {
    EnvironmentIndex { index: usize::MAX, env_type: SEXPType::BaseEnv }
}

pub fn parse_empty_environment_body() -> EnvironmentIndex {
    EnvironmentIndex { index: usize::MAX, env_type: SEXPType::EmptyEnv }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_global_env() {
        let e = parse_global_environment_body();
        assert_eq!(e.env_type, SEXPType::GlobalEnv);
    }

    #[test]
    fn test_base_env() {
        let e = parse_base_environment_body();
        assert_eq!(e.env_type, SEXPType::BaseEnv);
    }

    #[test]
    fn test_empty_env() {
        let e = parse_empty_environment_body();
        assert_eq!(e.env_type, SEXPType::EmptyEnv);
    }
}
