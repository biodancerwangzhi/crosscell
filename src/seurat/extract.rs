//! Seurat 对象提取辅助函数

use crate::rds::{RObject, RdsFile};
use crate::seurat::error::{Result, SeuratError};

/// S4 对象中提取槽
pub fn extract_slot<'a>(s4_obj: &'a RObject, slot_name: &str, _file: &'a RdsFile) -> Result<&'a RObject> {
    match s4_obj {
        RObject::S4Object(s4) => {
            s4.attributes.get(slot_name)
                .ok_or_else(|| SeuratError::MissingSlot(slot_name.to_string()))
        }
        _ => Err(SeuratError::NotS4Object),
    }
}

/// S4 对象中提取可选槽
pub fn extract_slot_optional<'a>(s4_obj: &'a RObject, slot_name: &str, _file: &'a RdsFile) -> Result<Option<&'a RObject>> {
    match s4_obj {
        RObject::S4Object(s4) => Ok(s4.attributes.get(slot_name)),
        _ => Err(SeuratError::NotS4Object),
    }
}

/// 验证 S4 对象的类
pub fn verify_s4_class(s4_obj: &RObject, expected_class: &str) -> Result<()> {
    match s4_obj {
        RObject::S4Object(s4) => {
            if s4.class_name == expected_class {
                Ok(())
            } else {
                Err(SeuratError::WrongClass {
                    expected: expected_class.to_string(),
                    actual: vec![s4.class_name.clone()],
                })
            }
        }
        _ => Err(SeuratError::NotS4Object),
    }
}

/// 提取 active assay 名称
pub fn extract_active_assay(seurat: &RObject, file: &RdsFile) -> Result<String> {
    let active_assay = extract_slot(seurat, "active.assay", file)?;
    match active_assay {
        RObject::StringVector(sv) if !sv.data.is_empty() => Ok(sv.data[0].clone()),
        _ => Err(SeuratError::InvalidActiveAssay),
    }
}

/// 提取 project name
pub fn extract_project_name(seurat: &RObject, file: &RdsFile) -> Result<String> {
    let project_name = extract_slot(seurat, "project.name", file)?;
    match project_name {
        RObject::StringVector(sv) if !sv.data.is_empty() => Ok(sv.data[0].clone()),
        _ => Err(SeuratError::InvalidProjectName),
    }
}

/// 从列表中提取命名元素
pub fn extract_from_list<'a>(list: &'a RObject, name: &str, _file: &'a RdsFile) -> Result<&'a RObject> {
    match list {
        RObject::GenericVector(gv) => {
            // Check if it has names attribute
            if let Some(names) = gv.attributes.get_names() {
                for (i, n) in names.iter().enumerate() {
                    if n == name {
                        return gv.data.get(i)
                            .ok_or_else(|| SeuratError::MissingSlot(name.to_string()));
                    }
                }
            }
            Err(SeuratError::MissingSlot(name.to_string()))
        }
        _ => Err(SeuratError::InvalidAssaysStructure),
    }
}

/// 检RObject 是否S4 对象且类名匹
pub fn is_s4_class(obj: &RObject, class_name: &str) -> bool {
    match obj {
        RObject::S4Object(s4) => s4.class_name == class_name,
        _ => false,
    }
}

/// 获取 RObject 的类型名
pub fn type_name(obj: &RObject) -> &'static str {
    obj.type_name()
}

/// RObject 获取属
pub fn get_attribute<'a>(obj: &'a RObject, name: &str) -> Option<&'a RObject> {
    obj.attributes().and_then(|attrs| attrs.get(name))
}

/// 检查列表是否有 names 属性并返回命名列表迭代
/// 支持 List、PairList，以及包.Data slot S4 对象（如 Seurat Assays 容器
pub fn get_named_list_items<'a>(obj: &'a RObject) -> Option<Vec<(&'a str, &'a RObject)>> {
    match obj {
        RObject::GenericVector(gv) => {
            if let Some(names) = gv.attributes.get_names() {
                Some(names.iter()
                    .zip(gv.data.iter())
                    .map(|(n, v)| (n.as_str(), v))
                    .collect())
            } else {
                None
            }
        }
        RObject::PairList(pl) => {
            Some(pl.tag_names.iter()
                .zip(pl.data.iter())
                .map(|(n, v)| (n.as_str(), v))
                .collect())
        }
        // S4 对象可能包含 .Data slot（当 S4 类继承自 list 时）
        // 例如 Seurat assays slot 是一S4 Assays 容器
        RObject::S4Object(s4) => {
            // 尝试.Data slot 获取 named list
            if let Some(data_slot) = s4.attributes.get(".Data") {
                return get_named_list_items(data_slot);
            }
            // 尝试names + 属性中构建 named list
            // 有些 S4 对象直接attributes 中存储命名的子对
            if let Some(names_obj) = s4.attributes.get("names") {
                if let RObject::StringVector(sv) = names_obj {
                    // 查找对应.Data 列表
                    // 如果没有 .Data，尝试用 names 和其他属性匹
                    let names = &sv.data;
                    // 检查是否有names 对应的值存储在某个 list 属性中
                    for attr_name in &s4.attributes.names {
                        if let Some(list_obj) = s4.attributes.get(attr_name) {
                            if let RObject::GenericVector(gv) = list_obj {
                                if gv.data.len() == names.len() {
                                    return Some(names.iter()
                                        .zip(gv.data.iter())
                                        .map(|(n, v)| (n.as_str(), v))
                                        .collect());
                                }
                            }
                        }
                    }
                }
            }
            None
        }
        _ => None,
    }
}
