//! SingleCellExperiment 对象提取辅助函数

use crate::rds::{RObject, RdsFile};
use crate::sce::error::{Result, SceError};

/// SCE 相关的类名
pub const SCE_CLASS: &str = "SingleCellExperiment";
pub const SUMMARIZED_EXPERIMENT_CLASS: &str = "SummarizedExperiment";
pub const RANGED_SUMMARIZED_EXPERIMENT_CLASS: &str = "RangedSummarizedExperiment";

/// 检查对象是否是 SingleCellExperiment
pub fn is_sce_object(obj: &RObject) -> bool {
    match obj {
        RObject::S4Object(s4) => {
            s4.class_name == SCE_CLASS
                || s4.class_name == SUMMARIZED_EXPERIMENT_CLASS
                || s4.class_name == RANGED_SUMMARIZED_EXPERIMENT_CLASS
        }
        _ => false,
    }
}

/// 验证对象是 SingleCellExperiment
pub fn ensure_sce_class(obj: &RObject) -> Result<()> {
    match obj {
        RObject::S4Object(s4) => {
            if s4.class_name == SCE_CLASS
                || s4.class_name == SUMMARIZED_EXPERIMENT_CLASS
                || s4.class_name == RANGED_SUMMARIZED_EXPERIMENT_CLASS
            {
                Ok(())
            } else {
                Err(SceError::NotSceObject(s4.class_name.clone()))
            }
        }
        _ => Err(SceError::NotS4Object),
    }
}

/// 从 S4 对象中提取槽
pub fn extract_slot<'a>(s4_obj: &'a RObject, slot_name: &str, _file: &'a RdsFile) -> Result<&'a RObject> {
    match s4_obj {
        RObject::S4Object(s4) => {
            s4.attributes.get(slot_name)
                .ok_or_else(|| SceError::MissingSlot(slot_name.to_string()))
        }
        _ => Err(SceError::NotS4Object),
    }
}

/// 从 S4 对象中提取可选槽
pub fn extract_slot_optional<'a>(s4_obj: &'a RObject, slot_name: &str, _file: &'a RdsFile) -> Result<Option<&'a RObject>> {
    match s4_obj {
        RObject::S4Object(s4) => Ok(s4.attributes.get(slot_name)),
        _ => Err(SceError::NotS4Object),
    }
}

/// 获取 RObject 的类型名称
pub fn type_name(obj: &RObject) -> &'static str {
    obj.type_name()
}

/// 从 RObject 获取属性
pub fn get_attribute<'a>(obj: &'a RObject, name: &str) -> Option<&'a RObject> {
    obj.attributes().and_then(|attrs| attrs.get(name))
}

/// 检查列表是否有 names 属性并返回命名列表迭代器
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
        _ => None,
    }
}

/// 提取整数向量
pub fn extract_integer_vector(robj: &RObject) -> Result<Vec<i32>> {
    match robj {
        RObject::IntegerVector(iv) => Ok(iv.data.clone()),
        _ => Err(SceError::ParseError(format!(
            "Expected integer vector, got {}",
            robj.type_name()
        ))),
    }
}

/// 提取实数向量
pub fn extract_real_vector(robj: &RObject) -> Result<Vec<f64>> {
    match robj {
        RObject::DoubleVector(dv) => Ok(dv.data.clone()),
        _ => Err(SceError::ParseError(format!(
            "Expected real vector, got {}",
            robj.type_name()
        ))),
    }
}

/// 提取字符串向量
pub fn extract_string_vector(robj: &RObject) -> Result<Vec<String>> {
    match robj {
        RObject::StringVector(sv) => Ok(sv.data.clone()),
        _ => Err(SceError::ParseError(format!(
            "Expected string vector, got {}",
            robj.type_name()
        ))),
    }
}

/// 提取字符串标量
pub fn extract_string_scalar(robj: &RObject) -> Result<String> {
    match robj {
        RObject::StringVector(sv) => {
            if sv.data.is_empty() {
                Err(SceError::ParseError("Empty string vector".to_string()))
            } else {
                Ok(sv.data[0].clone())
            }
        }
        _ => Err(SceError::ParseError(format!(
            "Expected string, got {}",
            robj.type_name()
        ))),
    }
}
