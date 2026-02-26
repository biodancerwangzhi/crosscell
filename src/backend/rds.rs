// RDS backend implementation for CrossCell
// Wraps the RDS reader/writer to provide StorageBackend interface
//
// This module uses rds::RObject directly.

use super::{
    AttributeOp, DatasetOp, DynArray, GroupOp, ScalarType, SelectInfo, StorageBackend, StoreOp,
    Value,
};
use crate::error::{CrossCellError, Result};
use crate::rds::{parse_rds, write_rds as rds_write, GenericVector, RObject, RdsFile};
use std::path::{Path, PathBuf};

/// RDS backend
pub struct RDSBackend;

/// RDS store (represents an RDS file)
pub struct RDSStore {
    path: PathBuf,
    file: RdsFile,
    modified: bool,
}

/// RDS group (represents a list node in RDS)
pub struct RDSGroup {
    #[allow(dead_code)]
    path: Vec<String>,
    data: RObject,
}

/// RDS dataset (represents a vector in RDS)
pub struct RDSDataset {
    #[allow(dead_code)]
    path: Vec<String>,
    data: RObject,
}

impl StorageBackend for RDSBackend {
    const NAME: &'static str = "rds";

    type Store = RDSStore;
    type Group = RDSGroup;
    type Dataset = RDSDataset;

    fn new<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        // Create a new RDS file with an empty named list as root
        let mut file = RdsFile::default();
        file.object = RObject::GenericVector(GenericVector::default());

        Ok(RDSStore {
            path: path.as_ref().to_path_buf(),
            file,
            modified: true,
        })
    }

    fn open<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        let file = parse_rds(path.as_ref()).map_err(|e| {
            CrossCellError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                e.to_string(),
            ))
        })?;

        Ok(RDSStore {
            path: path.as_ref().to_path_buf(),
            file,
            modified: false,
        })
    }

    fn open_rw<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        // Same as open, but mark as modifiable
        let file = parse_rds(path.as_ref()).map_err(|e| {
            CrossCellError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                e.to_string(),
            ))
        })?;

        Ok(RDSStore {
            path: path.as_ref().to_path_buf(),
            file,
            modified: false,
        })
    }
}

impl StoreOp<RDSBackend> for RDSStore {
    fn filename(&self) -> PathBuf {
        self.path.clone()
    }

    fn close(self) -> Result<()> {
        if self.modified {
            rds_write(&self.file, &self.path).map_err(|e| {
                CrossCellError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    e.to_string(),
                ))
            })?;
        }
        Ok(())
    }
}

impl GroupOp<RDSBackend> for RDSStore {
    fn list(&self) -> Result<Vec<String>> {
        list_from_robject(&self.file.object)
    }

    fn new_group(&self, _name: &str) -> Result<RDSGroup> {
        // RDS groups are created on-demand
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support creating groups directly".to_string(),
        ))
    }

    fn open_group(&self, name: &str) -> Result<RDSGroup> {
        let data = get_item_from_robject(&self.file.object, name)?;
        Ok(RDSGroup {
            path: vec![name.to_string()],
            data: data.clone(),
        })
    }

    fn new_dataset(&self, _name: &str, _shape: &[usize], _dtype: ScalarType) -> Result<RDSDataset> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support creating datasets directly".to_string(),
        ))
    }

    fn open_dataset(&self, name: &str) -> Result<RDSDataset> {
        let data = get_item_from_robject(&self.file.object, name)?;
        Ok(RDSDataset {
            path: vec![name.to_string()],
            data: data.clone(),
        })
    }

    fn exists(&self, name: &str) -> Result<bool> {
        exists_in_robject(&self.file.object, name)
    }

    fn delete(&self, _name: &str) -> Result<()> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support deletion".to_string(),
        ))
    }
}

impl GroupOp<RDSBackend> for RDSGroup {
    fn list(&self) -> Result<Vec<String>> {
        list_from_robject(&self.data)
    }

    fn new_group(&self, _name: &str) -> Result<RDSGroup> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support creating groups".to_string(),
        ))
    }

    fn open_group(&self, name: &str) -> Result<RDSGroup> {
        let data = get_item_from_robject(&self.data, name)?;
        let mut path = self.path.clone();
        path.push(name.to_string());

        Ok(RDSGroup {
            path,
            data: data.clone(),
        })
    }

    fn new_dataset(&self, _name: &str, _shape: &[usize], _dtype: ScalarType) -> Result<RDSDataset> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support creating datasets".to_string(),
        ))
    }

    fn open_dataset(&self, name: &str) -> Result<RDSDataset> {
        let data = get_item_from_robject(&self.data, name)?;
        let mut path = self.path.clone();
        path.push(name.to_string());

        Ok(RDSDataset {
            path,
            data: data.clone(),
        })
    }

    fn exists(&self, name: &str) -> Result<bool> {
        exists_in_robject(&self.data, name)
    }

    fn delete(&self, _name: &str) -> Result<()> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support deletion".to_string(),
        ))
    }
}

impl DatasetOp<RDSBackend> for RDSDataset {
    fn dtype(&self) -> Result<ScalarType> {
        match &self.data {
            RObject::IntegerVector(_) => Ok(ScalarType::I32),
            RObject::DoubleVector(_) => Ok(ScalarType::F64),
            RObject::StringVector(_) => Ok(ScalarType::String),
            RObject::LogicalVector(_) => Ok(ScalarType::Bool),
            _ => Err(CrossCellError::UnsupportedType {
                type_name: format!("Unsupported RDS type: {}", self.data.type_name()),
            }),
        }
    }

    fn shape(&self) -> Vec<usize> {
        match &self.data {
            RObject::IntegerVector(vec) => vec![vec.data.len()],
            RObject::DoubleVector(vec) => vec![vec.data.len()],
            RObject::StringVector(vec) => vec![vec.data.len()],
            RObject::LogicalVector(vec) => vec![vec.data.len()],
            _ => vec![],
        }
    }

    fn read_slice(&self, _selection: &[SelectInfo]) -> Result<DynArray> {
        // For now, read all data
        match &self.data {
            RObject::IntegerVector(vec) => Ok(DynArray::I32(vec.data.clone())),
            RObject::DoubleVector(vec) => Ok(DynArray::F64(vec.data.clone())),
            RObject::StringVector(vec) => Ok(DynArray::String(vec.data.clone())),
            RObject::LogicalVector(vec) => {
                // Convert i32 to bool
                let bools: Vec<bool> = vec.data.iter().map(|&v| v != 0).collect();
                Ok(DynArray::Bool(bools))
            }
            _ => Err(CrossCellError::UnsupportedType {
                type_name: format!(
                    "Unsupported RDS type for reading: {}",
                    self.data.type_name()
                ),
            }),
        }
    }

    fn write_slice(&self, _data: &DynArray, _selection: &[SelectInfo]) -> Result<()> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support writing".to_string(),
        ))
    }
}

impl AttributeOp<RDSBackend> for RDSGroup {
    fn get_attr(&self, name: &str) -> Result<Value> {
        if let Some(attrs) = self.data.attributes() {
            if let Some(attr_value) = attrs.get(name) {
                robject_to_value(attr_value)
            } else {
                Err(CrossCellError::NotFound(format!(
                    "Attribute '{}' not found",
                    name
                )))
            }
        } else {
            Err(CrossCellError::NotFound("No attributes found".to_string()))
        }
    }

    fn set_attr(&mut self, _name: &str, _value: &Value) -> Result<()> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS backend does not support setting attributes".to_string(),
        ))
    }

    fn list_attrs(&self) -> Result<Vec<String>> {
        if let Some(attrs) = self.data.attributes() {
            Ok(attrs.names.clone())
        } else {
            Ok(Vec::new())
        }
    }

    fn has_attr(&self, name: &str) -> Result<bool> {
        if let Some(attrs) = self.data.attributes() {
            Ok(attrs.get(name).is_some())
        } else {
            Ok(false)
        }
    }
}

impl AttributeOp<RDSBackend> for RDSDataset {
    fn get_attr(&self, name: &str) -> Result<Value> {
        if let Some(attrs) = self.data.attributes() {
            if let Some(attr_value) = attrs.get(name) {
                robject_to_value(attr_value)
            } else {
                Err(CrossCellError::NotFound(format!(
                    "Attribute '{}' not found",
                    name
                )))
            }
        } else {
            Err(CrossCellError::NotFound("No attributes found".to_string()))
        }
    }

    fn set_attr(&mut self, _name: &str, _value: &Value) -> Result<()> {
        Err(CrossCellError::UnsupportedOperation(
            "RDS datasets do not support setting attributes".to_string(),
        ))
    }

    fn list_attrs(&self) -> Result<Vec<String>> {
        if let Some(attrs) = self.data.attributes() {
            Ok(attrs.names.clone())
        } else {
            Ok(Vec::new())
        }
    }

    fn has_attr(&self, name: &str) -> Result<bool> {
        if let Some(attrs) = self.data.attributes() {
            Ok(attrs.get(name).is_some())
        } else {
            Ok(false)
        }
    }
}

// Helper functions for working with RObject

/// List items from an RObject (List or PairList)
fn list_from_robject(obj: &RObject) -> Result<Vec<String>> {
    match obj {
        RObject::GenericVector(gv) => {
            if let Some(names) = gv.attributes.get_names() {
                Ok(names.to_vec())
            } else {
                // Generate numeric names for unnamed list
                Ok((0..gv.data.len()).map(|i| i.to_string()).collect())
            }
        }
        RObject::PairList(pl) => Ok(pl.tag_names.clone()),
        RObject::S4Object(s4) => {
            // For S4 objects, list the slot names
            Ok(s4.attributes.names.clone())
        }
        _ => Err(CrossCellError::InvalidFormat(format!(
            "Expected list, got {}",
            obj.type_name()
        ))),
    }
}

/// Get an item from an RObject by name or index
fn get_item_from_robject<'a>(obj: &'a RObject, name: &str) -> Result<&'a RObject> {
    match obj {
        RObject::GenericVector(gv) => {
            // First try to find by name
            if let Some(names) = gv.attributes.get_names() {
                if let Some(idx) = names.iter().position(|n| n == name) {
                    return gv.data.get(idx).ok_or_else(|| {
                        CrossCellError::NotFound(format!("Item '{}' not found", name))
                    });
                }
            }
            // Try to parse as index
            if let Ok(idx) = name.parse::<usize>() {
                if idx < gv.data.len() {
                    return Ok(&gv.data[idx]);
                }
            }
            Err(CrossCellError::NotFound(format!(
                "Item '{}' not found",
                name
            )))
        }
        RObject::PairList(pl) => {
            if let Some(idx) = pl.tag_names.iter().position(|n| n == name) {
                return pl
                    .data
                    .get(idx)
                    .ok_or_else(|| CrossCellError::NotFound(format!("Item '{}' not found", name)));
            }
            Err(CrossCellError::NotFound(format!(
                "Item '{}' not found",
                name
            )))
        }
        RObject::S4Object(s4) => {
            // For S4 objects, get the slot
            s4.attributes
                .get(name)
                .ok_or_else(|| CrossCellError::NotFound(format!("Slot '{}' not found", name)))
        }
        _ => Err(CrossCellError::InvalidFormat(format!(
            "Expected list, got {}",
            obj.type_name()
        ))),
    }
}

/// Check if an item exists in an RObject
fn exists_in_robject(obj: &RObject, name: &str) -> Result<bool> {
    match obj {
        RObject::GenericVector(gv) => {
            if let Some(names) = gv.attributes.get_names() {
                if names.iter().any(|n| n == name) {
                    return Ok(true);
                }
            }
            // Try to parse as index
            if let Ok(idx) = name.parse::<usize>() {
                return Ok(idx < gv.data.len());
            }
            Ok(false)
        }
        RObject::PairList(pl) => Ok(pl.tag_names.iter().any(|n| n == name)),
        RObject::S4Object(s4) => Ok(s4.attributes.get(name).is_some()),
        _ => Ok(false),
    }
}

/// Convert RObject to Value
fn robject_to_value(rval: &RObject) -> Result<Value> {
    match rval {
        RObject::IntegerVector(vec) if vec.data.len() == 1 => Ok(Value::I32(vec.data[0])),
        RObject::DoubleVector(vec) if vec.data.len() == 1 => Ok(Value::F64(vec.data[0])),
        RObject::StringVector(vec) if vec.data.len() == 1 => Ok(Value::String(vec.data[0].clone())),
        RObject::LogicalVector(vec) if vec.data.len() == 1 => Ok(Value::Bool(vec.data[0] != 0)),
        RObject::IntegerVector(vec) => Ok(Value::Array(DynArray::I32(vec.data.clone()))),
        RObject::DoubleVector(vec) => Ok(Value::Array(DynArray::F64(vec.data.clone()))),
        RObject::StringVector(vec) => Ok(Value::Array(DynArray::String(vec.data.clone()))),
        RObject::LogicalVector(vec) => {
            let bools: Vec<bool> = vec.data.iter().map(|&v| v != 0).collect();
            Ok(Value::Array(DynArray::Bool(bools)))
        }
        RObject::Null => {
            // Null is represented as an empty string array
            Ok(Value::Array(DynArray::String(Vec::new())))
        }
        _ => Err(CrossCellError::UnsupportedType {
            type_name: format!("Unsupported RObject type: {}", rval.type_name()),
        }),
    }
}
