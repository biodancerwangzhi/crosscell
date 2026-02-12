// Memory backend implementation for CrossCell
// Pure in-memory storage for fast testing

use std::path::{Path, PathBuf};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use crate::error::{Result, CrossCellError};
use super::{
    StorageBackend, StoreOp, GroupOp, DatasetOp, AttributeOp,
    ScalarType, SelectInfo, DynArray, Value,
};

/// Memory backend
pub struct MemoryBackend;

/// Memory node (can be a group or dataset)
#[derive(Debug, Clone)]
enum MemoryNode {
    Group {
        children: HashMap<String, MemoryNode>,
        attributes: HashMap<String, Value>,
    },
    Dataset {
        data: DynArray,
        dtype: ScalarType,
        shape: Vec<usize>,
        attributes: HashMap<String, Value>,
    },
}

/// Memory store (root of the tree)
pub struct MemoryStore {
    path: PathBuf,
    root: Arc<RwLock<MemoryNode>>,
}

/// Memory group
pub struct MemoryGroup {
    #[allow(dead_code)]
    path: Vec<String>,
    node: Arc<RwLock<MemoryNode>>,
}

/// Memory dataset
pub struct MemoryDataset {
    #[allow(dead_code)]
    path: Vec<String>,
    node: Arc<RwLock<MemoryNode>>,
}

impl StorageBackend for MemoryBackend {
    const NAME: &'static str = "memory";
    
    type Store = MemoryStore;
    type Group = MemoryGroup;
    type Dataset = MemoryDataset;
    
    fn new<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        Ok(MemoryStore {
            path: path.as_ref().to_path_buf(),
            root: Arc::new(RwLock::new(MemoryNode::Group {
                children: HashMap::new(),
                attributes: HashMap::new(),
            })),
        })
    }
    
    fn open<P: AsRef<Path>>(_path: P) -> Result<Self::Store> {
        Err(CrossCellError::UnsupportedOperation("Memory backend cannot open existing files".to_string()))
    }
    
    fn open_rw<P: AsRef<Path>>(_path: P) -> Result<Self::Store> {
        Err(CrossCellError::UnsupportedOperation("Memory backend cannot open existing files".to_string()))
    }
}

impl StoreOp<MemoryBackend> for MemoryStore {
    fn filename(&self) -> PathBuf {
        self.path.clone()
    }
    
    fn close(self) -> Result<()> {
        // Nothing to do for memory backend
        Ok(())
    }
}

impl GroupOp<MemoryBackend> for MemoryStore {
    fn list(&self) -> Result<Vec<String>> {
        let node = self.root.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => Ok(children.keys().cloned().collect()),
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
    
    fn new_group(&self, name: &str) -> Result<MemoryGroup> {
        let mut node = self.root.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                let new_group = MemoryNode::Group {
                    children: HashMap::new(),
                    attributes: HashMap::new(),
                };
                let group_arc = Arc::new(RwLock::new(new_group.clone()));
                children.insert(name.to_string(), new_group);
                
                Ok(MemoryGroup {
                    path: vec![name.to_string()],
                    node: group_arc,
                })
            },
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
    
    fn open_group(&self, name: &str) -> Result<MemoryGroup> {
        let node = self.root.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => {
                if let Some(child) = children.get(name) {
                    match child {
                        MemoryNode::Group { .. } => {
                            Ok(MemoryGroup {
                                path: vec![name.to_string()],
                                node: Arc::new(RwLock::new(child.clone())),
                            })
                        },
                        _ => Err(CrossCellError::InvalidFormat(format!("'{}' is not a group", name))),
                    }
                } else {
                    Err(CrossCellError::NotFound(format!("Group '{}' not found", name)))
                }
            },
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
    
    fn new_dataset(&self, name: &str, shape: &[usize], dtype: ScalarType) -> Result<MemoryDataset> {
        let mut node = self.root.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                // Create empty data based on dtype
                let data = create_empty_array(dtype, shape)?;
                
                let new_dataset = MemoryNode::Dataset {
                    data: data.clone(),
                    dtype,
                    shape: shape.to_vec(),
                    attributes: HashMap::new(),
                };
                let dataset_arc = Arc::new(RwLock::new(new_dataset.clone()));
                children.insert(name.to_string(), new_dataset);
                
                Ok(MemoryDataset {
                    path: vec![name.to_string()],
                    node: dataset_arc,
                })
            },
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
    
    fn open_dataset(&self, name: &str) -> Result<MemoryDataset> {
        let node = self.root.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => {
                if let Some(child) = children.get(name) {
                    match child {
                        MemoryNode::Dataset { .. } => {
                            Ok(MemoryDataset {
                                path: vec![name.to_string()],
                                node: Arc::new(RwLock::new(child.clone())),
                            })
                        },
                        _ => Err(CrossCellError::InvalidFormat(format!("'{}' is not a dataset", name))),
                    }
                } else {
                    Err(CrossCellError::NotFound(format!("Dataset '{}' not found", name)))
                }
            },
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
    
    fn exists(&self, name: &str) -> Result<bool> {
        let node = self.root.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => Ok(children.contains_key(name)),
            _ => Ok(false),
        }
    }
    
    fn delete(&self, name: &str) -> Result<()> {
        let mut node = self.root.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                children.remove(name);
                Ok(())
            },
            _ => Err(CrossCellError::InvalidFormat("Root is not a group".to_string())),
        }
    }
}

impl GroupOp<MemoryBackend> for MemoryGroup {
    fn list(&self) -> Result<Vec<String>> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => Ok(children.keys().cloned().collect()),
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn new_group(&self, name: &str) -> Result<MemoryGroup> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                let new_group = MemoryNode::Group {
                    children: HashMap::new(),
                    attributes: HashMap::new(),
                };
                let group_arc = Arc::new(RwLock::new(new_group.clone()));
                children.insert(name.to_string(), new_group);
                
                let mut path = self.path.clone();
                path.push(name.to_string());
                
                Ok(MemoryGroup {
                    path,
                    node: group_arc,
                })
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn open_group(&self, name: &str) -> Result<MemoryGroup> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => {
                if let Some(child) = children.get(name) {
                    match child {
                        MemoryNode::Group { .. } => {
                            let mut path = self.path.clone();
                            path.push(name.to_string());
                            
                            Ok(MemoryGroup {
                                path,
                                node: Arc::new(RwLock::new(child.clone())),
                            })
                        },
                        _ => Err(CrossCellError::InvalidFormat(format!("'{}' is not a group", name))),
                    }
                } else {
                    Err(CrossCellError::NotFound(format!("Group '{}' not found", name)))
                }
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn new_dataset(&self, name: &str, shape: &[usize], dtype: ScalarType) -> Result<MemoryDataset> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                let data = create_empty_array(dtype, shape)?;
                
                let new_dataset = MemoryNode::Dataset {
                    data: data.clone(),
                    dtype,
                    shape: shape.to_vec(),
                    attributes: HashMap::new(),
                };
                let dataset_arc = Arc::new(RwLock::new(new_dataset.clone()));
                children.insert(name.to_string(), new_dataset);
                
                let mut path = self.path.clone();
                path.push(name.to_string());
                
                Ok(MemoryDataset {
                    path,
                    node: dataset_arc,
                })
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn open_dataset(&self, name: &str) -> Result<MemoryDataset> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => {
                if let Some(child) = children.get(name) {
                    match child {
                        MemoryNode::Dataset { .. } => {
                            let mut path = self.path.clone();
                            path.push(name.to_string());
                            
                            Ok(MemoryDataset {
                                path,
                                node: Arc::new(RwLock::new(child.clone())),
                            })
                        },
                        _ => Err(CrossCellError::InvalidFormat(format!("'{}' is not a dataset", name))),
                    }
                } else {
                    Err(CrossCellError::NotFound(format!("Dataset '{}' not found", name)))
                }
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn exists(&self, name: &str) -> Result<bool> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { children, .. } => Ok(children.contains_key(name)),
            _ => Ok(false),
        }
    }
    
    fn delete(&self, name: &str) -> Result<()> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Group { children, .. } => {
                children.remove(name);
                Ok(())
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
}

impl DatasetOp<MemoryBackend> for MemoryDataset {
    fn dtype(&self) -> Result<ScalarType> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { dtype, .. } => Ok(*dtype),
            _ => Err(CrossCellError::InvalidFormat("Not a dataset".to_string())),
        }
    }
    
    fn shape(&self) -> Vec<usize> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { shape, .. } => shape.clone(),
            _ => vec![],
        }
    }
    
    fn read_slice(&self, _selection: &[SelectInfo]) -> Result<DynArray> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { data, .. } => Ok(data.clone()),
            _ => Err(CrossCellError::InvalidFormat("Not a dataset".to_string())),
        }
    }
    
    fn write_slice(&self, data: &DynArray, _selection: &[SelectInfo]) -> Result<()> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Dataset { data: ref mut stored_data, .. } => {
                *stored_data = data.clone();
                Ok(())
            },
            _ => Err(CrossCellError::InvalidFormat("Not a dataset".to_string())),
        }
    }
}

impl AttributeOp<MemoryBackend> for MemoryGroup {
    fn get_attr(&self, name: &str) -> Result<Value> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { attributes, .. } => {
                attributes.get(name)
                    .cloned()
                    .ok_or_else(|| CrossCellError::NotFound(format!("Attribute '{}' not found", name)))
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn set_attr(&mut self, name: &str, value: &Value) -> Result<()> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Group { attributes, .. } => {
                attributes.insert(name.to_string(), value.clone());
                Ok(())
            },
            _ => Err(CrossCellError::InvalidFormat("Not a group".to_string())),
        }
    }
    
    fn list_attrs(&self) -> Result<Vec<String>> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { attributes, .. } => Ok(attributes.keys().cloned().collect()),
            _ => Ok(Vec::new()),
        }
    }
    
    fn has_attr(&self, name: &str) -> Result<bool> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Group { attributes, .. } => Ok(attributes.contains_key(name)),
            _ => Ok(false),
        }
    }
}

impl AttributeOp<MemoryBackend> for MemoryDataset {
    fn get_attr(&self, name: &str) -> Result<Value> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { attributes, .. } => {
                attributes.get(name)
                    .cloned()
                    .ok_or_else(|| CrossCellError::NotFound(format!("Attribute '{}' not found", name)))
            },
            _ => Err(CrossCellError::InvalidFormat("Not a dataset".to_string())),
        }
    }
    
    fn set_attr(&mut self, name: &str, value: &Value) -> Result<()> {
        let mut node = self.node.write().unwrap();
        match &mut *node {
            MemoryNode::Dataset { attributes, .. } => {
                attributes.insert(name.to_string(), value.clone());
                Ok(())
            },
            _ => Err(CrossCellError::InvalidFormat("Not a dataset".to_string())),
        }
    }
    
    fn list_attrs(&self) -> Result<Vec<String>> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { attributes, .. } => Ok(attributes.keys().cloned().collect()),
            _ => Ok(Vec::new()),
        }
    }
    
    fn has_attr(&self, name: &str) -> Result<bool> {
        let node = self.node.read().unwrap();
        match &*node {
            MemoryNode::Dataset { attributes, .. } => Ok(attributes.contains_key(name)),
            _ => Ok(false),
        }
    }
}

// Helper function to create empty array
fn create_empty_array(dtype: ScalarType, shape: &[usize]) -> Result<DynArray> {
    let total_size: usize = shape.iter().product();
    
    match dtype {
        ScalarType::I32 => Ok(DynArray::I32(vec![0; total_size])),
        ScalarType::I64 => Ok(DynArray::I64(vec![0; total_size])),
        ScalarType::F32 => Ok(DynArray::F32(vec![0.0; total_size])),
        ScalarType::F64 => Ok(DynArray::F64(vec![0.0; total_size])),
        ScalarType::Bool => Ok(DynArray::Bool(vec![false; total_size])),
        ScalarType::String => Ok(DynArray::String(vec![String::new(); total_size])),
        _ => Err(CrossCellError::UnsupportedType { type_name: format!("Unsupported dtype: {:?}", dtype) }),
    }
}
