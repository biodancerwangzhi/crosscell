// HDF5 backend implementation for CrossCell
// Wraps the hdf5 crate to provide StorageBackend interface

use super::{
    AttributeOp, DatasetOp, DynArray, GroupOp, ScalarType, SelectInfo, StorageBackend, StoreOp,
    Value,
};
use crate::error::{CrossCellError, Result};
use hdf5::{Dataset as H5Dataset, File, Group as H5Group};
use std::path::{Path, PathBuf};

/// HDF5 backend
pub struct H5Backend;

/// HDF5 store (file handle)
pub struct H5Store(File);

/// HDF5 group
pub struct H5GroupWrapper(H5Group);

/// HDF5 dataset
pub struct H5DatasetWrapper(H5Dataset);

impl StorageBackend for H5Backend {
    const NAME: &'static str = "hdf5";

    type Store = H5Store;
    type Group = H5GroupWrapper;
    type Dataset = H5DatasetWrapper;

    fn new<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        let file = File::create(path.as_ref())
            .map_err(|e| CrossCellError::Internal(format!("Failed to create HDF5 file: {}", e)))?;
        Ok(H5Store(file))
    }

    fn open<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        let file = File::open(path.as_ref())
            .map_err(|e| CrossCellError::Internal(format!("Failed to open HDF5 file: {}", e)))?;
        Ok(H5Store(file))
    }

    fn open_rw<P: AsRef<Path>>(path: P) -> Result<Self::Store> {
        let file = File::open_rw(path.as_ref()).map_err(|e| {
            CrossCellError::Internal(format!("Failed to open HDF5 file in RW mode: {}", e))
        })?;
        Ok(H5Store(file))
    }
}

impl StoreOp<H5Backend> for H5Store {
    fn filename(&self) -> PathBuf {
        PathBuf::from(self.0.filename())
    }

    fn close(self) -> Result<()> {
        self.0
            .close()
            .map_err(|e| CrossCellError::Internal(format!("Failed to close HDF5 file: {}", e)))
    }
}

impl GroupOp<H5Backend> for H5Store {
    fn list(&self) -> Result<Vec<String>> {
        let mut items = Vec::new();
        for name in self
            .0
            .member_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list group members: {}", e)))?
        {
            items.push(name);
        }
        Ok(items)
    }

    fn new_group(&self, name: &str) -> Result<H5GroupWrapper> {
        let group = self.0.create_group(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to create group '{}': {}", name, e))
        })?;
        Ok(H5GroupWrapper(group))
    }

    fn open_group(&self, name: &str) -> Result<H5GroupWrapper> {
        let group = self.0.group(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to open group '{}': {}", name, e))
        })?;
        Ok(H5GroupWrapper(group))
    }

    fn new_dataset(
        &self,
        name: &str,
        shape: &[usize],
        dtype: ScalarType,
    ) -> Result<H5DatasetWrapper> {
        // Convert shape to hdf5 format
        let h5_shape: Vec<usize> = shape.to_vec();

        // Create dataset based on dtype
        let dataset = match dtype {
            ScalarType::I32 => self.0.new_dataset::<i32>().shape(h5_shape).create(name),
            ScalarType::I64 => self.0.new_dataset::<i64>().shape(h5_shape).create(name),
            ScalarType::F32 => self.0.new_dataset::<f32>().shape(h5_shape).create(name),
            ScalarType::F64 => self.0.new_dataset::<f64>().shape(h5_shape).create(name),
            ScalarType::Bool => self.0.new_dataset::<bool>().shape(h5_shape).create(name),
            _ => {
                return Err(CrossCellError::UnsupportedType {
                    type_name: format!("Unsupported dtype: {:?}", dtype),
                })
            }
        }
        .map_err(|e| {
            CrossCellError::Internal(format!("Failed to create dataset '{}': {}", name, e))
        })?;

        Ok(H5DatasetWrapper(dataset))
    }

    fn open_dataset(&self, name: &str) -> Result<H5DatasetWrapper> {
        let dataset = self.0.dataset(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to open dataset '{}': {}", name, e))
        })?;
        Ok(H5DatasetWrapper(dataset))
    }

    fn exists(&self, name: &str) -> Result<bool> {
        Ok(self.0.link_exists(name))
    }

    fn delete(&self, name: &str) -> Result<()> {
        self.0
            .unlink(name)
            .map_err(|e| CrossCellError::Internal(format!("Failed to delete '{}': {}", name, e)))
    }
}

impl GroupOp<H5Backend> for H5GroupWrapper {
    fn list(&self) -> Result<Vec<String>> {
        let mut items = Vec::new();
        for name in self
            .0
            .member_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list group members: {}", e)))?
        {
            items.push(name);
        }
        Ok(items)
    }

    fn new_group(&self, name: &str) -> Result<H5GroupWrapper> {
        let group = self.0.create_group(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to create group '{}': {}", name, e))
        })?;
        Ok(H5GroupWrapper(group))
    }

    fn open_group(&self, name: &str) -> Result<H5GroupWrapper> {
        let group = self.0.group(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to open group '{}': {}", name, e))
        })?;
        Ok(H5GroupWrapper(group))
    }

    fn new_dataset(
        &self,
        name: &str,
        shape: &[usize],
        dtype: ScalarType,
    ) -> Result<H5DatasetWrapper> {
        let h5_shape: Vec<usize> = shape.to_vec();

        let dataset = match dtype {
            ScalarType::I32 => self.0.new_dataset::<i32>().shape(h5_shape).create(name),
            ScalarType::I64 => self.0.new_dataset::<i64>().shape(h5_shape).create(name),
            ScalarType::F32 => self.0.new_dataset::<f32>().shape(h5_shape).create(name),
            ScalarType::F64 => self.0.new_dataset::<f64>().shape(h5_shape).create(name),
            ScalarType::Bool => self.0.new_dataset::<bool>().shape(h5_shape).create(name),
            _ => {
                return Err(CrossCellError::UnsupportedType {
                    type_name: format!("Unsupported dtype: {:?}", dtype),
                })
            }
        }
        .map_err(|e| {
            CrossCellError::Internal(format!("Failed to create dataset '{}': {}", name, e))
        })?;

        Ok(H5DatasetWrapper(dataset))
    }

    fn open_dataset(&self, name: &str) -> Result<H5DatasetWrapper> {
        let dataset = self.0.dataset(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to open dataset '{}': {}", name, e))
        })?;
        Ok(H5DatasetWrapper(dataset))
    }

    fn exists(&self, name: &str) -> Result<bool> {
        Ok(self.0.link_exists(name))
    }

    fn delete(&self, name: &str) -> Result<()> {
        self.0
            .unlink(name)
            .map_err(|e| CrossCellError::Internal(format!("Failed to delete '{}': {}", name, e)))
    }
}

impl DatasetOp<H5Backend> for H5DatasetWrapper {
    fn dtype(&self) -> Result<ScalarType> {
        let dtype = self
            .0
            .dtype()
            .map_err(|e| CrossCellError::Internal(format!("Failed to get dtype: {}", e)))?;

        // Map HDF5 dtype to ScalarType
        if dtype.is::<i32>() {
            Ok(ScalarType::I32)
        } else if dtype.is::<i64>() {
            Ok(ScalarType::I64)
        } else if dtype.is::<f32>() {
            Ok(ScalarType::F32)
        } else if dtype.is::<f64>() {
            Ok(ScalarType::F64)
        } else if dtype.is::<bool>() {
            Ok(ScalarType::Bool)
        } else {
            Err(CrossCellError::UnsupportedType {
                type_name: "Unsupported HDF5 dtype".to_string(),
            })
        }
    }

    fn shape(&self) -> Vec<usize> {
        self.0.shape()
    }

    fn read_slice(&self, _selection: &[SelectInfo]) -> Result<DynArray> {
        // For now, read all data (slicing will be implemented later)
        let dtype = self.dtype()?;

        match dtype {
            ScalarType::I32 => {
                let data: Vec<i32> = self
                    .0
                    .read_raw()
                    .map_err(|e| CrossCellError::Internal(format!("Failed to read data: {}", e)))?;
                Ok(DynArray::I32(data))
            }
            ScalarType::I64 => {
                let data: Vec<i64> = self
                    .0
                    .read_raw()
                    .map_err(|e| CrossCellError::Internal(format!("Failed to read data: {}", e)))?;
                Ok(DynArray::I64(data))
            }
            ScalarType::F32 => {
                let data: Vec<f32> = self
                    .0
                    .read_raw()
                    .map_err(|e| CrossCellError::Internal(format!("Failed to read data: {}", e)))?;
                Ok(DynArray::F32(data))
            }
            ScalarType::F64 => {
                let data: Vec<f64> = self
                    .0
                    .read_raw()
                    .map_err(|e| CrossCellError::Internal(format!("Failed to read data: {}", e)))?;
                Ok(DynArray::F64(data))
            }
            ScalarType::Bool => {
                let data: Vec<bool> = self
                    .0
                    .read_raw()
                    .map_err(|e| CrossCellError::Internal(format!("Failed to read data: {}", e)))?;
                Ok(DynArray::Bool(data))
            }
            _ => Err(CrossCellError::UnsupportedType {
                type_name: format!("Unsupported dtype: {:?}", dtype),
            }),
        }
    }

    fn write_slice(&self, data: &DynArray, _selection: &[SelectInfo]) -> Result<()> {
        // For now, write all data (slicing will be implemented later)
        match data {
            DynArray::I32(vec) => {
                self.0.write(vec).map_err(|e| {
                    CrossCellError::Internal(format!("Failed to write data: {}", e))
                })?;
            }
            DynArray::I64(vec) => {
                self.0.write(vec).map_err(|e| {
                    CrossCellError::Internal(format!("Failed to write data: {}", e))
                })?;
            }
            DynArray::F32(vec) => {
                self.0.write(vec).map_err(|e| {
                    CrossCellError::Internal(format!("Failed to write data: {}", e))
                })?;
            }
            DynArray::F64(vec) => {
                self.0.write(vec).map_err(|e| {
                    CrossCellError::Internal(format!("Failed to write data: {}", e))
                })?;
            }
            DynArray::Bool(vec) => {
                self.0.write(vec).map_err(|e| {
                    CrossCellError::Internal(format!("Failed to write data: {}", e))
                })?;
            }
            _ => {
                return Err(CrossCellError::UnsupportedType {
                    type_name: "Unsupported data type".to_string(),
                })
            }
        }
        Ok(())
    }
}

impl AttributeOp<H5Backend> for H5GroupWrapper {
    fn get_attr(&self, name: &str) -> Result<Value> {
        let attr = self.0.attr(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to get attribute '{}': {}", name, e))
        })?;

        // Try to read as different types
        if let Ok(val) = attr.read_scalar::<i32>() {
            return Ok(Value::I32(val));
        }
        if let Ok(val) = attr.read_scalar::<i64>() {
            return Ok(Value::I64(val));
        }
        if let Ok(val) = attr.read_scalar::<f32>() {
            return Ok(Value::F32(val));
        }
        if let Ok(val) = attr.read_scalar::<f64>() {
            return Ok(Value::F64(val));
        }
        if let Ok(val) = attr.read_scalar::<bool>() {
            return Ok(Value::Bool(val));
        }

        Err(CrossCellError::UnsupportedType {
            type_name: format!("Unsupported attribute type for '{}'", name),
        })
    }

    fn set_attr(&mut self, name: &str, value: &Value) -> Result<()> {
        match value {
            Value::I32(v) => {
                self.0
                    .new_attr::<i32>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::I64(v) => {
                self.0
                    .new_attr::<i64>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::F32(v) => {
                self.0
                    .new_attr::<f32>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::F64(v) => {
                self.0
                    .new_attr::<f64>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::Bool(v) => {
                self.0
                    .new_attr::<bool>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            _ => {
                return Err(CrossCellError::UnsupportedType {
                    type_name: "Unsupported attribute type".to_string(),
                })
            }
        }
        Ok(())
    }

    fn list_attrs(&self) -> Result<Vec<String>> {
        let mut attrs = Vec::new();
        for name in self
            .0
            .attr_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list attributes: {}", e)))?
        {
            attrs.push(name);
        }
        Ok(attrs)
    }

    fn has_attr(&self, name: &str) -> Result<bool> {
        Ok(self
            .0
            .attr_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list attributes: {}", e)))?
            .contains(&name.to_string()))
    }
}

impl AttributeOp<H5Backend> for H5DatasetWrapper {
    fn get_attr(&self, name: &str) -> Result<Value> {
        let attr = self.0.attr(name).map_err(|e| {
            CrossCellError::Internal(format!("Failed to get attribute '{}': {}", name, e))
        })?;

        if let Ok(val) = attr.read_scalar::<i32>() {
            return Ok(Value::I32(val));
        }
        if let Ok(val) = attr.read_scalar::<i64>() {
            return Ok(Value::I64(val));
        }
        if let Ok(val) = attr.read_scalar::<f32>() {
            return Ok(Value::F32(val));
        }
        if let Ok(val) = attr.read_scalar::<f64>() {
            return Ok(Value::F64(val));
        }
        if let Ok(val) = attr.read_scalar::<bool>() {
            return Ok(Value::Bool(val));
        }

        Err(CrossCellError::UnsupportedType {
            type_name: format!("Unsupported attribute type for '{}'", name),
        })
    }

    fn set_attr(&mut self, name: &str, value: &Value) -> Result<()> {
        match value {
            Value::I32(v) => {
                self.0
                    .new_attr::<i32>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::I64(v) => {
                self.0
                    .new_attr::<i64>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::F32(v) => {
                self.0
                    .new_attr::<f32>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::F64(v) => {
                self.0
                    .new_attr::<f64>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            Value::Bool(v) => {
                self.0
                    .new_attr::<bool>()
                    .create(name)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to create attribute: {}", e))
                    })?
                    .write_scalar(v)
                    .map_err(|e| {
                        CrossCellError::Internal(format!("Failed to write attribute: {}", e))
                    })?;
            }
            _ => {
                return Err(CrossCellError::UnsupportedType {
                    type_name: "Unsupported attribute type".to_string(),
                })
            }
        }
        Ok(())
    }

    fn list_attrs(&self) -> Result<Vec<String>> {
        let mut attrs = Vec::new();
        for name in self
            .0
            .attr_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list attributes: {}", e)))?
        {
            attrs.push(name);
        }
        Ok(attrs)
    }

    fn has_attr(&self, name: &str) -> Result<bool> {
        Ok(self
            .0
            .attr_names()
            .map_err(|e| CrossCellError::Internal(format!("Failed to list attributes: {}", e)))?
            .contains(&name.to_string()))
    }
}
