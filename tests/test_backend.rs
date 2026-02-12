// Tests for backend abstraction layer

use crosscell::backend::{
    StorageBackend, StoreOp, GroupOp, DatasetOp, AttributeOp,
    ScalarType, DynArray, Value,
    memory::MemoryBackend,
};

#[test]
fn test_memory_backend_basic() {
    // Create a new memory store
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create a group
    let group = store.new_group("test_group").unwrap();
    
    // Verify group exists
    assert!(store.exists("test_group").unwrap());
    
    // List should contain the group
    let items = store.list().unwrap();
    assert_eq!(items.len(), 1);
    assert_eq!(items[0], "test_group");
    
    // Close the store
    store.close().unwrap();
}

#[test]
fn test_memory_backend_dataset() {
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create a dataset
    let dataset = store.new_dataset("test_data", &[10], ScalarType::I32).unwrap();
    
    // Check dtype and shape
    assert_eq!(dataset.dtype().unwrap(), ScalarType::I32);
    assert_eq!(dataset.shape(), vec![10]);
    
    // Write data
    let data = DynArray::I32(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    dataset.write_all(&data).unwrap();
    
    // Read data back
    let read_data = dataset.read_all().unwrap();
    match read_data {
        DynArray::I32(vec) => {
            assert_eq!(vec, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        },
        _ => panic!("Expected I32 array"),
    }
    
    store.close().unwrap();
}

#[test]
fn test_memory_backend_attributes() {
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create a group
    let mut group = store.new_group("test_group").unwrap();
    
    // Set attributes
    group.set_attr("version", &Value::I32(1)).unwrap();
    group.set_attr("name", &Value::String("test".to_string())).unwrap();
    group.set_attr("pi", &Value::F64(3.14159)).unwrap();
    
    // Check attributes exist
    assert!(group.has_attr("version").unwrap());
    assert!(group.has_attr("name").unwrap());
    assert!(group.has_attr("pi").unwrap());
    
    // List attributes
    let attrs = group.list_attrs().unwrap();
    assert_eq!(attrs.len(), 3);
    
    // Get attributes
    match group.get_attr("version").unwrap() {
        Value::I32(v) => assert_eq!(v, 1),
        _ => panic!("Expected I32"),
    }
    
    match group.get_attr("name").unwrap() {
        Value::String(s) => assert_eq!(s, "test"),
        _ => panic!("Expected String"),
    }
    
    match group.get_attr("pi").unwrap() {
        Value::F64(v) => assert!((v - 3.14159).abs() < 1e-6),
        _ => panic!("Expected F64"),
    }
    
    store.close().unwrap();
}

#[test]
fn test_memory_backend_nested_groups() {
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create nested groups
    let group1 = store.new_group("level1").unwrap();
    let group2 = group1.new_group("level2").unwrap();
    let _group3 = group2.new_group("level3").unwrap();
    
    // Verify structure
    assert!(store.exists("level1").unwrap());
    assert!(group1.exists("level2").unwrap());
    assert!(group2.exists("level3").unwrap());
    
    store.close().unwrap();
}

#[test]
fn test_memory_backend_multiple_datasets() {
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create multiple datasets with different types
    let ds1 = store.new_dataset("int_data", &[5], ScalarType::I32).unwrap();
    let ds2 = store.new_dataset("float_data", &[5], ScalarType::F64).unwrap();
    let ds3 = store.new_dataset("bool_data", &[5], ScalarType::Bool).unwrap();
    
    // Write data
    ds1.write_all(&DynArray::I32(vec![1, 2, 3, 4, 5])).unwrap();
    ds2.write_all(&DynArray::F64(vec![1.1, 2.2, 3.3, 4.4, 5.5])).unwrap();
    ds3.write_all(&DynArray::Bool(vec![true, false, true, false, true])).unwrap();
    
    // Verify all datasets exist
    assert!(store.exists("int_data").unwrap());
    assert!(store.exists("float_data").unwrap());
    assert!(store.exists("bool_data").unwrap());
    
    // List datasets
    let items = store.list().unwrap();
    assert_eq!(items.len(), 3);
    
    store.close().unwrap();
}

#[test]
fn test_memory_backend_delete() {
    let store = MemoryBackend::new("test.mem").unwrap();
    
    // Create items
    store.new_group("group1").unwrap();
    store.new_dataset("data1", &[10], ScalarType::I32).unwrap();
    
    // Verify they exist
    assert!(store.exists("group1").unwrap());
    assert!(store.exists("data1").unwrap());
    
    // Delete one
    store.delete("group1").unwrap();
    
    // Verify deletion
    assert!(!store.exists("group1").unwrap());
    assert!(store.exists("data1").unwrap());
    
    store.close().unwrap();
}
