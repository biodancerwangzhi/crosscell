//! RDS basic tests

use crosscell::rds::{IntegerVector, DoubleVector, StringVector, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use tempfile::TempDir;

#[test]
fn test_rds_module_imports() {
    let _obj = RObject::Null;
    let _obj2 = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3], ..Default::default() });
    let _obj3 = RObject::DoubleVector(DoubleVector { data: vec![1.0, 2.0, 3.0], ..Default::default() });
    let mut sv = StringVector::default();
    sv.add("hello".to_string(), StringEncoding::Utf8);
    let _obj4 = RObject::StringVector(sv);
}

#[test]
fn test_write_integer_vector() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("test_int.rds");
    let obj = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3, 4, 5], ..Default::default() });
    write_rds(&file_path, &obj).expect("Failed to write RDS file");
    assert!(file_path.exists());
    assert!(std::fs::metadata(&file_path).unwrap().len() > 0);
}

#[test]
fn test_read_integer_vector() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("test_int.rds");
    let original = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3, 4, 5], ..Default::default() });
    write_rds(&file_path, &original).expect("Failed to write RDS file");
    let result = read_rds(&file_path).expect("Failed to read RDS file");
    match result {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data.len(), 5);
            assert_eq!(iv.data, vec![1, 2, 3, 4, 5]);
        }
        _ => panic!("Expected IntegerVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_roundtrip_integer() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("roundtrip_int.rds");
    let original = RObject::IntegerVector(IntegerVector { data: vec![10, 20, 30, 40, 50], ..Default::default() });
    write_rds(&file_path, &original).expect("Failed to write");
    let result = read_rds(&file_path).expect("Failed to read");
    match (&original, &result) {
        (RObject::IntegerVector(orig), RObject::IntegerVector(rt)) => {
            assert_eq!(orig.data, rt.data);
        }
        _ => panic!("Type mismatch in roundtrip"),
    }
}

#[test]
fn test_write_real_vector() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("test_real.rds");
    let obj = RObject::DoubleVector(DoubleVector { data: vec![1.1, 2.2, 3.3, 4.4, 5.5], ..Default::default() });
    write_rds(&file_path, &obj).expect("Failed to write RDS file");
    assert!(file_path.exists());
}

#[test]
fn test_roundtrip_real() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("roundtrip_real.rds");
    let original = RObject::DoubleVector(DoubleVector { data: vec![1.1, 2.2, 3.3], ..Default::default() });
    write_rds(&file_path, &original).expect("Failed to write");
    let result = read_rds(&file_path).expect("Failed to read");
    match (&original, &result) {
        (RObject::DoubleVector(orig), RObject::DoubleVector(rt)) => {
            assert_eq!(orig.data.len(), rt.data.len());
            for (a, b) in orig.data.iter().zip(rt.data.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
        }
        _ => panic!("Type mismatch in roundtrip"),
    }
}

#[test]
fn test_write_string_vector() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("test_string.rds");
    let mut sv = StringVector::default();
    sv.add("hello".to_string(), StringEncoding::Utf8);
    sv.add("world".to_string(), StringEncoding::Utf8);
    sv.add("rust".to_string(), StringEncoding::Utf8);
    let obj = RObject::StringVector(sv);
    write_rds(&file_path, &obj).expect("Failed to write RDS file");
    assert!(file_path.exists());
}

#[test]
fn test_roundtrip_string() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("roundtrip_string.rds");
    let mut sv = StringVector::default();
    sv.add("alpha".to_string(), StringEncoding::Utf8);
    sv.add("beta".to_string(), StringEncoding::Utf8);
    sv.add("gamma".to_string(), StringEncoding::Utf8);
    let original = RObject::StringVector(sv);
    write_rds(&file_path, &original).expect("Failed to write");
    let result = read_rds(&file_path).expect("Failed to read");
    match (&original, &result) {
        (RObject::StringVector(orig), RObject::StringVector(rt)) => {
            assert_eq!(orig.data, rt.data);
        }
        _ => panic!("Type mismatch in roundtrip"),
    }
}

#[test]
fn test_empty_vector() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("empty.rds");
    let original = RObject::IntegerVector(IntegerVector::default());
    write_rds(&file_path, &original).expect("Failed to write");
    let result = read_rds(&file_path).expect("Failed to read");
    match result {
        RObject::IntegerVector(iv) => assert_eq!(iv.data.len(), 0),
        _ => panic!("Expected IntegerVector"),
    }
}

#[test]
fn test_null_object() {
    let temp_dir = TempDir::new().unwrap();
    let file_path = temp_dir.path().join("null.rds");
    write_rds(&file_path, &RObject::Null).expect("Failed to write");
    let result = read_rds(&file_path).expect("Failed to read");
    assert!(matches!(result, RObject::Null));
}
