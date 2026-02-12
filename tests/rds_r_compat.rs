//! RDS R 兼容性测试
//!
//! 测试 Rust 生成的 RDS 文件能否被 R 读取
//! 测试 R 生成的 RDS 文件能否被 Rust 读取

use crosscell::rds::{IntegerVector, DoubleVector, StringVector, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

#[test]
fn test_read_r_integer() {
    let path = Path::new("tests/data/r_integer.rds");
    if !path.exists() {
        eprintln!("Skipping: {} does not exist", path.display());
        return;
    }
    let result = read_rds(path).expect("Failed to read R integer file");
    match result {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data.len(), 5);
            assert_eq!(iv.data, vec![1, 2, 3, 4, 5]);
        }
        _ => panic!("Expected IntegerVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_read_r_real() {
    let path = Path::new("tests/data/r_real.rds");
    if !path.exists() { return; }
    let result = read_rds(path).expect("Failed to read R real file");
    match result {
        RObject::DoubleVector(dv) => {
            assert_eq!(dv.data.len(), 5);
            let expected = vec![1.1, 2.2, 3.3, 4.4, 5.5];
            for (a, b) in dv.data.iter().zip(expected.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
        }
        _ => panic!("Expected DoubleVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_read_r_string() {
    let path = Path::new("tests/data/r_string.rds");
    if !path.exists() { return; }
    let result = read_rds(path).expect("Failed to read R string file");
    match result {
        RObject::StringVector(sv) => {
            assert_eq!(sv.data.len(), 4);
            assert_eq!(sv.data[0], "hello");
            assert_eq!(sv.data[1], "world");
        }
        _ => panic!("Expected StringVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_read_r_logical() {
    let path = Path::new("tests/data/r_logical.rds");
    if !path.exists() { return; }
    let result = read_rds(path).expect("Failed to read R logical file");
    match result {
        RObject::LogicalVector(lv) => {
            assert_eq!(lv.data.len(), 4);
            // R logical: TRUE=1, FALSE=0
            assert_eq!(lv.data, vec![1, 0, 1, 0]);
        }
        _ => panic!("Expected LogicalVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_read_r_empty() {
    let path = Path::new("tests/data/r_empty_int.rds");
    if !path.exists() { return; }
    let result = read_rds(path).expect("Failed to read R empty file");
    match result {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data.len(), 0);
        }
        _ => panic!("Expected IntegerVector, got {:?}", result.type_name()),
    }
}

#[test]
fn test_read_r_null() {
    let path = Path::new("tests/data/r_null.rds");
    if !path.exists() { return; }
    let result = read_rds(path).expect("Failed to read R NULL file");
    assert!(matches!(result, RObject::Null));
}

#[test]
fn generate_rust_test_files() {
    use std::fs;
    fs::create_dir_all("tests/data").expect("Failed to create tests/data directory");

    let int_obj = RObject::IntegerVector(IntegerVector { data: vec![10, 20, 30, 40, 50], ..Default::default() });
    write_rds(Path::new("tests/data/rust_integer.rds"), &int_obj).expect("Failed to write");

    let real_obj = RObject::DoubleVector(DoubleVector { data: vec![1.5, 2.5, 3.5, 4.5, 5.5], ..Default::default() });
    write_rds(Path::new("tests/data/rust_real.rds"), &real_obj).expect("Failed to write");

    let mut sv = StringVector::default();
    sv.add("alpha".to_string(), StringEncoding::Utf8);
    sv.add("beta".to_string(), StringEncoding::Utf8);
    sv.add("gamma".to_string(), StringEncoding::Utf8);
    write_rds(Path::new("tests/data/rust_string.rds"), &RObject::StringVector(sv)).expect("Failed to write");

    let empty_obj = RObject::IntegerVector(IntegerVector::default());
    write_rds(Path::new("tests/data/rust_empty.rds"), &empty_obj).expect("Failed to write");

    write_rds(Path::new("tests/data/rust_null.rds"), &RObject::Null).expect("Failed to write");
}

#[test]
fn test_roundtrip_integer() {
    let r_path = Path::new("tests/data/r_integer.rds");
    let rust_path = Path::new("tests/data/roundtrip_int.rds");
    if !r_path.exists() { return; }

    let obj = read_rds(r_path).expect("Failed to read R file");
    write_rds(rust_path, &obj).expect("Failed to write Rust file");
    let obj2 = read_rds(rust_path).expect("Failed to read Rust file");

    match (&obj, &obj2) {
        (RObject::IntegerVector(v1), RObject::IntegerVector(v2)) => {
            assert_eq!(v1.data, v2.data);
        }
        _ => panic!("Type mismatch in roundtrip"),
    }
}

#[test]
fn test_roundtrip_real() {
    let r_path = Path::new("tests/data/r_real.rds");
    let rust_path = Path::new("tests/data/roundtrip_real.rds");
    if !r_path.exists() { return; }

    let obj = read_rds(r_path).expect("Failed to read");
    write_rds(rust_path, &obj).expect("Failed to write");
    let obj2 = read_rds(rust_path).expect("Failed to read");

    match (&obj, &obj2) {
        (RObject::DoubleVector(v1), RObject::DoubleVector(v2)) => {
            assert_eq!(v1.data.len(), v2.data.len());
            for (a, b) in v1.data.iter().zip(v2.data.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
        }
        _ => panic!("Type mismatch"),
    }
}

#[test]
fn test_roundtrip_string() {
    let r_path = Path::new("tests/data/r_string.rds");
    let rust_path = Path::new("tests/data/roundtrip_string.rds");
    if !r_path.exists() { return; }

    let obj = read_rds(r_path).expect("Failed to read");
    write_rds(rust_path, &obj).expect("Failed to write");
    let obj2 = read_rds(rust_path).expect("Failed to read");

    match (&obj, &obj2) {
        (RObject::StringVector(v1), RObject::StringVector(v2)) => {
            assert_eq!(v1.data, v2.data);
        }
        _ => panic!("Type mismatch"),
    }
}
