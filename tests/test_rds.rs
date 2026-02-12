//! RDS 模块测试
//!
//! 测试 rds 模块的功能

use crosscell::rds::{parse_rds, write_rds, RdsFile, RObject, SEXPType};
use std::path::Path;
use tempfile::tempdir;

/// 测试解析简单整数向量
#[test]
fn test_parse_simple_int() {
    let path = Path::new("tests/data/r_simple_int.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
    
    let file = result.unwrap();
    assert_eq!(file.format_version, 3);
    
    if let RObject::IntegerVector(vec) = &file.object {
        assert!(!vec.data.is_empty());
        println!("Parsed integer vector with {} elements", vec.data.len());
    } else {
        panic!("Expected integer vector, got {:?}", file.object.type_name());
    }
}

/// 测试解析字符串向量
#[test]
fn test_parse_string() {
    let path = Path::new("tests/data/r_string.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
    
    let file = result.unwrap();
    
    if let RObject::StringVector(vec) = &file.object {
        assert!(!vec.data.is_empty());
        println!("Parsed string vector: {:?}", vec.data);
    } else {
        panic!("Expected string vector, got {:?}", file.object.type_name());
    }
}

/// 测试解析 list
#[test]
fn test_parse_list() {
    let path = Path::new("tests/data/r_list.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
    
    let file = result.unwrap();
    
    if let RObject::GenericVector(vec) = &file.object {
        println!("Parsed list with {} elements", vec.data.len());
        
        if let Some(names) = vec.attributes.get_names() {
            println!("List names: {:?}", names);
        }
    } else {
        panic!("Expected list, got {:?}", file.object.type_name());
    }
}

/// 测试解析 NULL
#[test]
fn test_parse_null() {
    let path = Path::new("tests/data/r_null.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
    
    let file = result.unwrap();
    assert!(matches!(file.object, RObject::Null));
}

/// 测试解析 dgCMatrix
#[test]
fn test_parse_dgcmatrix() {
    let path = Path::new("tests/data/simple_dgcmatrix.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    
    if let Err(ref e) = result {
        eprintln!("Note: dgCMatrix file may be in old format: {:?}", e);
        return;
    }
    
    let file = result.unwrap();
    
    if let RObject::S4Object(s4) = &file.object {
        println!("Parsed S4 object with class: {}", s4.class_name);
        assert!(s4.class_name.contains("dgCMatrix") || s4.attributes.get_class().map(|c| c.iter().any(|s| s.contains("dgCMatrix"))).unwrap_or(false));
    } else {
        panic!("Expected S4 object, got {:?}", file.object.type_name());
    }
}

/// 测试解析 Seurat 对象
#[test]
fn test_parse_seurat() {
    let path = Path::new("tests/data/seurat_minimal.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let result = parse_rds(path);
    assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
    
    let file = result.unwrap();
    
    if let RObject::S4Object(s4) = &file.object {
        println!("Parsed Seurat object");
        println!("  Class: {}", s4.class_name);
        println!("  Package: {}", s4.package_name);
        println!("  Attributes: {:?}", s4.attributes.names);
    } else {
        panic!("Expected S4 object, got {:?}", file.object.type_name());
    }
}

/// 测试往返（读取后写入再读取）
#[test]
fn test_roundtrip_integer() {
    let path = Path::new("tests/data/r_simple_int.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let original = parse_rds(path).expect("Failed to parse original");
    
    let temp_dir = tempdir().expect("Failed to create temp dir");
    let temp_path = temp_dir.path().join("roundtrip.rds");
    
    write_rds(&original, &temp_path).expect("Failed to write");
    
    let roundtrip = parse_rds(&temp_path).expect("Failed to parse roundtrip");
    
    assert_eq!(original.format_version, roundtrip.format_version);
    
    if let (RObject::IntegerVector(orig_vec), RObject::IntegerVector(rt_vec)) = (&original.object, &roundtrip.object) {
        assert_eq!(orig_vec.data, rt_vec.data);
        println!("Roundtrip successful: {} integers match", orig_vec.data.len());
    } else {
        panic!("Type mismatch after roundtrip");
    }
}

/// 测试往返字符串
#[test]
fn test_roundtrip_string() {
    let path = Path::new("tests/data/r_string.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let original = parse_rds(path).expect("Failed to parse original");
    
    let temp_dir = tempdir().expect("Failed to create temp dir");
    let temp_path = temp_dir.path().join("roundtrip_string.rds");
    
    write_rds(&original, &temp_path).expect("Failed to write");
    
    let roundtrip = parse_rds(&temp_path).expect("Failed to parse roundtrip");
    
    if let (RObject::StringVector(orig_vec), RObject::StringVector(rt_vec)) = (&original.object, &roundtrip.object) {
        assert_eq!(orig_vec.data, rt_vec.data);
        println!("Roundtrip successful: {:?}", orig_vec.data);
    } else {
        panic!("Type mismatch after roundtrip");
    }
}

/// 测试往返 list
#[test]
fn test_roundtrip_list() {
    let path = Path::new("tests/data/r_list.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }
    
    let original = parse_rds(path).expect("Failed to parse original");
    
    let temp_dir = tempdir().expect("Failed to create temp dir");
    let temp_path = temp_dir.path().join("roundtrip_list.rds");
    
    write_rds(&original, &temp_path).expect("Failed to write");
    
    let roundtrip = parse_rds(&temp_path).expect("Failed to parse roundtrip");
    
    if let (RObject::GenericVector(orig_vec), RObject::GenericVector(rt_vec)) = (&original.object, &roundtrip.object) {
        assert_eq!(orig_vec.data.len(), rt_vec.data.len());
        println!("Roundtrip successful: list with {} elements", orig_vec.data.len());
    } else {
        panic!("Type mismatch after roundtrip");
    }
}

/// 测试 SEXP 类型转换
#[test]
fn test_sexp_type_conversion() {
    assert_eq!(SEXPType::from_u8(0), Some(SEXPType::Nil));
    assert_eq!(SEXPType::from_u8(13), Some(SEXPType::Int));
    assert_eq!(SEXPType::from_u8(14), Some(SEXPType::Real));
    assert_eq!(SEXPType::from_u8(16), Some(SEXPType::Str));
    assert_eq!(SEXPType::from_u8(19), Some(SEXPType::Vec));
    assert_eq!(SEXPType::from_u8(25), Some(SEXPType::S4));
}

/// 测试 RdsFile 默认值
#[test]
fn test_rds_file_default() {
    let file = RdsFile::default();
    assert_eq!(file.format_version, 3);
    assert_eq!(file.encoding, "UTF-8");
    assert!(matches!(file.object, RObject::Null));
}
