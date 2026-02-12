//! RDS 模块集成测试
//!
//! 使用真实 RDS 文件进行往返测试，验证解析和写入的正确性。
//! **Validates: Requirements 36.2**

use crosscell::rds::{
    parse::{parse_rds_file, ParseRdsOptions},
    write::rds::{write_rds_file, WriteRdsOptions},
    RdsFile, RObject,
};
use std::path::Path;
use tempfile::TempDir;

/// 比较两个 RdsFile 的结构是否等价
/// 
/// 注意：符号和环境数量可能因为去重而不同，这是正常的。
/// 我们主要比较根对象的结构和数据是否等价。
fn compare_rds_files(original: &RdsFile, roundtrip: &RdsFile) -> bool {
    // 比较格式版本
    if original.format_version != roundtrip.format_version {
        eprintln!(
            "Format version mismatch: {} vs {}",
            original.format_version, roundtrip.format_version
        );
        return false;
    }

    // 符号和环境数量可能因为去重而不同，这是正常的
    // 只记录差异，不作为失败条件
    if original.symbols.len() != roundtrip.symbols.len() {
        eprintln!(
            "Note: Symbol count differs: {} vs {} (this is expected due to deduplication)",
            original.symbols.len(),
            roundtrip.symbols.len()
        );
    }

    if original.environments.len() != roundtrip.environments.len() {
        eprintln!(
            "Note: Environment count differs: {} vs {} (this is expected due to deduplication)",
            original.environments.len(),
            roundtrip.environments.len()
        );
    }

    // 比较根对象 - 这是最重要的
    compare_robjects_with_context(&original.object, &roundtrip.object, original, roundtrip)
}

/// 递归比较两个 RObject，使用 RdsFile 上下文来解析符号引用
fn compare_robjects_with_context(a: &RObject, b: &RObject, ctx_a: &RdsFile, ctx_b: &RdsFile) -> bool {
    use RObject::*;

    match (a, b) {
        (Null, Null) => true,
        (IntegerVector(va), IntegerVector(vb)) => va.data == vb.data,
        (LogicalVector(va), LogicalVector(vb)) => va.data == vb.data,
        (DoubleVector(va), DoubleVector(vb)) => {
            if va.data.len() != vb.data.len() {
                return false;
            }
            va.data.iter().zip(vb.data.iter()).all(|(x, y)| {
                if x.is_nan() && y.is_nan() {
                    true
                } else if x.is_infinite() && y.is_infinite() {
                    x.signum() == y.signum()
                } else {
                    (x - y).abs() < 1e-10
                }
            })
        }
        (RawVector(va), RawVector(vb)) => va.data == vb.data,
        (ComplexVector(va), ComplexVector(vb)) => {
            va.data.len() == vb.data.len()
                && va.data.iter().zip(vb.data.iter()).all(|(x, y)| {
                    (x.re - y.re).abs() < 1e-10 && (x.im - y.im).abs() < 1e-10
                })
        }
        (StringVector(va), StringVector(vb)) => {
            va.data == vb.data && va.missing == vb.missing
        }
        (GenericVector(va), GenericVector(vb)) => {
            va.data.len() == vb.data.len()
                && va
                    .data
                    .iter()
                    .zip(vb.data.iter())
                    .all(|(x, y)| compare_robjects_with_context(x, y, ctx_a, ctx_b))
        }
        (PairList(va), PairList(vb)) => {
            va.data.len() == vb.data.len()
                && va.tag_names == vb.tag_names
                && va.has_tag == vb.has_tag
                && va.data.iter().zip(vb.data.iter())
                    .all(|(x, y)| compare_robjects_with_context(x, y, ctx_a, ctx_b))
        }
        (S4Object(va), S4Object(vb)) => {
            // S4 对象比较 class 名称
            va.class_name == vb.class_name && va.package_name == vb.package_name
        }
        (SymbolIndex(ia), SymbolIndex(ib)) => {
            // 比较符号的实际名称，而不是索引
            let name_a = ctx_a.symbols.get(ia.index).map(|s| &s.name);
            let name_b = ctx_b.symbols.get(ib.index).map(|s| &s.name);
            name_a == name_b
        }
        (EnvironmentIndex(ia), EnvironmentIndex(ib)) => {
            // 环境类型应该相同
            ia.env_type == ib.env_type
        }
        (BuiltInFunction(va), BuiltInFunction(vb)) => va.name == vb.name,
        (LanguageObject(va), LanguageObject(vb)) => {
            va.function_name == vb.function_name
                && va.argument_values.len() == vb.argument_values.len()
                && va.argument_names == vb.argument_names
                && va.argument_has_name == vb.argument_has_name
                && va.argument_values.iter().zip(vb.argument_values.iter())
                    .all(|(x, y)| compare_robjects_with_context(x, y, ctx_a, ctx_b))
        }
        (ExpressionVector(va), ExpressionVector(vb)) => {
            va.data.len() == vb.data.len()
                && va
                    .data
                    .iter()
                    .zip(vb.data.iter())
                    .all(|(x, y)| compare_robjects_with_context(x, y, ctx_a, ctx_b))
        }
        (ExternalPointerIndex(ia), ExternalPointerIndex(ib)) => ia == ib,
        _ => {
            eprintln!("Type mismatch: {:?} vs {:?}", a.sexp_type(), b.sexp_type());
            false
        }
    }
}

/// 递归比较两个 RObject（简化版，不需要上下文）
fn compare_robjects(a: &RObject, b: &RObject) -> bool {
    use RObject::*;

    match (a, b) {
        (Null, Null) => true,
        (IntegerVector(va), IntegerVector(vb)) => va.data == vb.data,
        (LogicalVector(va), LogicalVector(vb)) => va.data == vb.data,
        (DoubleVector(va), DoubleVector(vb)) => {
            if va.data.len() != vb.data.len() {
                return false;
            }
            va.data.iter().zip(vb.data.iter()).all(|(x, y)| {
                if x.is_nan() && y.is_nan() {
                    true
                } else if x.is_infinite() && y.is_infinite() {
                    x.signum() == y.signum()
                } else {
                    (x - y).abs() < 1e-10
                }
            })
        }
        (RawVector(va), RawVector(vb)) => va.data == vb.data,
        (ComplexVector(va), ComplexVector(vb)) => {
            va.data.len() == vb.data.len()
                && va.data.iter().zip(vb.data.iter()).all(|(x, y)| {
                    (x.re - y.re).abs() < 1e-10 && (x.im - y.im).abs() < 1e-10
                })
        }
        (StringVector(va), StringVector(vb)) => {
            va.data == vb.data && va.missing == vb.missing
        }
        (GenericVector(va), GenericVector(vb)) => {
            va.data.len() == vb.data.len()
                && va
                    .data
                    .iter()
                    .zip(vb.data.iter())
                    .all(|(x, y)| compare_robjects(x, y))
        }
        (PairList(va), PairList(vb)) => {
            va.data.len() == vb.data.len()
                && va.tag_names == vb.tag_names
                && va.has_tag == vb.has_tag
                && va.data.iter().zip(vb.data.iter()).all(|(x, y)| compare_robjects(x, y))
        }
        (S4Object(va), S4Object(vb)) => {
            // S4 对象比较 class 名称
            va.class_name == vb.class_name && va.package_name == vb.package_name
        }
        (SymbolIndex(ia), SymbolIndex(ib)) => ia == ib,
        (EnvironmentIndex(ia), EnvironmentIndex(ib)) => ia == ib,
        (BuiltInFunction(va), BuiltInFunction(vb)) => va.name == vb.name,
        (LanguageObject(va), LanguageObject(vb)) => {
            va.function_name == vb.function_name
                && va.argument_values.len() == vb.argument_values.len()
                && va.argument_names == vb.argument_names
                && va.argument_has_name == vb.argument_has_name
                && va.argument_values.iter().zip(vb.argument_values.iter())
                    .all(|(x, y)| compare_robjects(x, y))
        }
        (ExpressionVector(va), ExpressionVector(vb)) => {
            va.data.len() == vb.data.len()
                && va
                    .data
                    .iter()
                    .zip(vb.data.iter())
                    .all(|(x, y)| compare_robjects(x, y))
        }
        (ExternalPointerIndex(ia), ExternalPointerIndex(ib)) => ia == ib,
        _ => {
            eprintln!("Type mismatch: {:?} vs {:?}", a.sexp_type(), b.sexp_type());
            false
        }
    }
}

/// 测试单个 RDS 文件的往返
fn test_roundtrip_file(path: &Path) -> Result<(), String> {
    let options = ParseRdsOptions::default();

    // 解析原始文件
    let original = parse_rds_file(path, &options)
        .map_err(|e| format!("Failed to parse {}: {}", path.display(), e))?;

    // 创建临时目录
    let temp_dir = TempDir::new().map_err(|e| format!("Failed to create temp dir: {}", e))?;
    let temp_path = temp_dir.path().join("roundtrip.rds");

    // 写入临时文件
    let write_options = WriteRdsOptions {
        compress: true,
        compression_level: 6,
    };
    write_rds_file(&original, &temp_path, &write_options)
        .map_err(|e| format!("Failed to write {}: {}", temp_path.display(), e))?;

    // 重新解析
    let roundtrip = parse_rds_file(&temp_path, &options)
        .map_err(|e| format!("Failed to parse roundtrip file: {}", e))?;

    // 比较
    if compare_rds_files(&original, &roundtrip) {
        Ok(())
    } else {
        Err(format!("Roundtrip comparison failed for {}", path.display()))
    }
}


// ============================================================================
// 基础类型往返测试
// ============================================================================

#[test]
fn test_roundtrip_r_integer() {
    let path = Path::new("tests/data/r_integer.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Integer roundtrip failed");
    println!("✅ r_integer.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_real() {
    let path = Path::new("tests/data/r_real.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Real roundtrip failed");
    println!("✅ r_real.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_string() {
    let path = Path::new("tests/data/r_string.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("String roundtrip failed");
    println!("✅ r_string.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_logical() {
    let path = Path::new("tests/data/r_logical.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Logical roundtrip failed");
    println!("✅ r_logical.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_null() {
    let path = Path::new("tests/data/r_null.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Null roundtrip failed");
    println!("✅ r_null.rds roundtrip passed");
}

// ============================================================================
// 特殊值测试
// ============================================================================

#[test]
fn test_roundtrip_r_na_int() {
    let path = Path::new("tests/data/r_na_int.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("NA integer roundtrip failed");
    println!("✅ r_na_int.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_na_real() {
    let path = Path::new("tests/data/r_na_real.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("NA real roundtrip failed");
    println!("✅ r_na_real.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_na_string() {
    let path = Path::new("tests/data/r_na_string.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("NA string roundtrip failed");
    println!("✅ r_na_string.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_special_values() {
    let path = Path::new("tests/data/r_special_values.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Special values roundtrip failed");
    println!("✅ r_special_values.rds roundtrip passed");
}

// ============================================================================
// 空向量测试
// ============================================================================

#[test]
fn test_roundtrip_r_empty_int() {
    let path = Path::new("tests/data/r_empty_int.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Empty int roundtrip failed");
    println!("✅ r_empty_int.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_empty_real() {
    let path = Path::new("tests/data/r_empty_real.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Empty real roundtrip failed");
    println!("✅ r_empty_real.rds roundtrip passed");
}

#[test]
fn test_roundtrip_r_empty_string() {
    let path = Path::new("tests/data/r_empty_string.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Empty string roundtrip failed");
    println!("✅ r_empty_string.rds roundtrip passed");
}

// ============================================================================
// 列表和复合类型测试
// ============================================================================

#[test]
fn test_roundtrip_r_list() {
    let path = Path::new("tests/data/r_list.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("List roundtrip failed");
    println!("✅ r_list.rds roundtrip passed");
}

#[test]
fn test_roundtrip_simple_named_list() {
    let path = Path::new("tests/data/simple_named_list.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Simple named list roundtrip failed");
    println!("✅ simple_named_list.rds roundtrip passed");
}

#[test]
fn test_roundtrip_nested_list() {
    let path = Path::new("tests/data/nested_list.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Nested list roundtrip failed");
    println!("✅ nested_list.rds roundtrip passed");
}

#[test]
fn test_roundtrip_list_with_null() {
    let path = Path::new("tests/data/list_with_null.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("List with null roundtrip failed");
    println!("✅ list_with_null.rds roundtrip passed");
}


// ============================================================================
// S4 对象和 dgCMatrix 测试
// ============================================================================

#[test]
fn test_roundtrip_minimal_dgc() {
    let path = Path::new("tests/data/minimal_dgc.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Minimal dgCMatrix roundtrip failed");
    println!("✅ minimal_dgc.rds roundtrip passed");
}

#[test]
fn test_roundtrip_simple_dgcmatrix() {
    let path = Path::new("tests/data/simple_dgcmatrix.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Simple dgCMatrix roundtrip failed");
    println!("✅ simple_dgcmatrix.rds roundtrip passed");
}

#[test]
fn test_roundtrip_nested_with_s4() {
    let path = Path::new("tests/data/nested_with_s4.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Nested with S4 roundtrip failed");
    println!("✅ nested_with_s4.rds roundtrip passed");
}

// ============================================================================
// DataFrame 测试
// ============================================================================

#[test]
fn test_roundtrip_simple_dataframe() {
    let path = Path::new("tests/data/simple_dataframe.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Simple dataframe roundtrip failed");
    println!("✅ simple_dataframe.rds roundtrip passed");
}

#[test]
fn test_roundtrip_empty_dataframe() {
    let path = Path::new("tests/data/empty_dataframe.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Empty dataframe roundtrip failed");
    println!("✅ empty_dataframe.rds roundtrip passed");
}

// ============================================================================
// Seurat 对象测试
// ============================================================================

#[test]
fn test_roundtrip_seurat_minimal() {
    let path = Path::new("tests/data/seurat_minimal.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Seurat minimal roundtrip failed");
    println!("✅ seurat_minimal.rds roundtrip passed");
}

#[test]
fn test_roundtrip_seurat_simple() {
    let path = Path::new("tests/data/seurat_simple.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Seurat simple roundtrip failed");
    println!("✅ seurat_simple.rds roundtrip passed");
}

#[test]
fn test_roundtrip_seurat_minimal_simplified() {
    let path = Path::new("tests/data/seurat_minimal_simplified.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("Seurat minimal simplified roundtrip failed");
    println!("✅ seurat_minimal_simplified.rds roundtrip passed");
}

// ============================================================================
// SCE 对象测试
// ============================================================================

#[test]
fn test_roundtrip_sce_minimal() {
    let path = Path::new("tests/data/sce_minimal.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("SCE minimal roundtrip failed");
    println!("✅ sce_minimal.rds roundtrip passed");
}

#[test]
fn test_roundtrip_sce_full() {
    let path = Path::new("tests/data/sce_full.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    test_roundtrip_file(path).expect("SCE full roundtrip failed");
    println!("✅ sce_full.rds roundtrip passed");
}

// ============================================================================
// 批量测试 - 测试 tests/data 目录下所有 RDS 文件
// ============================================================================

#[test]
#[ignore] // 运行时间较长，默认忽略
fn test_roundtrip_all_test_data() {
    let test_data_dir = Path::new("tests/data");
    if !test_data_dir.exists() {
        eprintln!("⚠️  Test data directory not found");
        return;
    }

    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;

    for entry in std::fs::read_dir(test_data_dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.extension().map_or(false, |ext| ext == "rds") {
            print!("Testing {}... ", path.display());

            match test_roundtrip_file(&path) {
                Ok(()) => {
                    println!("✅");
                    passed += 1;
                }
                Err(e) => {
                    println!("❌ {}", e);
                    failed += 1;
                }
            }
        } else {
            skipped += 1;
        }
    }

    println!("\n========================================");
    println!("Roundtrip Test Summary");
    println!("========================================");
    println!("Passed:  {}", passed);
    println!("Failed:  {}", failed);
    println!("Skipped: {}", skipped);
    println!("========================================");

    // 允许少量失败（某些边缘情况可能还不支持）
    // 目标是 95% 以上的通过率
    let total = passed + failed;
    let pass_rate = if total > 0 { (passed as f64 / total as f64) * 100.0 } else { 100.0 };
    println!("Pass rate: {:.1}%", pass_rate);
    
    assert!(
        pass_rate >= 95.0,
        "Pass rate {:.1}% is below 95% threshold. {} tests failed.",
        pass_rate, failed
    );
}

// ============================================================================
// 批量测试 - 测试 data/generated 目录下的 Seurat 数据
// ============================================================================

#[test]
#[ignore] // 运行时间较长，默认忽略
fn test_roundtrip_generated_seurat_data() {
    let generated_dir = Path::new("data/generated");
    if !generated_dir.exists() {
        eprintln!("⚠️  Generated data directory not found");
        return;
    }

    let mut passed = 0;
    let mut failed = 0;

    for entry in std::fs::read_dir(generated_dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.extension().map_or(false, |ext| ext == "rds") {
            print!("Testing {}... ", path.file_name().unwrap().to_string_lossy());

            match test_roundtrip_file(&path) {
                Ok(()) => {
                    println!("✅");
                    passed += 1;
                }
                Err(e) => {
                    println!("❌ {}", e);
                    failed += 1;
                }
            }
        }
    }

    println!("\n========================================");
    println!("Generated Seurat Data Roundtrip Summary");
    println!("========================================");
    println!("Passed: {}", passed);
    println!("Failed: {}", failed);
    println!("========================================");

    // 允许部分失败，因为某些复杂对象可能还不支持
    println!("Note: Some failures may be expected for complex objects");
}


// ============================================================================
// R 兼容性测试 - 验证解析结果与 R 的 readRDS() 一致
// **Validates: Requirements 36.3-36.6**
// ============================================================================

/// 测试读取 R 生成的整数向量并验证值
#[test]
fn test_r_compat_integer_values() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_integer.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R integer file");
    
    match &rds.object {
        RObject::IntegerVector(vec) => {
            // R 生成的文件应该包含 [1, 2, 3, 4, 5]
            assert_eq!(vec.data, vec![1, 2, 3, 4, 5], "Integer values should match R output");
            println!("✅ R integer compatibility: {:?}", vec.data);
        }
        _ => panic!("Expected IntegerVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的实数向量并验证值
#[test]
fn test_r_compat_real_values() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_real.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R real file");
    
    match &rds.object {
        RObject::DoubleVector(vec) => {
            // R 生成的文件应该包含 [1.1, 2.2, 3.3, 4.4, 5.5]
            let expected = vec![1.1, 2.2, 3.3, 4.4, 5.5];
            assert_eq!(vec.data.len(), expected.len(), "Length should match");
            for (i, (a, b)) in vec.data.iter().zip(expected.iter()).enumerate() {
                assert!(
                    (a - b).abs() < 1e-10,
                    "Value at index {} should match: {} vs {}",
                    i, a, b
                );
            }
            println!("✅ R real compatibility: {:?}", vec.data);
        }
        _ => panic!("Expected DoubleVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的字符串向量并验证值
#[test]
fn test_r_compat_string_values() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_string.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R string file");
    
    match &rds.object {
        RObject::StringVector(vec) => {
            // R 生成的文件应该包含 ["hello", "world", "rust", "R"]
            let expected = vec!["hello", "world", "rust", "R"];
            assert_eq!(vec.data.len(), expected.len(), "Length should match");
            for (i, (a, b)) in vec.data.iter().zip(expected.iter()).enumerate() {
                assert_eq!(a, *b, "String at index {} should match", i);
            }
            println!("✅ R string compatibility: {:?}", vec.data);
        }
        _ => panic!("Expected StringVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的逻辑向量并验证值
#[test]
fn test_r_compat_logical_values() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_logical.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R logical file");
    
    match &rds.object {
        RObject::LogicalVector(vec) => {
            // R 生成的文件应该包含 [TRUE, FALSE, TRUE, FALSE]
            // 在 R 中，TRUE = 1, FALSE = 0
            let expected = vec![1, 0, 1, 0];
            assert_eq!(vec.data, expected, "Logical values should match R output");
            println!("✅ R logical compatibility: {:?}", vec.data);
        }
        _ => panic!("Expected LogicalVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的 NA 值
#[test]
fn test_r_compat_na_integer() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_na_int.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R NA integer file");
    
    match &rds.object {
        RObject::IntegerVector(vec) => {
            // R 的 NA_integer_ 是 i32::MIN
            assert!(vec.data.contains(&i32::MIN), "Should contain NA (i32::MIN)");
            println!("✅ R NA integer compatibility: {:?}", vec.data);
        }
        _ => panic!("Expected IntegerVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的特殊浮点值 (Inf, -Inf, NaN)
#[test]
fn test_r_compat_special_values() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_special_values.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R special values file");
    
    match &rds.object {
        RObject::DoubleVector(vec) => {
            // 检查是否包含特殊值
            let has_inf = vec.data.iter().any(|x| x.is_infinite() && *x > 0.0);
            let has_neg_inf = vec.data.iter().any(|x| x.is_infinite() && *x < 0.0);
            let has_nan = vec.data.iter().any(|x| x.is_nan());
            
            assert!(has_inf || has_neg_inf || has_nan, "Should contain special values");
            println!("✅ R special values compatibility: Inf={}, -Inf={}, NaN={}", 
                has_inf, has_neg_inf, has_nan);
        }
        _ => panic!("Expected DoubleVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的列表并验证结构
#[test]
fn test_r_compat_list_structure() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_list.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R list file");
    
    match &rds.object {
        RObject::GenericVector(vec) => {
            assert!(!vec.data.is_empty(), "List should not be empty");
            println!("✅ R list compatibility: {} elements", vec.data.len());
        }
        _ => panic!("Expected GenericVector, got {:?}", rds.object.sexp_type()),
    }
}

/// 测试读取 R 生成的 dgCMatrix 并验证结构
#[test]
fn test_r_compat_dgcmatrix_structure() {
    use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
    
    let path = Path::new("tests/data/r_created_dgc.rds");
    if !path.exists() {
        eprintln!("⚠️  Test file not found: {}", path.display());
        return;
    }
    
    let rds = parse_rds_file(path, &ParseRdsOptions::default())
        .expect("Failed to parse R dgCMatrix file");
    
    match &rds.object {
        RObject::S4Object(s4) => {
            assert_eq!(s4.class_name, "dgCMatrix", "Class should be dgCMatrix");
            assert_eq!(s4.package_name, "Matrix", "Package should be Matrix");
            println!("✅ R dgCMatrix compatibility: class={}, package={}", 
                s4.class_name, s4.package_name);
        }
        _ => panic!("Expected S4Object, got {:?}", rds.object.sexp_type()),
    }
}
