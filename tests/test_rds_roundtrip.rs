//! RDS 往返测试：读取 → 写入 → 重新读取
//!
//! 验证 RDS 写入模块能正确写出可被重新解析的文件。
//! 使用 data/generated 目录下的真实 Seurat 数据集。

use crosscell::rds::{
    parse::{parse_rds_file, ParseRdsOptions},
    write::rds::{write_rds_file, WriteRdsOptions},
    RObject,
};
use std::path::Path;
use std::time::Instant;

/// 递归比较两个 RObject 的结构和数据
/// 返回 (match, description)
fn compare_robjects(a: &RObject, b: &RObject, path: &str) -> Vec<String> {
    let mut diffs = Vec::new();

    match (a, b) {
        (RObject::Null, RObject::Null) => {}

        (RObject::IntegerVector(va), RObject::IntegerVector(vb)) => {
            if va.data != vb.data {
                diffs.push(format!(
                    "{}: IntegerVector data mismatch ({} vs {} elements)",
                    path,
                    va.data.len(),
                    vb.data.len()
                ));
            }
        }
        (RObject::LogicalVector(va), RObject::LogicalVector(vb)) => {
            if va.data != vb.data {
                diffs.push(format!("{}: LogicalVector data mismatch", path));
            }
        }
        (RObject::DoubleVector(va), RObject::DoubleVector(vb)) => {
            if va.data.len() != vb.data.len() {
                diffs.push(format!(
                    "{}: DoubleVector length mismatch ({} vs {})",
                    path,
                    va.data.len(),
                    vb.data.len()
                ));
            } else {
                for (i, (x, y)) in va.data.iter().zip(vb.data.iter()).enumerate() {
                    if x.is_nan() && y.is_nan() {
                        continue;
                    }
                    if x != y {
                        diffs.push(format!(
                            "{}: DoubleVector[{}] mismatch ({} vs {})",
                            path, i, x, y
                        ));
                        break;
                    }
                }
            }
        }
        (RObject::RawVector(va), RObject::RawVector(vb)) => {
            if va.data != vb.data {
                diffs.push(format!("{}: RawVector data mismatch", path));
            }
        }
        (RObject::ComplexVector(va), RObject::ComplexVector(vb)) => {
            if va.data.len() != vb.data.len() {
                diffs.push(format!("{}: ComplexVector length mismatch", path));
            }
        }
        (RObject::StringVector(va), RObject::StringVector(vb)) => {
            if va.data != vb.data {
                let max_show = 3;
                let mismatches: Vec<_> = va
                    .data
                    .iter()
                    .zip(vb.data.iter())
                    .enumerate()
                    .filter(|(_, (a, b))| a != b)
                    .take(max_show)
                    .map(|(i, (a, b))| format!("[{}]: {:?} vs {:?}", i, a, b))
                    .collect();
                diffs.push(format!(
                    "{}: StringVector data mismatch ({} vs {} elements, first diffs: {})",
                    path,
                    va.data.len(),
                    vb.data.len(),
                    mismatches.join(", ")
                ));
            }
            if va.missing != vb.missing {
                diffs.push(format!("{}: StringVector missing flags mismatch", path));
            }
        }

        (RObject::GenericVector(va), RObject::GenericVector(vb)) => {
            if va.data.len() != vb.data.len() {
                diffs.push(format!(
                    "{}: GenericVector length mismatch ({} vs {})",
                    path,
                    va.data.len(),
                    vb.data.len()
                ));
            } else {
                for (i, (ca, cb)) in va.data.iter().zip(vb.data.iter()).enumerate() {
                    let child_path = format!("{}[{}]", path, i);
                    diffs.extend(compare_robjects(ca, cb, &child_path));
                    if diffs.len() > 20 {
                        break;
                    } // 限制差异数量
                }
            }
        }

        (RObject::PairList(va), RObject::PairList(vb)) => {
            if va.data.len() != vb.data.len() {
                diffs.push(format!(
                    "{}: PairList length mismatch ({} vs {})",
                    path,
                    va.data.len(),
                    vb.data.len()
                ));
            }
        }

        (RObject::S4Object(sa), RObject::S4Object(sb)) => {
            if sa.class_name != sb.class_name {
                diffs.push(format!(
                    "{}: S4 class mismatch ({} vs {})",
                    path, sa.class_name, sb.class_name
                ));
            }
            if sa.package_name != sb.package_name {
                diffs.push(format!(
                    "{}: S4 package mismatch ({} vs {})",
                    path, sa.package_name, sb.package_name
                ));
            }
        }

        (RObject::SymbolIndex(_), RObject::SymbolIndex(_)) => {}
        (RObject::EnvironmentIndex(_), RObject::EnvironmentIndex(_)) => {}
        (RObject::BuiltInFunction(fa), RObject::BuiltInFunction(fb)) => {
            if fa.name != fb.name {
                diffs.push(format!(
                    "{}: BuiltIn name mismatch ({} vs {})",
                    path, fa.name, fb.name
                ));
            }
        }
        (RObject::LanguageObject(_), RObject::LanguageObject(_)) => {}
        (RObject::ExpressionVector(_), RObject::ExpressionVector(_)) => {}
        (RObject::ExternalPointerIndex(_), RObject::ExternalPointerIndex(_)) => {}

        _ => {
            diffs.push(format!(
                "{}: type mismatch ({} vs {})",
                path,
                type_name(a),
                type_name(b)
            ));
        }
    }

    diffs
}

fn type_name(obj: &RObject) -> &'static str {
    match obj {
        RObject::Null => "Null",
        RObject::IntegerVector(_) => "IntegerVector",
        RObject::LogicalVector(_) => "LogicalVector",
        RObject::DoubleVector(_) => "DoubleVector",
        RObject::RawVector(_) => "RawVector",
        RObject::ComplexVector(_) => "ComplexVector",
        RObject::StringVector(_) => "StringVector",
        RObject::GenericVector(_) => "GenericVector",
        RObject::PairList(_) => "PairList",
        RObject::S4Object(_) => "S4Object",
        RObject::SymbolIndex(_) => "SymbolIndex",
        RObject::EnvironmentIndex(_) => "EnvironmentIndex",
        RObject::BuiltInFunction(_) => "BuiltInFunction",
        RObject::LanguageObject(_) => "LanguageObject",
        RObject::ExpressionVector(_) => "ExpressionVector",
        RObject::ExternalPointerIndex(_) => "ExternalPointerIndex",
    }
}

/// 对单个文件执行往返测试
fn roundtrip_single_file(input_path: &Path) -> Result<Vec<String>, String> {
    let filename = input_path.file_name().unwrap().to_string_lossy();
    let output_path = input_path.with_extension("roundtrip.rds");

    // 1. 读取原始文件
    let options = ParseRdsOptions::default();
    let original = parse_rds_file(input_path, &options).map_err(|e| format!("读取失败: {}", e))?;

    // 2. 写入到新文件
    let write_opts = WriteRdsOptions {
        compress: true,
        compression_level: 6,
    };
    write_rds_file(&original, &output_path, &write_opts).map_err(|e| format!("写入失败: {}", e))?;

    // 3. 重新读取写入的文件
    let reloaded =
        parse_rds_file(&output_path, &options).map_err(|e| format!("重新读取失败: {}", e))?;

    // 4. 比较
    let diffs = compare_robjects(&original.object, &reloaded.object, &filename);

    // 清理临时文件
    let _ = std::fs::remove_file(&output_path);

    Ok(diffs)
}

#[test]
fn test_roundtrip_all_generated_rds() {
    let generated_dir = Path::new("data/generated");
    if !generated_dir.exists() {
        eprintln!("⚠️  data/generated 目录不存在，跳过往返测试");
        return;
    }

    println!("\n========================================");
    println!("🔄 RDS 往返测试 (读取 → 写入 → 重新读取)");
    println!("========================================\n");

    let mut passed = 0;
    let mut write_failed = 0;
    let mut diff_failed = 0;
    let mut errors: Vec<(String, String)> = Vec::new();

    let total_start = Instant::now();

    let mut rds_files: Vec<_> = std::fs::read_dir(generated_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "rds"))
        .collect();
    rds_files.sort_by_key(|e| e.file_name());

    println!("找到 {} 个 RDS 文件\n", rds_files.len());

    for entry in rds_files {
        let path = entry.path();
        let filename = path.file_name().unwrap().to_string_lossy().to_string();

        let start = Instant::now();
        match roundtrip_single_file(&path) {
            Ok(diffs) => {
                let elapsed = start.elapsed();
                if diffs.is_empty() {
                    println!("✅ {} ({:.2?})", filename, elapsed);
                    passed += 1;
                } else {
                    println!(
                        "⚠️  {} ({:.2?}) - {} 个差异",
                        filename,
                        elapsed,
                        diffs.len()
                    );
                    for d in &diffs {
                        println!("   {}", d);
                    }
                    diff_failed += 1;
                    errors.push((filename, format!("{} diffs", diffs.len())));
                }
            }
            Err(e) => {
                let elapsed = start.elapsed();
                println!("❌ {} ({:.2?})", filename, elapsed);
                println!("   {}", e);
                write_failed += 1;
                errors.push((filename, e));
            }
        }
    }

    let total_elapsed = total_start.elapsed();

    println!("\n========================================");
    println!("📊 往返测试结果");
    println!("========================================");
    println!("✅ 完全匹配: {}", passed);
    println!("⚠️  数据差异: {}", diff_failed);
    println!("❌ 写入/读取失败: {}", write_failed);
    println!("⏱️  总耗时: {:.2?}", total_elapsed);

    if !errors.is_empty() {
        println!("\n问题文件:");
        for (f, e) in &errors {
            println!("  - {}: {}", f, e);
        }
    }

    // 写入失败是硬错误
    assert_eq!(write_failed, 0, "有 {} 个文件写入/读取失败", write_failed);
}
