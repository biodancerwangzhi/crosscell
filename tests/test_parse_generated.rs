//! 测试解析 data/generated 目录下的所有 RDS 文件

use crosscell::rds::parse::{parse_rds_file, ParseRdsOptions};
use std::path::Path;

/// 测试解析单个文件
fn test_parse_file(path: &Path) -> Result<String, String> {
    let options = ParseRdsOptions::default();
    
    match parse_rds_file(path, &options) {
        Ok(rds) => {
            let obj_type = format!("{:?}", rds.object.sexp_type());
            Ok(format!("✅ {} - {}", path.file_name().unwrap().to_string_lossy(), obj_type))
        }
        Err(e) => {
            Err(format!("❌ {} - {}", path.file_name().unwrap().to_string_lossy(), e))
        }
    }
}

#[test]
fn test_parse_all_generated_rds() {
    let generated_dir = Path::new("data/generated");
    if !generated_dir.exists() {
        eprintln!("⚠️  Generated data directory not found");
        return;
    }

    let mut passed = 0;
    let mut failed = 0;
    let mut results: Vec<String> = Vec::new();

    for entry in std::fs::read_dir(generated_dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.extension().map_or(false, |ext| ext == "rds") {
            match test_parse_file(&path) {
                Ok(msg) => {
                    results.push(msg);
                    passed += 1;
                }
                Err(msg) => {
                    results.push(msg);
                    failed += 1;
                }
            }
        }
    }

    // 排序结果
    results.sort();
    
    println!("\n========================================");
    println!("Parse Test Results for data/generated/");
    println!("========================================\n");
    
    for result in &results {
        println!("{}", result);
    }
    
    println!("\n========================================");
    println!("Summary: {} passed, {} failed", passed, failed);
    println!("========================================");
}
