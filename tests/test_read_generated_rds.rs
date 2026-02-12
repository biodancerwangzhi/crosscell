//! 测试读取 data/generated 目录下的所有 RDS 文件
//!
//! 使用 crosscell 的 rds 模块解析真实的 Seurat 数据集

use crosscell::rds::{
    parse::{parse_rds_file, ParseRdsOptions},
    RdsFile, RObject,
};
use std::path::Path;
use std::time::Instant;

/// 获取 RObject 的简要描述
fn describe_robject(obj: &RObject) -> String {
    match obj {
        RObject::Null => "Null".to_string(),
        RObject::IntegerVector(v) => format!("IntegerVector[{}]", v.data.len()),
        RObject::LogicalVector(v) => format!("LogicalVector[{}]", v.data.len()),
        RObject::DoubleVector(v) => format!("DoubleVector[{}]", v.data.len()),
        RObject::RawVector(v) => format!("RawVector[{}]", v.data.len()),
        RObject::ComplexVector(v) => format!("ComplexVector[{}]", v.data.len()),
        RObject::StringVector(v) => format!("StringVector[{}]", v.data.len()),
        RObject::GenericVector(v) => format!("GenericVector[{}]", v.data.len()),
        RObject::PairList(v) => format!("PairList[{}]", v.data.len()),
        RObject::S4Object(s4) => format!("S4Object({}::{})", s4.package_name, s4.class_name),
        RObject::SymbolIndex(idx) => format!("Symbol(#{})", idx.index),
        RObject::EnvironmentIndex(idx) => format!("Environment({:?})", idx.env_type),
        RObject::BuiltInFunction(f) => format!("BuiltIn({})", f.name),
        RObject::LanguageObject(lang) => format!("Language({})", lang.function_name),
        RObject::ExpressionVector(v) => format!("ExpressionVector[{}]", v.data.len()),
        RObject::ExternalPointerIndex(idx) => format!("ExternalPointer(#{})", idx.index),
    }
}

/// 打印 RdsFile 的摘要信息
fn print_rds_summary(rds: &RdsFile, filename: &str) {
    println!("  📄 {}", filename);
    println!("     Format: v{}", rds.format_version);
    println!("     Symbols: {}", rds.symbols.len());
    println!("     Environments: {}", rds.environments.len());
    println!("     Root: {}", describe_robject(&rds.object));
}

/// 读取单个 RDS 文件
fn read_rds_file(path: &Path) -> Result<RdsFile, String> {
    let options = ParseRdsOptions::default();
    parse_rds_file(path, &options)
        .map_err(|e| format!("{}", e))
}

#[test]
fn test_read_all_generated_rds() {
    let generated_dir = Path::new("data/generated");
    if !generated_dir.exists() {
        eprintln!("⚠️  data/generated 目录不存在");
        return;
    }

    println!("\n========================================");
    println!("📂 读取 data/generated 目录下的 RDS 文件");
    println!("========================================\n");

    let mut passed = 0;
    let mut failed = 0;
    let mut failed_files: Vec<(String, String)> = Vec::new();

    let total_start = Instant::now();

    // 收集所有 RDS 文件
    let mut rds_files: Vec<_> = std::fs::read_dir(generated_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "rds"))
        .collect();
    
    // 按文件名排序
    rds_files.sort_by_key(|e| e.file_name());

    println!("找到 {} 个 RDS 文件\n", rds_files.len());

    for entry in rds_files {
        let path = entry.path();
        let filename = path.file_name().unwrap().to_string_lossy().to_string();

        let start = Instant::now();
        match read_rds_file(&path) {
            Ok(rds) => {
                let elapsed = start.elapsed();
                println!("✅ {} ({:.2?})", filename, elapsed);
                print_rds_summary(&rds, &filename);
                println!();
                passed += 1;
            }
            Err(e) => {
                let elapsed = start.elapsed();
                println!("❌ {} ({:.2?})", filename, elapsed);
                println!("   Error: {}", e);
                println!();
                failed_files.push((filename.clone(), e));
                failed += 1;
            }
        }
    }

    let total_elapsed = total_start.elapsed();

    println!("========================================");
    println!("📊 读取结果汇总");
    println!("========================================");
    println!("✅ 成功: {}", passed);
    println!("❌ 失败: {}", failed);
    println!("⏱️  总耗时: {:.2?}", total_elapsed);
    println!("========================================");

    if !failed_files.is_empty() {
        println!("\n❌ 失败文件列表:");
        for (filename, error) in &failed_files {
            println!("  - {}: {}", filename, error);
        }
    }

    // 计算通过率
    let total = passed + failed;
    let pass_rate = if total > 0 { (passed as f64 / total as f64) * 100.0 } else { 100.0 };
    println!("\n📈 通过率: {:.1}%", pass_rate);
}

/// 单独测试 Seurat V4 文件
#[test]
fn test_read_seurat_v4_pbmc3k() {
    let path = Path::new("data/generated/seurat_v4_pbmc3k_raw.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V4 pbmc3k...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v4_pbmc3k_raw.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}

/// 单独测试 Seurat V5 文件
#[test]
fn test_read_seurat_v5_pbmc3k() {
    let path = Path::new("data/generated/seurat_v5_pbmc3k_raw.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V5 pbmc3k...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v5_pbmc3k_raw.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}

/// 单独测试 Seurat V5 processed 文件
#[test]
fn test_read_seurat_v5_pbmc3k_processed() {
    let path = Path::new("data/generated/seurat_v5_pbmc3k_processed.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V5 pbmc3k processed...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v5_pbmc3k_processed.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}

/// 单独测试 Seurat V4 processed 文件
#[test]
fn test_read_seurat_v4_pbmc3k_processed() {
    let path = Path::new("data/generated/seurat_v4_pbmc3k_processed.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V4 pbmc3k processed...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v4_pbmc3k_processed.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}


/// 单独测试 Seurat V4 bmcite processed 文件（用于调试）
#[test]
fn test_read_seurat_v4_bmcite_processed() {
    let path = Path::new("data/generated/seurat_v4_bmcite_processed.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V4 bmcite processed...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v4_bmcite_processed.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}

/// 单独测试 Seurat V4 bmcite raw 文件（用于调试）
#[test]
fn test_read_seurat_v4_bmcite_raw() {
    let path = Path::new("data/generated/seurat_v4_bmcite_raw.rds");
    if !path.exists() {
        eprintln!("⚠️  文件不存在: {}", path.display());
        return;
    }

    println!("\n测试读取 Seurat V4 bmcite raw...");
    let start = Instant::now();
    
    match read_rds_file(path) {
        Ok(rds) => {
            let elapsed = start.elapsed();
            println!("✅ 读取成功 ({:.2?})", elapsed);
            print_rds_summary(&rds, "seurat_v4_bmcite_raw.rds");
        }
        Err(e) => {
            println!("❌ 读取失败: {}", e);
        }
    }
}
