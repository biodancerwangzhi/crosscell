//! 测试读取新生成的测试 Seurat 文件

use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_read_minimal_seurat_file() {
    println!("\n========================================");
    println!("Test: Read Minimal Seurat File");
    println!("========================================\n");
    
    let path = Path::new("tests/data/test_minimal_seurat.rds");
    
    if !path.exists() {
        eprintln!("⚠️  Test file not found, skipping");
        return;
    }
    
    println!("📖 Reading RDS file...");
    let obj = read_rds(path).expect("Failed to read RDS");
    
    println!("✓ Successfully read RDS file");
    println!("  Object type: {}", obj.type_name());
    println!("  Debug: {:?}", obj);
    
    println!("\n✅ Test PASSED");
}
