//! 测试读取简单命名列表

use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_read_simple_named_list() {
    println!("\n========================================");
    println!("Test: Read Simple Named List");
    println!("========================================\n");
    
    let path = Path::new("tests/data/r_simple_named_list.rds");
    
    if !path.exists() {
        eprintln!("⚠️  Test file not found, skipping");
        return;
    }
    
    println!("📖 Reading RDS file...");
    eprintln!("DEBUG: About to call read_rds");
    
    match read_rds(path) {
        Ok(obj) => {
            println!("✓ Successfully read RDS file");
            println!("  Object type: {}", obj.type_name());
            println!("  Debug: {:?}", obj);
            println!("\n✅ Test PASSED");
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("Test failed: {:?}", e);
        }
    }
}
