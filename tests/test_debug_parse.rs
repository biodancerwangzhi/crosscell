//! Debug test for RDS parsing

use crosscell::read_rds;
use std::path::Path;

#[test]
fn test_debug_seurat_simple() {
    // Enable debug output
    std::env::set_var("RUST_BACKTRACE", "1");
    
    let path = Path::new("tests/data/seurat_simple.rds");
    
    println!("\n========================================");
    println!("Debug: Reading seurat_simple.rds");
    println!("========================================\n");
    
    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}

#[test]
fn test_debug_minimal_seurat() {
    let path = Path::new("tests/data/test_minimal_seurat.rds");
    
    if !path.exists() {
        eprintln!("Skipping: {} not found", path.display());
        return;
    }
    
    println!("\n========================================");
    println!("Debug: Reading test_minimal_seurat.rds");
    println!("========================================\n");
    
    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
