//! Test reading simple dgCMatrix

use crosscell::read_rds;
use std::path::Path;

#[test]
fn test_read_simple_dgc() {
    let path = Path::new("tests/data/test_dgc_simple.rds");
    
    if !path.exists() {
        eprintln!("Skipping: {} not found", path.display());
        return;
    }
    
    println!("\n========================================");
    println!("Reading simple dgCMatrix");
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
