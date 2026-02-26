use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_read_seurat_simplified() {
    let path = Path::new("tests/data/seurat_minimal_simplified.rds");

    println!("\n========================================");
    println!("Reading Seurat simplified RDS");
    println!("========================================\n");

    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            println!("\nFirst 100 chars of debug output:");
            let debug_str = format!("{:#?}", robj);
            println!("{}", &debug_str[..debug_str.len().min(500)]);
            if debug_str.len() > 500 {
                println!("... (truncated)");
            }
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
