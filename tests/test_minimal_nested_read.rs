use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_minimal_nested() {
    let path = Path::new("tests/data/minimal_nested.rds");
    
    println!("\n========================================");
    println!("Reading minimal nested list");
    println!("========================================\n");
    
    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            println!("{:#?}", robj);
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
