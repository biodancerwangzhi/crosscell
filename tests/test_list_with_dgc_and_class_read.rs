use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_list_with_dgc_and_class() {
    let path = Path::new("tests/data/list_with_dgc_and_class.rds");
    
    println!("\n========================================");
    println!("Reading list with dgCMatrix and class");
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
