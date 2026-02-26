use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_nested_null() {
    let path = Path::new("tests/data/nested_null.rds");

    println!("\n========================================");
    println!("Reading nested list with NULL");
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
