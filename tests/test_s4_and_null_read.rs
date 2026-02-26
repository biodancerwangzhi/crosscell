use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_s4_and_null() {
    let path = Path::new("tests/data/s4_and_null.rds");

    println!("\n========================================");
    println!("Reading list with S4 and NULL");
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
