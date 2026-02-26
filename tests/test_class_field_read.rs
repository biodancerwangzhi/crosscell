use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_class_field() {
    let path = Path::new("tests/data/class_field.rds");

    println!("\n========================================");
    println!("Reading list with 'class' field name");
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
