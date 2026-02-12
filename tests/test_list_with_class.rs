use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_list_with_class() {
    let path = Path::new("tests/data/list_with_class.rds");
    
    println!("\n========================================");
    println!("Reading list with class attribute");
    println!("========================================\n");
    
    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            println!("\nObject structure:");
            println!("{:#?}", robj);
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
