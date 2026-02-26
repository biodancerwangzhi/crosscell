use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_empty_dataframe() {
    let path = Path::new("tests/data/empty_dataframe.rds");

    println!("\n========================================");
    println!("Reading empty data.frame");
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
