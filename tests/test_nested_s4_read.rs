use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_nested_with_s4() {
    let path = Path::new("tests/data/nested_with_s4.rds");

    match read_rds(path) {
        Ok(robj) => {
            println!("Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            match &robj {
                RObject::GenericVector(gv) => {
                    println!("GenericVector with {} elements", gv.data.len());
                    if let Some(names) = gv.attributes.get_names() {
                        for (name, value) in names.iter().zip(gv.data.iter()) {
                            println!("  - {}: {}", name, value.type_name());
                        }
                    }
                }
                _ => {
                    println!("Type: {}", robj.type_name());
                }
            }
        }
        Err(e) => {
            eprintln!("Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
