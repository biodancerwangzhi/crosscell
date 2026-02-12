use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_nested_list() {
    let path = Path::new("tests/data/nested_list.rds");

    match read_rds(path) {
        Ok(robj) => {
            println!("Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            match &robj {
                RObject::GenericVector(gv) => {
                    println!("GenericVector with {} elements", gv.data.len());
                    if let Some(names) = gv.attributes.get_names() {
                        for name in names {
                            println!("  - {}", name);
                        }
                    }
                }
                _ => {
                    println!("Not a GenericVector, it's {}", robj.type_name());
                }
            }
        }
        Err(e) => {
            eprintln!("Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
