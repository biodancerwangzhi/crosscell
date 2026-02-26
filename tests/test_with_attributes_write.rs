//! Test writing objects with attributes

use crosscell::rds::{IntegerVector, StringEncoding, StringVector};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

#[test]
fn test_simple_with_attributes() {
    // Create an integer vector with a class attribute
    let mut iv = IntegerVector {
        data: vec![1, 2, 3],
        ..Default::default()
    };
    let mut class_sv = StringVector::default();
    class_sv.add("myclass".to_string(), StringEncoding::Utf8);
    iv.attributes.add(
        "class".to_string(),
        RObject::StringVector(class_sv),
        StringEncoding::Utf8,
    );

    let obj = RObject::IntegerVector(iv);

    let path = "tests/data/test_with_attributes.rds";
    write_rds(Path::new(path), &obj).expect("Failed to write");

    let read_obj = read_rds(Path::new(path)).expect("Failed to read");

    match read_obj {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data.len(), 3);
            assert_eq!(iv.data[0], 1);
            // Check class attribute
            let class = iv.attributes.get_class();
            assert!(class.is_some(), "Should have class attribute");
            assert_eq!(class.unwrap()[0], "myclass");
        }
        _ => panic!("Expected IntegerVector, got {:?}", read_obj.type_name()),
    }
}
