//! Test writing named list (GenericVector with names attribute) objects

use crosscell::rds::{GenericVector, IntegerVector, StringEncoding, StringVector};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

fn make_named_list(items: Vec<(&str, RObject)>) -> RObject {
    let mut gv = GenericVector::default();
    let mut names_sv = StringVector::default();
    for (name, val) in items {
        gv.data.push(val);
        names_sv.add(name.to_string(), StringEncoding::Utf8);
    }
    gv.attributes.add(
        "names".to_string(),
        RObject::StringVector(names_sv),
        StringEncoding::Utf8,
    );
    RObject::GenericVector(gv)
}

fn make_string_vec(strings: &[&str]) -> RObject {
    let mut sv = StringVector::default();
    for s in strings {
        sv.add(s.to_string(), StringEncoding::Utf8);
    }
    RObject::StringVector(sv)
}

#[test]
fn test_write_simple_named_list() {
    let named_list = make_named_list(vec![
        (
            "a",
            RObject::IntegerVector(IntegerVector {
                data: vec![1, 2, 3],
                ..Default::default()
            }),
        ),
        ("b", make_string_vec(&["hello"])),
    ]);

    let path = "tests/data/test_named_list_write.rds";
    let result = write_rds(Path::new(path), &named_list);
    assert!(
        result.is_ok(),
        "Failed to write named list: {:?}",
        result.err()
    );

    let read_result = read_rds(Path::new(path));
    match read_result {
        Ok(obj) => {
            println!("Read back object: {}", obj.type_name());
            match obj {
                RObject::GenericVector(gv) => {
                    println!("  - {} items", gv.data.len());
                }
                _ => {
                    println!("Got type: {}", obj.type_name());
                }
            }
        }
        Err(e) => panic!("Failed to read back: {:?}", e),
    }
}

#[test]
fn test_write_nested_named_list() {
    let inner = make_named_list(vec![
        (
            "x",
            RObject::IntegerVector(IntegerVector {
                data: vec![1, 2],
                ..Default::default()
            }),
        ),
        (
            "y",
            RObject::IntegerVector(IntegerVector {
                data: vec![3, 4],
                ..Default::default()
            }),
        ),
    ]);
    let outer = make_named_list(vec![("data", inner), ("name", make_string_vec(&["test"]))]);

    let path = "tests/data/test_nested_named_list_write.rds";
    let result = write_rds(Path::new(path), &outer);
    assert!(
        result.is_ok(),
        "Failed to write nested named list: {:?}",
        result.err()
    );

    let read_result = read_rds(Path::new(path));
    assert!(
        read_result.is_ok(),
        "Failed to read back: {:?}",
        read_result.err()
    );
}

#[test]
fn test_write_named_list_with_class_attribute() {
    let mut gv = GenericVector::default();
    let mut names_sv = StringVector::default();

    gv.data.push(RObject::IntegerVector(IntegerVector {
        data: vec![1, 2, 3],
        ..Default::default()
    }));
    names_sv.add("a".to_string(), StringEncoding::Utf8);

    gv.data.push(make_string_vec(&["hello"]));
    names_sv.add("b".to_string(), StringEncoding::Utf8);

    gv.attributes.add(
        "names".to_string(),
        RObject::StringVector(names_sv),
        StringEncoding::Utf8,
    );

    // Add class attribute
    let mut class_sv = StringVector::default();
    class_sv.add("MyClass".to_string(), StringEncoding::Utf8);
    class_sv.add("list".to_string(), StringEncoding::Utf8);
    gv.attributes.add(
        "class".to_string(),
        RObject::StringVector(class_sv),
        StringEncoding::Utf8,
    );

    let obj = RObject::GenericVector(gv);

    let path = "tests/data/test_named_list_with_class.rds";
    let result = write_rds(Path::new(path), &obj);
    assert!(result.is_ok(), "Failed to write: {:?}", result.err());

    let read_result = read_rds(Path::new(path));
    assert!(
        read_result.is_ok(),
        "Failed to read back: {:?}",
        read_result.err()
    );
}

#[test]
fn test_write_simple_integer() {
    let int_vec = RObject::IntegerVector(IntegerVector {
        data: vec![1, 2, 3, 4, 5],
        ..Default::default()
    });

    let path = "tests/data/test_simple_int.rds";
    let result = write_rds(Path::new(path), &int_vec);
    assert!(
        result.is_ok(),
        "Failed to write integer vector: {:?}",
        result.err()
    );

    let read_result = read_rds(Path::new(path));
    assert!(
        read_result.is_ok(),
        "Failed to read back: {:?}",
        read_result.err()
    );
}
