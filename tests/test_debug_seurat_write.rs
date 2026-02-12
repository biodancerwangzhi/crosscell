//! Debug test for Seurat writing

use crosscell::rds::{S4Object, IntegerVector, DoubleVector, StringVector, GenericVector, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

fn make_named_list(items: Vec<(&str, RObject)>) -> RObject {
    let mut gv = GenericVector::default();
    let mut names_sv = StringVector::default();
    for (name, val) in items {
        gv.data.push(val);
        names_sv.add(name.to_string(), StringEncoding::Utf8);
    }
    gv.attributes.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);
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
fn test_write_simple_s4_with_nested_list() {
    let mut s4 = S4Object::default();
    s4.class_name = "TestClass".to_string();
    s4.class_encoding = StringEncoding::Utf8;

    s4.attributes.add("i".to_string(),
        RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3], ..Default::default() }),
        StringEncoding::Utf8);
    s4.attributes.add("nested".to_string(),
        make_named_list(vec![("key", make_string_vec(&["value"]))]),
        StringEncoding::Utf8);

    let obj = RObject::S4Object(s4);

    let output_path = "tests/data/test_debug_s4_nested.rds";
    let result = write_rds(Path::new(output_path), &obj);
    assert!(result.is_ok(), "Failed to write S4: {:?}", result.err());

    let read_result = read_rds(Path::new(output_path));
    assert!(read_result.is_ok(), "Failed to read S4: {:?}", read_result.err());
}

#[test]
fn test_write_named_list_with_s4() {
    let mut s4 = S4Object::default();
    s4.class_name = "dgCMatrix".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.attributes.add("x".to_string(),
        RObject::DoubleVector(DoubleVector { data: vec![1.0, 2.0], ..Default::default() }),
        StringEncoding::Utf8);

    let named_list = make_named_list(vec![
        ("counts", RObject::S4Object(s4)),
        ("features", make_string_vec(&["Gene1"])),
    ]);

    let output_path = "tests/data/test_debug_named_list_with_s4.rds";
    let result = write_rds(Path::new(output_path), &named_list);
    assert!(result.is_ok(), "Failed to write: {:?}", result.err());

    let read_result = read_rds(Path::new(output_path));
    assert!(read_result.is_ok(), "Failed to read: {:?}", read_result.err());
}
