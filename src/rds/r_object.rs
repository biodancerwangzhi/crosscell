//! R 对象数据结构 - 对应 rds2cpp 的 RObject.hpp

use num_complex::Complex64;
use super::sexp_type::SEXPType;
use super::string_encoding::StringEncoding;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SymbolIndex { pub index: usize }

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnvironmentIndex { pub index: usize, pub env_type: SEXPType }
impl Default for EnvironmentIndex {
    fn default() -> Self { Self { index: usize::MAX, env_type: SEXPType::GlobalEnv } }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExternalPointerIndex { pub index: usize }

#[derive(Debug, Clone, Default)]
pub struct Attributes {
    pub names: Vec<String>,
    pub encodings: Vec<StringEncoding>,
    pub values: Vec<RObject>,
}
impl PartialEq for Attributes {
    fn eq(&self, other: &Self) -> bool {
        self.names == other.names && self.encodings == other.encodings && self.values == other.values
    }
}
impl Eq for Attributes {}
impl Attributes {
    pub fn new() -> Self { Self::default() }
    pub fn add(&mut self, name: String, value: RObject, encoding: StringEncoding) {
        self.names.push(name); self.values.push(value); self.encodings.push(encoding);
    }
    pub fn is_empty(&self) -> bool { self.names.is_empty() }
    pub fn len(&self) -> usize { self.names.len() }
    pub fn get(&self, name: &str) -> Option<&RObject> {
        self.names.iter().position(|n| n == name).map(|i| &self.values[i])
    }
    /// 获取 "names" 属性（如果存在且是字符串向量）
    pub fn get_names(&self) -> Option<&[String]> {
        self.get("names").and_then(|obj| {
            if let RObject::StringVector(sv) = obj {
                Some(sv.data.as_slice())
            } else {
                None
            }
        })
    }
    /// 获取 "class" 属性
    pub fn get_class(&self) -> Option<&[String]> {
        self.get("class").and_then(|obj| {
            if let RObject::StringVector(sv) = obj {
                Some(sv.data.as_slice())
            } else {
                None
            }
        })
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntegerVector { pub data: Vec<i32>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LogicalVector { pub data: Vec<i32>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DoubleVector { pub data: Vec<f64>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RawVector { pub data: Vec<u8>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ComplexVector { pub data: Vec<Complex64>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StringVector {
    pub data: Vec<String>, pub encodings: Vec<StringEncoding>,
    pub missing: Vec<bool>, pub attributes: Attributes,
}
impl StringVector {
    pub fn add(&mut self, value: String, encoding: StringEncoding) {
        self.data.push(value); self.encodings.push(encoding); self.missing.push(false);
    }
    pub fn add_missing(&mut self) {
        self.data.push(String::new()); self.encodings.push(StringEncoding::None); self.missing.push(true);
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct GenericVector { pub data: Vec<RObject>, pub attributes: Attributes }
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PairList {
    pub data: Vec<RObject>, pub has_tag: Vec<bool>,
    pub tag_names: Vec<String>, pub tag_encodings: Vec<StringEncoding>,
    pub attributes: Attributes,
}
impl PairList {
    pub fn add_tagged(&mut self, tag: String, value: RObject, encoding: StringEncoding) {
        self.data.push(value); self.has_tag.push(true);
        self.tag_names.push(tag); self.tag_encodings.push(encoding);
    }
    pub fn add_untagged(&mut self, value: RObject) {
        self.data.push(value); self.has_tag.push(false);
        self.tag_names.push(String::new()); self.tag_encodings.push(StringEncoding::None);
    }
}
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct S4Object {
    pub class_name: String, pub class_encoding: StringEncoding,
    pub package_name: String, pub package_encoding: StringEncoding,
    pub attributes: Attributes,
}
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BuiltInFunction { pub name: String }
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LanguageObject {
    pub function_name: String, pub function_encoding: StringEncoding,
    pub argument_values: Vec<RObject>, pub argument_names: Vec<String>,
    pub argument_has_name: Vec<bool>, pub argument_encodings: Vec<StringEncoding>,
    pub attributes: Attributes,
}
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ExpressionVector { pub data: Vec<RObject>, pub attributes: Attributes }


#[derive(Debug, Clone, PartialEq)]
pub enum RObject {
    Null, SymbolIndex(SymbolIndex), EnvironmentIndex(EnvironmentIndex),
    ExternalPointerIndex(ExternalPointerIndex), IntegerVector(IntegerVector),
    LogicalVector(LogicalVector), DoubleVector(DoubleVector), RawVector(RawVector),
    ComplexVector(ComplexVector), StringVector(StringVector), GenericVector(GenericVector),
    PairList(PairList), S4Object(S4Object), BuiltInFunction(BuiltInFunction),
    LanguageObject(LanguageObject), ExpressionVector(ExpressionVector),
}
impl Default for RObject { fn default() -> Self { RObject::Null } }
impl RObject {
    pub fn sexp_type(&self) -> SEXPType {
        match self {
            RObject::Null => SEXPType::Nil, RObject::SymbolIndex(_) => SEXPType::Sym,
            RObject::EnvironmentIndex(env) => env.env_type,
            RObject::ExternalPointerIndex(_) => SEXPType::Extptr,
            RObject::IntegerVector(_) => SEXPType::Int, RObject::LogicalVector(_) => SEXPType::Lgl,
            RObject::DoubleVector(_) => SEXPType::Real, RObject::RawVector(_) => SEXPType::Raw,
            RObject::ComplexVector(_) => SEXPType::Cplx, RObject::StringVector(_) => SEXPType::Str,
            RObject::GenericVector(_) => SEXPType::Vec, RObject::PairList(_) => SEXPType::List,
            RObject::S4Object(_) => SEXPType::S4, RObject::BuiltInFunction(_) => SEXPType::Builtin,
            RObject::LanguageObject(_) => SEXPType::Lang, RObject::ExpressionVector(_) => SEXPType::Expr,
        }
    }
    pub fn is_null(&self) -> bool { matches!(self, RObject::Null) }
    pub fn is_atomic_vector(&self) -> bool {
        matches!(self, RObject::IntegerVector(_) | RObject::LogicalVector(_) |
            RObject::DoubleVector(_) | RObject::RawVector(_) |
            RObject::ComplexVector(_) | RObject::StringVector(_))
    }
    pub fn attributes(&self) -> Option<&Attributes> {
        match self {
            RObject::IntegerVector(v) => Some(&v.attributes),
            RObject::LogicalVector(v) => Some(&v.attributes),
            RObject::DoubleVector(v) => Some(&v.attributes),
            RObject::RawVector(v) => Some(&v.attributes),
            RObject::ComplexVector(v) => Some(&v.attributes),
            RObject::StringVector(v) => Some(&v.attributes),
            RObject::GenericVector(v) => Some(&v.attributes),
            RObject::PairList(v) => Some(&v.attributes),
            RObject::S4Object(v) => Some(&v.attributes),
            RObject::LanguageObject(v) => Some(&v.attributes),
            RObject::ExpressionVector(v) => Some(&v.attributes),
            _ => None,
        }
    }
    pub fn attributes_mut(&mut self) -> Option<&mut Attributes> {
        match self {
            RObject::IntegerVector(v) => Some(&mut v.attributes),
            RObject::LogicalVector(v) => Some(&mut v.attributes),
            RObject::DoubleVector(v) => Some(&mut v.attributes),
            RObject::RawVector(v) => Some(&mut v.attributes),
            RObject::ComplexVector(v) => Some(&mut v.attributes),
            RObject::StringVector(v) => Some(&mut v.attributes),
            RObject::GenericVector(v) => Some(&mut v.attributes),
            RObject::PairList(v) => Some(&mut v.attributes),
            RObject::S4Object(v) => Some(&mut v.attributes),
            RObject::LanguageObject(v) => Some(&mut v.attributes),
            RObject::ExpressionVector(v) => Some(&mut v.attributes),
            _ => None,
        }
    }
    pub fn type_name(&self) -> &'static str {
        match self {
            RObject::Null => "NULL",
            RObject::SymbolIndex(_) => "symbol",
            RObject::EnvironmentIndex(_) => "environment",
            RObject::ExternalPointerIndex(_) => "externalptr",
            RObject::IntegerVector(_) => "integer",
            RObject::LogicalVector(_) => "logical",
            RObject::DoubleVector(_) => "double",
            RObject::RawVector(_) => "raw",
            RObject::ComplexVector(_) => "complex",
            RObject::StringVector(_) => "character",
            RObject::GenericVector(_) => "list",
            RObject::PairList(_) => "pairlist",
            RObject::S4Object(_) => "S4",
            RObject::BuiltInFunction(_) => "builtin",
            RObject::LanguageObject(_) => "language",
            RObject::ExpressionVector(_) => "expression",
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn test_symbol_index_default() { assert_eq!(SymbolIndex::default().index, 0); }
    #[test] fn test_environment_index_default() {
        let idx = EnvironmentIndex::default();
        assert_eq!(idx.index, usize::MAX);
        assert_eq!(idx.env_type, SEXPType::GlobalEnv);
    }
    #[test] fn test_attributes_operations() {
        let mut attrs = Attributes::new();
        assert!(attrs.is_empty());
        attrs.add("names".to_string(), RObject::Null, StringEncoding::Utf8);
        assert!(!attrs.is_empty());
        assert!(attrs.get("names").is_some());
    }
    #[test] fn test_robject_sexp_type() {
        assert_eq!(RObject::Null.sexp_type(), SEXPType::Nil);
        assert_eq!(RObject::IntegerVector(IntegerVector::default()).sexp_type(), SEXPType::Int);
    }
}
