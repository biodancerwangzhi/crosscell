//! SingleCellExperiment extraction error types

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SceError {
    #[error("Not an S4 object")]
    NotS4Object,

    #[error("Not a SingleCellExperiment object: {0}")]
    NotSceObject(String),

    #[error("Missing slot: {0}")]
    MissingSlot(String),

    #[error("Wrong class: expected {expected}, got {actual:?}")]
    WrongClass {
        expected: String,
        actual: Vec<String>,
    },

    #[error("Missing assay: {0}")]
    MissingAssay(String),

    #[error("Invalid slot type for '{0}': expected {1}")]
    InvalidSlotType(&'static str, &'static str),

    #[error("Invalid dimensions")]
    InvalidDimensions,

    #[error("Invalid assays structure")]
    InvalidAssaysStructure,

    #[error("Invalid colData structure")]
    InvalidColDataStructure,

    #[error("Invalid rowData structure")]
    InvalidRowDataStructure,

    #[error("Invalid reducedDims structure")]
    InvalidReducedDimsStructure,

    #[error("Not a data.frame")]
    NotDataFrame,

    #[error("Not a matrix")]
    NotMatrix,

    #[error("Unsupported column type: {0}")]
    UnsupportedColumnType(String),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Validation error: {0}")]
    ValidationError(String),

    #[error("RDS error: {0}")]
    RdsError(#[from] crate::rds::RdsError),

    #[error("IO error: {0}")]
    IoError(String),

    #[error("Conversion error: {0}")]
    ConversionError(String),
}

pub type Result<T> = std::result::Result<T, SceError>;
