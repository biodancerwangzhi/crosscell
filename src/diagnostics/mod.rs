//! Data diagnostics and auto-fix module
//!
//! This module provides intelligent data problem detection and automatic fixing
//! capabilities, inspired by anndata_Seurat_Utilities.
//!
//! # Features
//! - Detect special characters in column names
//! - Detect empty columns
//! - Detect risky data types
//! - Auto-fix detected issues
//! - Generate audit logs
//! - User-friendly error messages

mod audit;
mod cleaner;
mod detector;
mod error_messages;

pub use audit::{AuditLog, Modification};
pub use cleaner::{CleaningOptions, CleaningReport, DataCleaner};
pub use detector::{detect_issues, DataIssue, DiagnosticReport, IssueSeverity, IssueType};
pub use error_messages::{format_error_report, format_error_with_suggestions, ErrorContext};
