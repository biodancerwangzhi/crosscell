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

mod detector;
mod cleaner;
mod audit;
mod error_messages;

pub use detector::{DataIssue, IssueSeverity, IssueType, DiagnosticReport, detect_issues};
pub use cleaner::{DataCleaner, CleaningReport, CleaningOptions};
pub use audit::{AuditLog, Modification};
pub use error_messages::{format_error_with_suggestions, format_error_report, ErrorContext};
