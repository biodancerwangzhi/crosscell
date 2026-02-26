//! Data issue detection module
//!
//! Detects potential problems in single-cell data that may cause conversion issues.

use crate::ir::{DataFrame, SingleCellData};
use std::collections::HashSet;

/// Severity level of detected issues
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IssueSeverity {
    /// Critical issues that will cause conversion failure
    Critical,
    /// Warning issues that may cause data loss or unexpected behavior
    Warning,
    /// Informational issues that are auto-fixable
    Info,
}

impl std::fmt::Display for IssueSeverity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IssueSeverity::Critical => write!(f, "CRITICAL"),
            IssueSeverity::Warning => write!(f, "WARNING"),
            IssueSeverity::Info => write!(f, "INFO"),
        }
    }
}

/// Type of detected issue
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IssueType {
    /// Column name contains special characters
    SpecialCharInName { column: String, chars: Vec<char> },
    /// Duplicate column names
    DuplicateColumnName { column: String, count: usize },
    /// Empty column (all NA/null values)
    EmptyColumn { column: String },
    /// Risky data type that may cause issues
    RiskyDataType { column: String, dtype: String },
    /// Very long column name
    LongColumnName { column: String, length: usize },
    /// Column name starts with number
    NumericStartName { column: String },
}

impl std::fmt::Display for IssueType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IssueType::SpecialCharInName { column, chars } => {
                write!(
                    f,
                    "Column '{}' contains special characters: {:?}",
                    column, chars
                )
            }
            IssueType::DuplicateColumnName { column, count } => {
                write!(f, "Column '{}' appears {} times", column, count)
            }
            IssueType::EmptyColumn { column } => {
                write!(f, "Column '{}' is empty (all null values)", column)
            }
            IssueType::RiskyDataType { column, dtype } => {
                write!(f, "Column '{}' uses risky dtype: {}", column, dtype)
            }
            IssueType::LongColumnName { column, length } => {
                write!(f, "Column '{}' name is too long ({} chars)", column, length)
            }
            IssueType::NumericStartName { column } => {
                write!(f, "Column '{}' starts with a number", column)
            }
        }
    }
}

/// A detected data issue
#[derive(Debug, Clone)]
pub struct DataIssue {
    pub severity: IssueSeverity,
    pub issue_type: IssueType,
    pub auto_fixable: bool,
    pub suggestion: String,
}

impl DataIssue {
    pub fn new(
        severity: IssueSeverity,
        issue_type: IssueType,
        auto_fixable: bool,
        suggestion: &str,
    ) -> Self {
        Self {
            severity,
            issue_type,
            auto_fixable,
            suggestion: suggestion.to_string(),
        }
    }
}

impl std::fmt::Display for DataIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let fix_status = if self.auto_fixable {
            "auto-fixable"
        } else {
            "manual fix required"
        };
        write!(
            f,
            "[{}] {} ({})",
            self.severity, self.issue_type, fix_status
        )
    }
}

/// Diagnostic report containing all detected issues
#[derive(Debug, Clone)]
pub struct DiagnosticReport {
    pub issues: Vec<DataIssue>,
    pub compatibility_score: u8,
    pub total_issues: usize,
    pub critical_count: usize,
    pub warning_count: usize,
    pub info_count: usize,
    pub auto_fixable_count: usize,
}

impl DiagnosticReport {
    pub fn new() -> Self {
        Self {
            issues: Vec::new(),
            compatibility_score: 100,
            total_issues: 0,
            critical_count: 0,
            warning_count: 0,
            info_count: 0,
            auto_fixable_count: 0,
        }
    }

    pub fn add_issue(&mut self, issue: DataIssue) {
        match issue.severity {
            IssueSeverity::Critical => {
                self.critical_count += 1;
                self.compatibility_score = self.compatibility_score.saturating_sub(20);
            }
            IssueSeverity::Warning => {
                self.warning_count += 1;
                self.compatibility_score = self.compatibility_score.saturating_sub(5);
            }
            IssueSeverity::Info => {
                self.info_count += 1;
                self.compatibility_score = self.compatibility_score.saturating_sub(1);
            }
        }
        if issue.auto_fixable {
            self.auto_fixable_count += 1;
        }
        self.total_issues += 1;
        self.issues.push(issue);
    }

    pub fn is_convertible(&self) -> bool {
        self.critical_count == 0
    }

    pub fn has_issues(&self) -> bool {
        self.total_issues > 0
    }
}

impl Default for DiagnosticReport {
    fn default() -> Self {
        Self::new()
    }
}

/// Special characters that may cause issues in column names
const SPECIAL_CHARS: &[char] = &[
    '/', ' ', ',', ';', ':', '\\', '|', '?', '*', '<', '>', '"', '\'', '`', '(', ')', '[', ']',
    '{', '}', '@', '#', '$', '%', '^', '&', '=', '+', '~',
];

/// Detect issues in single-cell data
pub fn detect_issues(data: &SingleCellData) -> DiagnosticReport {
    let mut report = DiagnosticReport::new();

    // Check cell metadata
    detect_dataframe_issues(&data.cell_metadata, "obs", &mut report);

    // Check gene metadata
    detect_dataframe_issues(&data.gene_metadata, "var", &mut report);

    report
}

/// Detect issues in a DataFrame
fn detect_dataframe_issues(df: &DataFrame, prefix: &str, report: &mut DiagnosticReport) {
    let mut seen_names: HashSet<String> = HashSet::new();
    let mut name_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();

    // Count column name occurrences
    for col_name in &df.columns {
        *name_counts.entry(col_name.clone()).or_insert(0) += 1;
    }

    for col_name in &df.columns {
        let full_name = format!("{}.{}", prefix, col_name);

        // Check for special characters
        let special_chars: Vec<char> = col_name
            .chars()
            .filter(|c| SPECIAL_CHARS.contains(c))
            .collect();

        if !special_chars.is_empty() {
            let sanitized = sanitize_column_name(col_name);
            report.add_issue(DataIssue::new(
                IssueSeverity::Info,
                IssueType::SpecialCharInName {
                    column: full_name.clone(),
                    chars: special_chars,
                },
                true,
                &format!("Will sanitize to '{}'", sanitized),
            ));
        }

        // Check for duplicate names
        if let Some(&count) = name_counts.get(col_name) {
            if count > 1 && !seen_names.contains(col_name) {
                report.add_issue(DataIssue::new(
                    IssueSeverity::Warning,
                    IssueType::DuplicateColumnName {
                        column: full_name.clone(),
                        count,
                    },
                    true,
                    "Will add _dup2, _dup3, ... suffixes",
                ));
                seen_names.insert(col_name.clone());
            }
        }

        // Check for column names starting with numbers
        if col_name
            .chars()
            .next()
            .map(|c| c.is_ascii_digit())
            .unwrap_or(false)
        {
            report.add_issue(DataIssue::new(
                IssueSeverity::Info,
                IssueType::NumericStartName {
                    column: full_name.clone(),
                },
                true,
                "Will prefix with 'X'",
            ));
        }

        // Check for very long column names (R has 10000 char limit, but shorter is better)
        if col_name.len() > 256 {
            report.add_issue(DataIssue::new(
                IssueSeverity::Warning,
                IssueType::LongColumnName {
                    column: full_name.clone(),
                    length: col_name.len(),
                },
                true,
                "Will truncate to 256 characters",
            ));
        }
    }
}

/// Sanitize a column name by replacing special characters with underscores
pub fn sanitize_column_name(name: &str) -> String {
    let mut result = String::with_capacity(name.len());
    let mut last_was_underscore = false;

    for c in name.chars() {
        if SPECIAL_CHARS.contains(&c) {
            if !last_was_underscore {
                result.push('_');
                last_was_underscore = true;
            }
        } else {
            result.push(c);
            last_was_underscore = false;
        }
    }

    // Remove trailing underscore
    if result.ends_with('_') {
        result.pop();
    }

    // Handle names starting with numbers
    if result
        .chars()
        .next()
        .map(|c| c.is_ascii_digit())
        .unwrap_or(false)
    {
        result = format!("X{}", result);
    }

    // Handle empty result
    if result.is_empty() {
        result = "unnamed".to_string();
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sanitize_column_name() {
        assert_eq!(sanitize_column_name("Kir3dl1/2"), "Kir3dl1_2");
        assert_eq!(sanitize_column_name("cell type"), "cell_type");
        assert_eq!(sanitize_column_name("a/b/c"), "a_b_c");
        assert_eq!(sanitize_column_name("normal_name"), "normal_name");
        assert_eq!(sanitize_column_name("123abc"), "X123abc");
        assert_eq!(sanitize_column_name("a  b"), "a_b");
    }

    #[test]
    fn test_issue_severity_display() {
        assert_eq!(format!("{}", IssueSeverity::Critical), "CRITICAL");
        assert_eq!(format!("{}", IssueSeverity::Warning), "WARNING");
        assert_eq!(format!("{}", IssueSeverity::Info), "INFO");
    }

    #[test]
    fn test_diagnostic_report_scoring() {
        let mut report = DiagnosticReport::new();
        assert_eq!(report.compatibility_score, 100);

        report.add_issue(DataIssue::new(
            IssueSeverity::Info,
            IssueType::SpecialCharInName {
                column: "test".to_string(),
                chars: vec!['/'],
            },
            true,
            "test",
        ));
        assert_eq!(report.compatibility_score, 99);

        report.add_issue(DataIssue::new(
            IssueSeverity::Warning,
            IssueType::DuplicateColumnName {
                column: "test".to_string(),
                count: 2,
            },
            true,
            "test",
        ));
        assert_eq!(report.compatibility_score, 94);

        report.add_issue(DataIssue::new(
            IssueSeverity::Critical,
            IssueType::EmptyColumn {
                column: "test".to_string(),
            },
            false,
            "test",
        ));
        assert_eq!(report.compatibility_score, 74);
    }
}
