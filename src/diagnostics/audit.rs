//! Audit log module
//!
//! Records all modifications made during data cleaning and conversion.

use chrono::Utc;
use serde::{Deserialize, Serialize};

/// A single modification record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Modification {
    /// Column was renamed
    ColumnRenamed { from: String, to: String },
    /// Column was dropped
    ColumnDropped { name: String, reason: String },
    /// Data type was converted
    TypeConverted { column: String, from: String, to: String },
    /// Empty column was removed
    EmptyColumnRemoved { name: String },
}

impl std::fmt::Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Modification::ColumnRenamed { from, to } => {
                write!(f, "Renamed: '{}' → '{}'", from, to)
            }
            Modification::ColumnDropped { name, reason } => {
                write!(f, "Dropped: '{}' ({})", name, reason)
            }
            Modification::TypeConverted { column, from, to } => {
                write!(f, "Type converted: '{}' ({} → {})", column, from, to)
            }
            Modification::EmptyColumnRemoved { name } => {
                write!(f, "Removed empty column: '{}'", name)
            }
        }
    }
}

/// Audit log containing all modifications
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuditLog {
    /// Timestamp of the conversion
    pub timestamp: String,
    /// Tool version
    pub tool_version: String,
    /// List of modifications
    pub modifications: Vec<Modification>,
}

impl AuditLog {
    pub fn new() -> Self {
        Self {
            timestamp: Utc::now().to_rfc3339(),
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            modifications: Vec::new(),
        }
    }

    pub fn add_modification(&mut self, modification: Modification) {
        self.modifications.push(modification);
    }

    pub fn is_empty(&self) -> bool {
        self.modifications.is_empty()
    }

    pub fn len(&self) -> usize {
        self.modifications.len()
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Deserialize from JSON string
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Generate a human-readable summary
    pub fn summary(&self) -> String {
        if self.modifications.is_empty() {
            return "No modifications were made.".to_string();
        }

        let mut summary = format!(
            "CrossCell Audit Log (v{})\nTimestamp: {}\n\nModifications ({}):\n",
            self.tool_version,
            self.timestamp,
            self.modifications.len()
        );

        for (i, modification) in self.modifications.iter().enumerate() {
            summary.push_str(&format!("  {}. {}\n", i + 1, modification));
        }

        summary
    }
}

impl Default for AuditLog {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_audit_log_creation() {
        let log = AuditLog::new();
        assert!(log.is_empty());
        assert_eq!(log.len(), 0);
    }

    #[test]
    fn test_add_modification() {
        let mut log = AuditLog::new();
        log.add_modification(Modification::ColumnRenamed {
            from: "cell/type".to_string(),
            to: "cell_type".to_string(),
        });
        
        assert!(!log.is_empty());
        assert_eq!(log.len(), 1);
    }

    #[test]
    fn test_json_serialization() {
        let mut log = AuditLog::new();
        log.add_modification(Modification::ColumnRenamed {
            from: "cell/type".to_string(),
            to: "cell_type".to_string(),
        });
        log.add_modification(Modification::ColumnDropped {
            name: "empty_col".to_string(),
            reason: "All null values".to_string(),
        });
        
        let json = log.to_json().unwrap();
        assert!(json.contains("ColumnRenamed"));
        assert!(json.contains("cell/type"));
        
        let parsed = AuditLog::from_json(&json).unwrap();
        assert_eq!(parsed.len(), 2);
    }

    #[test]
    fn test_summary() {
        let mut log = AuditLog::new();
        log.add_modification(Modification::ColumnRenamed {
            from: "cell/type".to_string(),
            to: "cell_type".to_string(),
        });
        
        let summary = log.summary();
        assert!(summary.contains("Renamed"));
        assert!(summary.contains("cell/type"));
        assert!(summary.contains("cell_type"));
    }

    #[test]
    fn test_modification_display() {
        let m1 = Modification::ColumnRenamed {
            from: "a".to_string(),
            to: "b".to_string(),
        };
        assert_eq!(format!("{}", m1), "Renamed: 'a' → 'b'");
        
        let m2 = Modification::ColumnDropped {
            name: "col".to_string(),
            reason: "empty".to_string(),
        };
        assert_eq!(format!("{}", m2), "Dropped: 'col' (empty)");
    }
}
