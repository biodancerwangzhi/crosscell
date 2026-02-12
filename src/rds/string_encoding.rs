//! String encoding types
//!
//! Defines R string encoding types (None, Latin1, UTF8, ASCII).
//! Corresponds to rds2cpp StringEncoding.hpp

/// String encoding type
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum StringEncoding {
    /// Unspecified encoding
    None = 0,
    /// ISO-8859-1 encoding
    Latin1 = 1,
    /// UTF-8 encoding (default)
    #[default]
    Utf8 = 2,
    /// ASCII encoding
    Ascii = 3,
}

impl StringEncoding {
    /// Create StringEncoding from u8 value
    #[inline]
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            0 => Some(StringEncoding::None),
            1 => Some(StringEncoding::Latin1),
            2 => Some(StringEncoding::Utf8),
            3 => Some(StringEncoding::Ascii),
            _ => None,
        }
    }

    /// Check if encoding is UTF-8 compatible
    #[inline]
    pub fn is_utf8_compatible(&self) -> bool {
        matches!(self, StringEncoding::Utf8 | StringEncoding::Ascii)
    }

    /// Get encoding name
    #[inline]
    pub fn name(&self) -> &'static str {
        match self {
            StringEncoding::None => "unknown",
            StringEncoding::Latin1 => "latin1",
            StringEncoding::Utf8 => "UTF-8",
            StringEncoding::Ascii => "ASCII",
        }
    }
}

impl std::fmt::Display for StringEncoding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encoding_values() {
        assert_eq!(StringEncoding::None as u8, 0);
        assert_eq!(StringEncoding::Latin1 as u8, 1);
        assert_eq!(StringEncoding::Utf8 as u8, 2);
        assert_eq!(StringEncoding::Ascii as u8, 3);
    }

    #[test]
    fn test_default_encoding() {
        assert_eq!(StringEncoding::default(), StringEncoding::Utf8);
    }

    #[test]
    fn test_from_u8() {
        assert_eq!(StringEncoding::from_u8(0), Some(StringEncoding::None));
        assert_eq!(StringEncoding::from_u8(1), Some(StringEncoding::Latin1));
        assert_eq!(StringEncoding::from_u8(2), Some(StringEncoding::Utf8));
        assert_eq!(StringEncoding::from_u8(3), Some(StringEncoding::Ascii));
        assert_eq!(StringEncoding::from_u8(4), None);
        assert_eq!(StringEncoding::from_u8(255), None);
    }

    #[test]
    fn test_is_utf8_compatible() {
        assert!(!StringEncoding::None.is_utf8_compatible());
        assert!(!StringEncoding::Latin1.is_utf8_compatible());
        assert!(StringEncoding::Utf8.is_utf8_compatible());
        assert!(StringEncoding::Ascii.is_utf8_compatible());
    }

    #[test]
    fn test_encoding_names() {
        assert_eq!(StringEncoding::None.name(), "unknown");
        assert_eq!(StringEncoding::Latin1.name(), "latin1");
        assert_eq!(StringEncoding::Utf8.name(), "UTF-8");
        assert_eq!(StringEncoding::Ascii.name(), "ASCII");
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", StringEncoding::Utf8), "UTF-8");
        assert_eq!(format!("{}", StringEncoding::Latin1), "latin1");
    }
}
