//! DNA raw data file parsers.
//!
//! Supports 23andMe, AncestryDNA, and VCF formats.
//! All parsers take `&str` content — no filesystem IO in core.

pub mod ancestry;
pub mod twentythreeandme;
pub mod vcf;

use crate::models::Variant;

/// Supported DNA file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    TwentyThreeAndMe,
    AncestryDna,
    Vcf,
}

/// Errors that can occur during DNA file parsing.
#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("Unrecognized file format")]
    UnrecognizedFormat,
    #[error("Invalid line {line}: {reason}")]
    InvalidLine { line: usize, reason: String },
}

/// Detect the file format by examining the content header.
pub fn detect_format(content: &str) -> Result<FileFormat, ParseError> {
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            if trimmed.starts_with("##fileformat=VCF") {
                return Ok(FileFormat::Vcf);
            }
            continue;
        }
        // Check first data line column count
        let cols: Vec<&str> = trimmed.split('\t').collect();
        return match cols.len() {
            4 => Ok(FileFormat::TwentyThreeAndMe),
            5 => Ok(FileFormat::AncestryDna),
            _ if cols.len() >= 8 => Ok(FileFormat::Vcf),
            _ => Err(ParseError::UnrecognizedFormat),
        };
    }
    Err(ParseError::UnrecognizedFormat)
}

/// Parse DNA content with auto-detected format.
pub fn parse_auto(content: &str) -> Result<Vec<Variant>, ParseError> {
    let format = detect_format(content)?;
    match format {
        FileFormat::TwentyThreeAndMe => twentythreeandme::parse(content),
        FileFormat::AncestryDna => ancestry::parse(content),
        FileFormat::Vcf => vcf::parse(content),
    }
}
