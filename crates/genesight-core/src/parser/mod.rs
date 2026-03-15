//! DNA raw data file parsers.
//!
//! Supports 23andMe, AncestryDNA, and VCF formats.
//! All parsers take `&str` content — no filesystem IO in core.

pub mod ancestry;
pub mod twentythreeandme;
pub mod vcf;

use crate::models::{GenomeAssembly, Variant};

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

/// A parsed DNA file including both variants and detected metadata.
#[derive(Debug, Clone)]
pub struct ParsedFile {
    /// Parsed genetic variants.
    pub variants: Vec<Variant>,
    /// Detected genome assembly from the file headers.
    pub assembly: GenomeAssembly,
}

/// Parse DNA content with auto-detected format, returning variants and metadata.
///
/// This is the recommended entry point when assembly information is needed.
/// For backwards compatibility, see [`parse_auto`] which discards metadata.
pub fn parse_auto_with_metadata(content: &str) -> Result<ParsedFile, ParseError> {
    let format = detect_format(content)?;
    let (variants, assembly) = match format {
        FileFormat::TwentyThreeAndMe => (
            twentythreeandme::parse(content)?,
            twentythreeandme::detect_assembly(content),
        ),
        FileFormat::AncestryDna => (
            ancestry::parse(content)?,
            ancestry::detect_assembly(content),
        ),
        FileFormat::Vcf => (vcf::parse(content)?, vcf::detect_assembly(content)),
    };
    Ok(ParsedFile { variants, assembly })
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
