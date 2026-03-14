//! Report generation in various output formats.

pub mod html;
pub mod json;
pub mod markdown;

use crate::models::Report;

/// Errors that can occur during report generation.
#[derive(Debug, thiserror::Error)]
pub enum ReportError {
    #[error("Serialization error: {0}")]
    Serialization(String),
}

/// Output format for reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Markdown,
    Json,
    Html,
}

/// Render a report in the specified format.
pub fn render(report: &Report, format: OutputFormat) -> Result<String, ReportError> {
    match format {
        OutputFormat::Markdown => markdown::render(report),
        OutputFormat::Json => json::render(report),
        OutputFormat::Html => html::render(report),
    }
}
