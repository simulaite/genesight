//! JSON report renderer.

use crate::models::Report;

use super::ReportError;

/// Render a report as JSON.
pub fn render(report: &Report) -> Result<String, ReportError> {
    serde_json::to_string_pretty(report).map_err(|e| ReportError::Serialization(e.to_string()))
}
