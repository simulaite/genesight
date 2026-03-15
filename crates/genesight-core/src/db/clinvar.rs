//! ClinVar database adapter.
//!
//! Queries clinically classified variants (pathogenic/benign) from the local
//! ClinVar SQLite table. Uses a temporary table for efficient batch lookups.
//!
//! Since ClinVar 2024, entries may carry a `classification_type` column
//! (`germline`, `somatic`, `oncogenicity`). When multiple entries exist for the
//! same rsID with different classification types, germline entries take
//! precedence because this is a consumer (germline) DNA analysis tool.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;
use crate::models::annotation::{ClinVarAnnotation, ClinVarClassificationType};

/// Batch-query ClinVar annotations for a list of rsIDs.
///
/// Creates a temporary table, inserts all requested rsIDs, then JOINs against
/// the `clinvar` table to retrieve matching annotations in a single query.
/// Returns a `HashMap` keyed by rsID.
///
/// When multiple ClinVar entries exist for the same rsID with different
/// classification types, germline entries are preferred over somatic or
/// oncogenicity entries (since this tool analyzes consumer germline DNA).
///
/// The `classification_type` column is read with a fallback: if the column
/// does not exist in the database (older schema), all entries default to
/// `Germline`.
pub fn query_batch(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, ClinVarAnnotation>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying ClinVar");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    // Detect whether the classification_type column exists for backward
    // compatibility with databases created before this column was added.
    let has_classification_type = has_column(conn, "clinvar", "classification_type");

    let query = if has_classification_type {
        "SELECT c.rsid, c.clinical_significance, c.review_status, c.conditions, \
         c.gene_symbol, c.classification_type \
         FROM clinvar c \
         INNER JOIN tmp_rsids t ON c.rsid = t.rsid"
    } else {
        "SELECT c.rsid, c.clinical_significance, c.review_status, c.conditions, \
         c.gene_symbol, NULL AS classification_type \
         FROM clinvar c \
         INNER JOIN tmp_rsids t ON c.rsid = t.rsid"
    };

    let mut stmt = conn.prepare(query)?;

    let mut results: HashMap<String, ClinVarAnnotation> = HashMap::new();

    let rows = stmt.query_map([], |row| {
        let rsid: String = row.get(0)?;
        let significance: String = row.get::<_, Option<String>>(1)?.unwrap_or_default();
        let review_stars: u8 = row.get::<_, Option<u8>>(2)?.unwrap_or(0);
        let conditions_json: Option<String> = row.get(3)?;
        let gene_symbol: Option<String> = row.get(4)?;
        let classification_type_str: Option<String> = row.get(5)?;
        Ok((
            rsid,
            significance,
            review_stars,
            conditions_json,
            gene_symbol,
            classification_type_str,
        ))
    })?;

    for row in rows {
        let (rsid, significance, review_stars, conditions_json, gene_symbol, ct_str) = row?;

        let conditions = parse_json_string_array(conditions_json.as_deref());
        let classification_type =
            ClinVarClassificationType::from_db_str(ct_str.as_deref().unwrap_or("germline"));

        let annotation = ClinVarAnnotation {
            significance,
            review_stars,
            conditions,
            gene_symbol,
            classification_type,
        };

        // When multiple entries exist for the same rsID, prefer germline.
        match results.get(&rsid) {
            Some(existing)
                if existing.classification_type == ClinVarClassificationType::Germline =>
            {
                // Already have a germline entry — keep it.
            }
            _ => {
                // Either no entry yet, or existing is non-germline and we may upgrade.
                results.insert(rsid, annotation);
            }
        }
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(found = results.len(), "ClinVar batch query complete");
    Ok(results)
}

/// Check whether a table has a specific column, for backward compatibility.
fn has_column(conn: &Connection, table: &str, column: &str) -> bool {
    let mut stmt = match conn.prepare(&format!("PRAGMA table_info({})", table)) {
        Ok(s) => s,
        Err(_) => return false,
    };
    stmt.query_map([], |row| row.get::<_, String>(1))
        .ok()
        .map(|rows| rows.flatten().any(|name| name == column))
        .unwrap_or(false)
}

/// Parse a JSON string containing an array of strings, returning an empty vec on failure.
fn parse_json_string_array(json: Option<&str>) -> Vec<String> {
    match json {
        Some(s) => serde_json::from_str::<Vec<String>>(s).unwrap_or_default(),
        None => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a ClinVar test DB with the classification_type column.
    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE clinvar (
                rsid TEXT,
                clinical_significance TEXT,
                review_status INTEGER,
                conditions TEXT,
                gene_symbol TEXT,
                classification_type TEXT DEFAULT 'germline'
            );
            INSERT INTO clinvar VALUES ('rs123', 'Pathogenic', 3, '[\"Breast cancer\"]', 'BRCA1', 'germline');
            INSERT INTO clinvar VALUES ('rs456', 'Benign', 1, '[\"Unspecified\"]', 'TP53', 'germline');",
        )
        .expect("setup");
        conn
    }

    /// Helper: create a ClinVar test DB without the classification_type column
    /// (simulates an older database schema for backward compatibility).
    fn setup_test_db_legacy() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE clinvar (
                rsid TEXT,
                clinical_significance TEXT,
                review_status INTEGER,
                conditions TEXT,
                gene_symbol TEXT
            );
            INSERT INTO clinvar VALUES ('rs123', 'Pathogenic', 3, '[\"Breast cancer\"]', 'BRCA1');
            INSERT INTO clinvar VALUES ('rs456', 'Benign', 1, '[\"Unspecified\"]', 'TP53');",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_matching_entries() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs123", "rs999"]).expect("query");
        assert_eq!(result.len(), 1);
        let ann = result.get("rs123").expect("rs123");
        assert_eq!(ann.significance, "Pathogenic");
        assert_eq!(ann.review_stars, 3);
        assert_eq!(ann.conditions, vec!["Breast cancer".to_string()]);
        assert_eq!(ann.gene_symbol.as_deref(), Some("BRCA1"));
        assert_eq!(ann.classification_type, ClinVarClassificationType::Germline);
    }

    #[test]
    fn query_batch_empty_input() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }

    #[test]
    fn query_batch_legacy_db_defaults_to_germline() {
        let conn = setup_test_db_legacy();
        let result = query_batch(&conn, &["rs123"]).expect("query");
        assert_eq!(result.len(), 1);
        let ann = result.get("rs123").expect("rs123");
        assert_eq!(ann.classification_type, ClinVarClassificationType::Germline);
    }

    #[test]
    fn query_batch_reads_somatic_classification() {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE clinvar (
                rsid TEXT,
                clinical_significance TEXT,
                review_status INTEGER,
                conditions TEXT,
                gene_symbol TEXT,
                classification_type TEXT DEFAULT 'germline'
            );
            INSERT INTO clinvar VALUES ('rs100', 'Pathogenic', 2, '[\"Lung cancer\"]', 'EGFR', 'somatic');",
        )
        .expect("setup");

        let result = query_batch(&conn, &["rs100"]).expect("query");
        let ann = result.get("rs100").expect("rs100");
        assert_eq!(ann.classification_type, ClinVarClassificationType::Somatic);
    }

    #[test]
    fn query_batch_germline_preferred_over_somatic() {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE clinvar (
                rsid TEXT,
                clinical_significance TEXT,
                review_status INTEGER,
                conditions TEXT,
                gene_symbol TEXT,
                classification_type TEXT DEFAULT 'germline'
            );
            INSERT INTO clinvar VALUES ('rs200', 'Pathogenic', 2, '[\"Tumor X\"]', 'TP53', 'somatic');
            INSERT INTO clinvar VALUES ('rs200', 'Likely pathogenic', 3, '[\"Li-Fraumeni syndrome\"]', 'TP53', 'germline');",
        )
        .expect("setup");

        let result = query_batch(&conn, &["rs200"]).expect("query");
        let ann = result.get("rs200").expect("rs200");
        assert_eq!(ann.classification_type, ClinVarClassificationType::Germline);
        assert_eq!(ann.significance, "Likely pathogenic");
    }

    #[test]
    fn parse_json_handles_null() {
        assert!(parse_json_string_array(None).is_empty());
    }

    #[test]
    fn parse_json_handles_malformed() {
        assert!(parse_json_string_array(Some("not json")).is_empty());
    }

    #[test]
    fn classification_type_from_db_str() {
        assert_eq!(
            ClinVarClassificationType::from_db_str("germline"),
            ClinVarClassificationType::Germline
        );
        assert_eq!(
            ClinVarClassificationType::from_db_str("somatic"),
            ClinVarClassificationType::Somatic
        );
        assert_eq!(
            ClinVarClassificationType::from_db_str("oncogenicity"),
            ClinVarClassificationType::Oncogenicity
        );
        assert_eq!(
            ClinVarClassificationType::from_db_str("SOMATIC"),
            ClinVarClassificationType::Somatic
        );
        // Unknown defaults to Germline
        assert_eq!(
            ClinVarClassificationType::from_db_str("unknown"),
            ClinVarClassificationType::Germline
        );
    }
}
