//! SNPedia database adapter.
//!
//! Queries SNP annotations, magnitude scores, and human-readable summaries
//! from the optional SNPedia SQLite database (`snpedia.db`).
//!
//! This database is separate from the main `genesight.db` because SNPedia
//! data is licensed under CC-BY-NC-SA 3.0.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;
use crate::models::annotation::SnpediaAnnotation;

/// Batch-query SNPedia annotations for a list of rsIDs.
///
/// The `conn` parameter must point to the separate `snpedia.db` database.
/// Returns a `HashMap` keyed by rsID.
pub fn query_batch(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, SnpediaAnnotation>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying SNPedia");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    let mut stmt = conn.prepare(
        "SELECT s.rsid, s.magnitude, s.repute, s.summary, s.genotype_descriptions \
         FROM snpedia s \
         INNER JOIN tmp_rsids t ON s.rsid = t.rsid",
    )?;

    let mut results = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            row.get::<_, Option<f64>>(1)?,
            row.get::<_, Option<String>>(2)?,
            row.get::<_, Option<String>>(3)?,
            row.get::<_, Option<String>>(4)?,
        ))
    })?;

    for row in rows {
        let (rsid, magnitude, repute, summary, genotype_json) = row?;

        let genotype_descriptions = parse_genotype_descriptions(genotype_json.as_deref());

        results.insert(
            rsid,
            SnpediaAnnotation {
                magnitude: magnitude.unwrap_or(0.0),
                repute,
                summary: summary.unwrap_or_default(),
                genotype_descriptions,
            },
        );
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(found = results.len(), "SNPedia batch query complete");
    Ok(results)
}

/// Parse genotype descriptions from a JSON string into a HashMap.
///
/// Expected format: `{"AA": "description", "AG": "description", ...}`
/// Returns `None` if the JSON is absent or malformed.
fn parse_genotype_descriptions(json: Option<&str>) -> Option<HashMap<String, String>> {
    let s = json?;
    serde_json::from_str::<HashMap<String, String>>(s).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            r#"CREATE TABLE snpedia (
                rsid TEXT PRIMARY KEY,
                magnitude REAL,
                repute TEXT,
                summary TEXT,
                genotype_descriptions TEXT
            );
            INSERT INTO snpedia VALUES (
                'rs4988235', 3.0, 'good', 'Lactase persistence',
                '{"TT": "Likely lactose tolerant", "CT": "Likely lactose tolerant", "CC": "Likely lactose intolerant"}'
            );
            INSERT INTO snpedia VALUES (
                'rs1805007', 4.0, 'bad', 'Red hair pigmentation',
                '{"TT": "Red hair likely", "CT": "Carrier for red hair"}'
            );"#,
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_annotations() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs4988235", "rs1805007"]).expect("query");
        assert_eq!(result.len(), 2);

        let ann = result.get("rs4988235").expect("rs4988235");
        assert!((ann.magnitude - 3.0).abs() < f64::EPSILON);
        assert_eq!(ann.repute.as_deref(), Some("good"));
        assert_eq!(ann.summary, "Lactase persistence");

        let descs = ann.genotype_descriptions.as_ref().expect("descriptions");
        assert_eq!(
            descs.get("TT").map(String::as_str),
            Some("Likely lactose tolerant")
        );
    }

    #[test]
    fn query_batch_empty_input() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }

    #[test]
    fn parse_genotype_descriptions_handles_none() {
        assert!(parse_genotype_descriptions(None).is_none());
    }

    #[test]
    fn parse_genotype_descriptions_handles_malformed() {
        assert!(parse_genotype_descriptions(Some("not json")).is_none());
    }
}
