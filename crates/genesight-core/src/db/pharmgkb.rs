//! PharmGKB database adapter.
//!
//! Queries pharmacogenomic annotations (drug-gene interactions)
//! from the local PharmGKB SQLite table.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;
use crate::models::annotation::PharmaAnnotation;

/// Batch-query pharmacogenomic annotations for a list of rsIDs.
///
/// Returns a `HashMap` keyed by rsID. If multiple rows exist for a single
/// rsID (e.g., interactions with different drugs), only the first with the
/// highest evidence level is kept.
pub fn query_batch(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, PharmaAnnotation>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying PharmGKB");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    // Order by evidence_level so stronger evidence comes first
    let mut stmt = conn.prepare(
        "SELECT p.rsid, p.gene_symbol, p.drug, p.phenotype_category, \
         p.evidence_level, p.clinical_recommendation \
         FROM pharmacogenomics p \
         INNER JOIN tmp_rsids t ON p.rsid = t.rsid \
         ORDER BY p.evidence_level ASC",
    )?;

    let mut results = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            PharmaAnnotation {
                gene: row.get::<_, Option<String>>(1)?.unwrap_or_default(),
                drug: row.get(2)?,
                phenotype_category: row.get(3)?,
                evidence_level: row
                    .get::<_, Option<String>>(4)?
                    .unwrap_or_else(|| "4".to_string()),
                clinical_recommendation: row.get(5)?,
            },
        ))
    })?;

    for row in rows {
        let (rsid, annotation) = row?;
        // Keep the first (highest evidence level due to ORDER BY)
        results.entry(rsid).or_insert(annotation);
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(found = results.len(), "PharmGKB batch query complete");
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE pharmacogenomics (
                rsid TEXT,
                drug TEXT NOT NULL,
                phenotype_category TEXT,
                evidence_level TEXT,
                clinical_recommendation TEXT,
                gene_symbol TEXT
            );
            INSERT INTO pharmacogenomics VALUES (
                'rs1065852', 'Codeine', 'Poor Metabolizer', '1A',
                'Consider alternative analgesic', 'CYP2D6'
            );
            INSERT INTO pharmacogenomics VALUES (
                'rs9923231', 'Warfarin', 'Intermediate Metabolizer', '1A',
                'Reduce dose', 'VKORC1'
            );",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_pharma_annotations() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs1065852", "rs9923231"]).expect("query");
        assert_eq!(result.len(), 2);

        let ann = result.get("rs1065852").expect("rs1065852");
        assert_eq!(ann.gene, "CYP2D6");
        assert_eq!(ann.drug, "Codeine");
        assert_eq!(ann.evidence_level, "1A");
    }

    #[test]
    fn query_batch_empty_input() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }
}
