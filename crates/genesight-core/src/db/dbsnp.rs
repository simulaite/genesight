//! dbSNP / gnomAD allele frequency database adapter.
//!
//! Queries allele frequencies from the local `frequencies` table,
//! which is populated from gnomAD or dbSNP data.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;
use crate::models::annotation::AlleleFrequency;

/// Batch-query allele frequencies for a list of rsIDs.
///
/// Returns a `HashMap` keyed by rsID. If multiple rows exist for a single
/// rsID (e.g., from different sources), only the first encountered row is kept.
pub fn query_batch(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, AlleleFrequency>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying allele frequencies");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    let mut stmt = conn.prepare(
        "SELECT f.rsid, f.af_total, f.af_afr, f.af_amr, f.af_eas, f.af_eur, f.af_sas, f.source \
         FROM frequencies f \
         INNER JOIN tmp_rsids t ON f.rsid = t.rsid",
    )?;

    let mut results = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            AlleleFrequency {
                af_total: row.get::<_, Option<f64>>(1)?.unwrap_or(0.0),
                af_afr: row.get(2)?,
                af_amr: row.get(3)?,
                af_eas: row.get(4)?,
                af_eur: row.get(5)?,
                af_sas: row.get(6)?,
                source: row
                    .get::<_, Option<String>>(7)?
                    .unwrap_or_else(|| "unknown".to_string()),
            },
        ))
    })?;

    for row in rows {
        let (rsid, freq) = row?;
        // Keep the first entry if duplicates exist (gnomAD preferred over dbSNP)
        results.entry(rsid).or_insert(freq);
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(found = results.len(), "frequency batch query complete");
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE frequencies (
                rsid TEXT,
                af_total REAL,
                af_afr REAL,
                af_amr REAL,
                af_eas REAL,
                af_eur REAL,
                af_sas REAL,
                source TEXT
            );
            INSERT INTO frequencies VALUES ('rs123', 0.25, 0.30, 0.20, 0.15, 0.28, 0.22, 'gnomad');
            INSERT INTO frequencies VALUES ('rs456', 0.01, NULL, NULL, NULL, 0.02, NULL, 'dbsnp');",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_frequencies() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs123", "rs456"]).expect("query");
        assert_eq!(result.len(), 2);

        let freq = result.get("rs123").expect("rs123");
        assert!((freq.af_total - 0.25).abs() < f64::EPSILON);
        assert_eq!(freq.source, "gnomad");

        let freq2 = result.get("rs456").expect("rs456");
        assert!((freq2.af_total - 0.01).abs() < f64::EPSILON);
        assert!(freq2.af_afr.is_none());
    }

    #[test]
    fn query_batch_empty_input() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }
}
