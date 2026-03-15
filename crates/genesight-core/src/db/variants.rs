//! Variants table adapter for fetching reference/alternate alleles.
//!
//! The `variants` table stores the ref/alt alleles for each rsID. This data
//! is essential for allele matching: determining whether a user actually
//! carries a variant allele vs. being homozygous reference.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;

/// Reference and alternate alleles for a variant position.
#[derive(Debug, Clone)]
pub struct VariantAlleles {
    /// Reference allele (e.g., "A", "G")
    pub ref_allele: String,
    /// Alternate allele (e.g., "T", "C")
    pub alt_allele: String,
}

/// Batch-query reference/alternate alleles for a list of rsIDs.
///
/// Returns a `HashMap` keyed by rsID. Only variants with both ref and alt
/// alleles populated are included.
pub fn query_batch_alleles(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, VariantAlleles>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    // Check if the variants table exists
    let has_table: bool = conn
        .query_row(
            "SELECT COUNT(*) > 0 FROM sqlite_master WHERE type='table' AND name='variants'",
            [],
            |row| row.get(0),
        )
        .unwrap_or(false);

    if !has_table {
        tracing::debug!("variants table not found, skipping allele lookup");
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying variant alleles");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    let mut stmt = conn.prepare(
        "SELECT v.rsid, v.ref_allele, v.alt_allele \
         FROM variants v \
         INNER JOIN tmp_rsids t ON v.rsid = t.rsid \
         WHERE v.ref_allele IS NOT NULL AND v.alt_allele IS NOT NULL",
    )?;

    let mut results = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            row.get::<_, String>(1)?,
            row.get::<_, String>(2)?,
        ))
    })?;

    for row in rows {
        let (rsid, ref_allele, alt_allele) = row?;
        if !ref_allele.is_empty() && !alt_allele.is_empty() {
            results.insert(
                rsid,
                VariantAlleles {
                    ref_allele,
                    alt_allele,
                },
            );
        }
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(
        found = results.len(),
        "variant alleles batch query complete"
    );
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE variants (
                rsid TEXT PRIMARY KEY,
                chromosome TEXT NOT NULL,
                position INTEGER NOT NULL,
                ref_allele TEXT,
                alt_allele TEXT
            );
            INSERT INTO variants VALUES ('rs123', '17', 43044295, 'A', 'G');
            INSERT INTO variants VALUES ('rs456', '7', 117559593, 'C', 'T');
            INSERT INTO variants VALUES ('rs789', '1', 100000, NULL, NULL);",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_alleles() {
        let conn = setup_test_db();
        let result = query_batch_alleles(&conn, &["rs123", "rs456"]).expect("query");
        assert_eq!(result.len(), 2);

        let alleles = result.get("rs123").expect("rs123");
        assert_eq!(alleles.ref_allele, "A");
        assert_eq!(alleles.alt_allele, "G");
    }

    #[test]
    fn null_alleles_are_excluded() {
        let conn = setup_test_db();
        let result = query_batch_alleles(&conn, &["rs789"]).expect("query");
        assert!(result.is_empty(), "NULL alleles should be excluded");
    }

    #[test]
    fn missing_rsid_returns_empty() {
        let conn = setup_test_db();
        let result = query_batch_alleles(&conn, &["rs999"]).expect("query");
        assert!(result.is_empty());
    }

    #[test]
    fn empty_input_returns_empty() {
        let conn = setup_test_db();
        let result = query_batch_alleles(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }

    #[test]
    fn no_variants_table_returns_empty() {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        let result = query_batch_alleles(&conn, &["rs123"]).expect("query");
        assert!(result.is_empty());
    }
}
