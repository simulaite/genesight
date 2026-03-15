//! GWAS Catalog database adapter.
//!
//! Queries genome-wide association study results for polygenic traits
//! from the local GWAS Catalog SQLite table. One rsID may have multiple
//! GWAS hits from different studies/traits.

use std::collections::HashMap;

use rusqlite::Connection;

use super::DbError;
use crate::models::annotation::GwasHit;

/// Batch-query GWAS Catalog hits for a list of rsIDs.
///
/// Returns a `HashMap` where each rsID maps to a `Vec<GwasHit>` because
/// a single SNP can be associated with multiple traits across studies.
pub fn query_batch(
    conn: &Connection,
    rsids: &[&str],
) -> Result<HashMap<String, Vec<GwasHit>>, DbError> {
    if rsids.is_empty() {
        return Ok(HashMap::new());
    }

    tracing::debug!(count = rsids.len(), "batch-querying GWAS Catalog");

    conn.execute_batch("CREATE TEMP TABLE IF NOT EXISTS tmp_rsids (rsid TEXT PRIMARY KEY)")?;
    conn.execute("DELETE FROM tmp_rsids", [])?;

    let mut insert = conn.prepare("INSERT OR IGNORE INTO tmp_rsids (rsid) VALUES (?1)")?;
    for rsid in rsids {
        insert.execute([rsid])?;
    }

    let mut stmt = conn.prepare(
        "SELECT g.rsid, g.trait, g.p_value, g.odds_ratio, g.beta, \
         g.risk_allele, g.risk_allele_frequency, g.pubmed_id, g.mapped_gene \
         FROM gwas g \
         INNER JOIN tmp_rsids t ON g.rsid = t.rsid",
    )?;

    let mut results: HashMap<String, Vec<GwasHit>> = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            GwasHit {
                trait_name: row.get(1)?,
                p_value: row.get::<_, Option<f64>>(2)?.unwrap_or(1.0),
                odds_ratio: row.get(3)?,
                beta: row.get(4)?,
                risk_allele: row.get(5)?,
                risk_allele_frequency: row.get(6)?,
                pubmed_id: row.get(7)?,
                mapped_gene: row.get(8)?,
            },
        ))
    })?;

    for row in rows {
        let (rsid, hit) = row?;
        results.entry(rsid).or_default().push(hit);
    }

    conn.execute("DROP TABLE IF EXISTS tmp_rsids", [])?;

    tracing::debug!(found = results.len(), "GWAS batch query complete");
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE gwas (
                rsid TEXT,
                trait TEXT NOT NULL,
                p_value REAL,
                odds_ratio REAL,
                beta REAL,
                risk_allele TEXT,
                risk_allele_frequency REAL,
                pubmed_id TEXT,
                mapped_gene TEXT
            );
            INSERT INTO gwas VALUES ('rs123', 'Type 2 Diabetes', 1e-9, 1.3, NULL, 'A', 0.35, '12345', 'TCF7L2');
            INSERT INTO gwas VALUES ('rs123', 'BMI', 5e-10, NULL, 0.05, 'G', 0.50, '67890', 'FTO');
            INSERT INTO gwas VALUES ('rs789', 'Height', 1e-20, NULL, 0.1, 'T', 0.45, '11111', 'HMGA2');",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn query_batch_returns_multiple_hits_per_rsid() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs123"]).expect("query");
        let hits = result.get("rs123").expect("rs123");
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn query_batch_multiple_rsids() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &["rs123", "rs789", "rs999"]).expect("query");
        assert_eq!(result.len(), 2);
        assert!(result.contains_key("rs123"));
        assert!(result.contains_key("rs789"));
    }

    #[test]
    fn query_batch_empty_input() {
        let conn = setup_test_db();
        let result = query_batch(&conn, &[]).expect("query");
        assert!(result.is_empty());
    }
}
