//! PGx allele definition loading from the database.
//!
//! Loads star-allele definitions, diplotype-phenotype mappings, and drug
//! recommendations from the `pgx_allele_definitions`, `pgx_diplotype_phenotypes`,
//! and `pgx_drug_recommendations` tables.

use std::collections::HashMap;

use rusqlite::Connection;

use crate::db::DbError;

/// A single star-allele defining variant.
#[derive(Debug, Clone)]
pub struct AlleleDefiningVariant {
    pub allele_name: String,
    pub rsid: String,
    pub alt_allele: String,
    pub function: String,
    pub activity_score: f64,
}

/// All allele definitions for a single gene.
#[derive(Debug, Clone)]
pub struct GeneAlleleDefinitions {
    /// Star allele name -> list of defining variants
    pub alleles: HashMap<String, Vec<AlleleDefiningVariant>>,
    /// All rsIDs that define any allele of this gene
    pub defining_rsids: Vec<String>,
}

/// Diplotype-to-phenotype mapping.
#[derive(Debug, Clone)]
pub struct DiplotypePhenotype {
    pub diplotype: String,
    pub phenotype: String,
    pub activity_score: Option<f64>,
}

/// Drug recommendation for a gene/phenotype combination.
#[derive(Debug, Clone)]
pub struct DrugRecommendation {
    pub drug: String,
    pub recommendation: String,
    pub evidence_level: String,
}

/// Load all PGx allele definitions from the database, grouped by gene.
pub fn load_allele_definitions(
    conn: &Connection,
) -> Result<HashMap<String, GeneAlleleDefinitions>, DbError> {
    // Check if the table exists
    let table_exists: bool = conn
        .prepare("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='pgx_allele_definitions'")?
        .query_row([], |row| row.get::<_, i64>(0))
        .map(|count| count > 0)?;

    if !table_exists {
        return Ok(HashMap::new());
    }

    let mut stmt = conn.prepare(
        "SELECT gene, allele_name, rsid, alt_allele, function, activity_score \
         FROM pgx_allele_definitions \
         ORDER BY gene, allele_name",
    )?;

    let mut gene_map: HashMap<String, GeneAlleleDefinitions> = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            AlleleDefiningVariant {
                allele_name: row.get(1)?,
                rsid: row.get(2)?,
                alt_allele: row.get(3)?,
                function: row.get(4)?,
                activity_score: row.get::<_, f64>(5)?,
            },
        ))
    })?;

    for row in rows {
        let (gene, variant) = row?;
        let entry = gene_map
            .entry(gene)
            .or_insert_with(|| GeneAlleleDefinitions {
                alleles: HashMap::new(),
                defining_rsids: Vec::new(),
            });

        if !entry.defining_rsids.contains(&variant.rsid) {
            entry.defining_rsids.push(variant.rsid.clone());
        }

        entry
            .alleles
            .entry(variant.allele_name.clone())
            .or_default()
            .push(variant);
    }

    Ok(gene_map)
}

/// Load diplotype-to-phenotype mappings for all genes.
pub fn load_diplotype_phenotypes(
    conn: &Connection,
) -> Result<HashMap<String, Vec<DiplotypePhenotype>>, DbError> {
    let table_exists: bool = conn
        .prepare("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='pgx_diplotype_phenotypes'")?
        .query_row([], |row| row.get::<_, i64>(0))
        .map(|count| count > 0)?;

    if !table_exists {
        return Ok(HashMap::new());
    }

    let mut stmt = conn.prepare(
        "SELECT gene, diplotype, phenotype, activity_score \
         FROM pgx_diplotype_phenotypes \
         ORDER BY gene",
    )?;

    let mut result: HashMap<String, Vec<DiplotypePhenotype>> = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            DiplotypePhenotype {
                diplotype: row.get(1)?,
                phenotype: row.get(2)?,
                activity_score: row.get(3)?,
            },
        ))
    })?;

    for row in rows {
        let (gene, dp) = row?;
        result.entry(gene).or_default().push(dp);
    }

    Ok(result)
}

/// Load drug recommendations for all gene/phenotype combinations.
pub fn load_drug_recommendations(
    conn: &Connection,
) -> Result<HashMap<(String, String), Vec<DrugRecommendation>>, DbError> {
    let table_exists: bool = conn
        .prepare("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='pgx_drug_recommendations'")?
        .query_row([], |row| row.get::<_, i64>(0))
        .map(|count| count > 0)?;

    if !table_exists {
        return Ok(HashMap::new());
    }

    let mut stmt = conn.prepare(
        "SELECT gene, phenotype, drug, recommendation, evidence_level \
         FROM pgx_drug_recommendations \
         ORDER BY gene, phenotype",
    )?;

    let mut result: HashMap<(String, String), Vec<DrugRecommendation>> = HashMap::new();

    let rows = stmt.query_map([], |row| {
        Ok((
            row.get::<_, String>(0)?,
            row.get::<_, String>(1)?,
            DrugRecommendation {
                drug: row.get(2)?,
                recommendation: row.get(3)?,
                evidence_level: row.get(4)?,
            },
        ))
    })?;

    for row in rows {
        let (gene, phenotype, rec) = row?;
        result.entry((gene, phenotype)).or_default().push(rec);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_from_empty_db() {
        let conn = Connection::open_in_memory().expect("open");
        let defs = load_allele_definitions(&conn).expect("load");
        assert!(defs.is_empty());
    }

    #[test]
    fn load_definitions_from_populated_db() {
        let conn = Connection::open_in_memory().expect("open");
        conn.execute_batch(
            "CREATE TABLE pgx_allele_definitions (
                gene TEXT NOT NULL,
                allele_name TEXT NOT NULL,
                rsid TEXT NOT NULL,
                alt_allele TEXT NOT NULL,
                function TEXT NOT NULL,
                activity_score REAL NOT NULL,
                PRIMARY KEY (gene, allele_name, rsid)
            );
            INSERT INTO pgx_allele_definitions VALUES ('CYP2C19', '*2', 'rs4244285', 'A', 'No Function', 0.0);
            INSERT INTO pgx_allele_definitions VALUES ('CYP2C19', '*17', 'rs12248560', 'T', 'Increased Function', 1.5);",
        ).expect("setup");

        let defs = load_allele_definitions(&conn).expect("load");
        assert_eq!(defs.len(), 1);
        let cyp2c19 = defs.get("CYP2C19").expect("CYP2C19");
        assert_eq!(cyp2c19.alleles.len(), 2);
        assert!(cyp2c19.defining_rsids.contains(&"rs4244285".to_string()));
        assert!(cyp2c19.defining_rsids.contains(&"rs12248560".to_string()));
    }
}
