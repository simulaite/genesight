//! Pharmacogenomic star allele calling and phenotype classification.
//!
//! This module implements CPIC-style star allele matching from genotype data
//! and translates diplotypes into metabolizer phenotypes with activity scores.
//!
//! # Overview
//!
//! The PGx pipeline works in three stages:
//! 1. **Star allele calling** — Match observed genotypes to known star allele definitions
//! 2. **Phenotype classification** — Translate diplotypes into phenotypes via activity scores
//! 3. **Drug recommendation lookup** — Map phenotypes to clinical recommendations
//!
//! # Limitations
//!
//! - Microarray data is unphased; diplotype inference is approximate
//! - Only single-variant star alleles are supported (no multi-variant haplotypes)
//! - Gene deletion/duplication (e.g., CYP2D6 CNV) is not detected

pub mod phenotype;

use std::collections::HashMap;

use rusqlite::Connection;

use crate::db::DbError;

/// A star allele definition from the database.
#[derive(Debug, Clone)]
pub struct AlleleDefinition {
    /// Gene symbol (e.g., "CYP2D6")
    pub gene: String,
    /// Star allele name (e.g., "*4")
    pub allele_name: String,
    /// rsID for the defining variant
    pub rsid: String,
    /// Alternate allele that defines this star allele
    pub alt_allele: String,
    /// Functional classification (e.g., "No Function")
    pub function: String,
    /// Numeric activity score
    pub activity_score: f64,
}

/// Result of star allele calling for a single gene.
#[derive(Debug, Clone)]
pub struct StarAlleleCall {
    /// Gene symbol
    pub gene: String,
    /// Called diplotype (e.g., "*1/*4")
    pub diplotype: String,
    /// Individual allele calls (two alleles)
    pub alleles: (String, String),
    /// Total activity score (sum of both alleles)
    pub activity_score: f64,
    /// Classified phenotype (e.g., "Poor Metabolizer")
    pub phenotype: String,
    /// Limitations or caveats about this call
    pub limitations: Vec<String>,
}

/// Star allele caller that matches genotypes to known allele definitions.
///
/// Loads allele definitions from the database and matches them against
/// observed variant genotypes to infer diplotypes and phenotypes.
pub struct StarAlleleCaller {
    /// Allele definitions grouped by gene
    definitions_by_gene: HashMap<String, Vec<AlleleDefinition>>,
}

impl StarAlleleCaller {
    /// Create a new `StarAlleleCaller` by loading definitions from the database.
    ///
    /// # Errors
    ///
    /// Returns `DbError` if the database query fails or the
    /// `pgx_allele_definitions` table does not exist.
    pub fn from_db(conn: &Connection) -> Result<Self, DbError> {
        let mut stmt = conn.prepare(
            "SELECT gene, allele_name, rsid, alt_allele, function, activity_score \
             FROM pgx_allele_definitions \
             ORDER BY gene, allele_name",
        )?;

        let mut definitions_by_gene: HashMap<String, Vec<AlleleDefinition>> = HashMap::new();

        let rows = stmt.query_map([], |row| {
            Ok(AlleleDefinition {
                gene: row.get(0)?,
                allele_name: row.get(1)?,
                rsid: row.get(2)?,
                alt_allele: row.get(3)?,
                function: row.get(4)?,
                activity_score: row.get(5)?,
            })
        })?;

        for row in rows {
            let def = row?;
            definitions_by_gene
                .entry(def.gene.clone())
                .or_default()
                .push(def);
        }

        tracing::info!(
            genes = definitions_by_gene.len(),
            "loaded PGx allele definitions"
        );

        Ok(Self {
            definitions_by_gene,
        })
    }

    /// Return the set of gene names that have allele definitions.
    pub fn supported_genes(&self) -> Vec<&str> {
        let mut genes: Vec<&str> = self
            .definitions_by_gene
            .keys()
            .map(|s| s.as_str())
            .collect();
        genes.sort();
        genes
    }

    /// Call star alleles for a specific gene given observed genotypes.
    ///
    /// `genotypes` maps rsID -> observed genotype string (e.g., "AG", "CC").
    /// Returns `None` if no allele definitions exist for the gene or if
    /// no relevant variants were observed.
    pub fn call_gene(
        &self,
        gene: &str,
        genotypes: &HashMap<String, String>,
    ) -> Option<StarAlleleCall> {
        let definitions = self.definitions_by_gene.get(gene)?;

        // Count how many alternate alleles the user carries for each star allele
        let mut allele_alt_counts: HashMap<&str, u8> = HashMap::new();
        let mut allele_scores: HashMap<&str, f64> = HashMap::new();
        let mut any_variant_observed = false;

        for def in definitions {
            if let Some(observed) = genotypes.get(&def.rsid) {
                any_variant_observed = true;
                let alt_char = def.alt_allele.chars().next()?;
                let count = observed.chars().filter(|&c| c == alt_char).count() as u8;

                // For alleles defined by multiple variants (e.g., TPMT *3A),
                // we track the minimum count across all defining variants
                let entry = allele_alt_counts.entry(&def.allele_name).or_insert(count);
                if count < *entry {
                    *entry = count;
                }
                allele_scores.insert(&def.allele_name, def.activity_score);
            }
        }

        if !any_variant_observed {
            return None;
        }

        // Determine the two alleles: for each star allele with count >= 1,
        // assign one copy; for count >= 2, assign two copies
        let mut called_alleles: Vec<(&str, f64)> = Vec::new();

        for (allele_name, count) in &allele_alt_counts {
            let score = allele_scores.get(allele_name).copied().unwrap_or(1.0);
            match count {
                2 => {
                    called_alleles.push((allele_name, score));
                    called_alleles.push((allele_name, score));
                }
                1 => {
                    called_alleles.push((allele_name, score));
                }
                _ => {}
            }
        }

        // Fill remaining slots with *1 (reference/normal function)
        let default_allele = "*1";
        let default_score = 1.0;
        while called_alleles.len() < 2 {
            called_alleles.push((default_allele, default_score));
        }

        // If more than 2 alleles were called (possible with multiple star alleles),
        // keep only the two with lowest activity scores (most impactful)
        called_alleles.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        called_alleles.truncate(2);

        // Sort alphabetically for canonical diplotype representation
        called_alleles.sort_by(|a, b| a.0.cmp(b.0));

        let allele1 = called_alleles[0].0.to_string();
        let allele2 = called_alleles[1].0.to_string();
        let total_score = called_alleles[0].1 + called_alleles[1].1;

        let diplotype = format!("{allele1}/{allele2}");
        let phenotype_result = phenotype::call_phenotype(gene, total_score);

        Some(StarAlleleCall {
            gene: gene.to_string(),
            diplotype,
            alleles: (allele1, allele2),
            activity_score: total_score,
            phenotype: phenotype_result.phenotype,
            limitations: phenotype_result.limitations,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_test_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            "CREATE TABLE pgx_allele_definitions (
                gene TEXT NOT NULL,
                allele_name TEXT NOT NULL,
                rsid TEXT NOT NULL,
                alt_allele TEXT NOT NULL,
                function TEXT NOT NULL,
                activity_score REAL NOT NULL
            );
            INSERT INTO pgx_allele_definitions VALUES
                ('CYP2D6', '*4', 'rs3892097', 'A', 'No Function', 0.0);
            INSERT INTO pgx_allele_definitions VALUES
                ('CYP2C19', '*2', 'rs4244285', 'A', 'No Function', 0.0);
            INSERT INTO pgx_allele_definitions VALUES
                ('CYP2C19', '*17', 'rs12248560', 'T', 'Increased Function', 1.5);",
        )
        .expect("setup");
        conn
    }

    #[test]
    fn from_db_loads_definitions() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");
        let genes = caller.supported_genes();
        assert!(genes.contains(&"CYP2D6"));
        assert!(genes.contains(&"CYP2C19"));
    }

    #[test]
    fn call_gene_homozygous_alt() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");

        let mut genotypes = HashMap::new();
        genotypes.insert("rs3892097".to_string(), "AA".to_string());

        let result = caller.call_gene("CYP2D6", &genotypes).expect("call");
        assert_eq!(result.diplotype, "*4/*4");
        assert!((result.activity_score - 0.0).abs() < f64::EPSILON);
        assert!(result.phenotype.contains("Poor"));
    }

    #[test]
    fn call_gene_heterozygous() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");

        let mut genotypes = HashMap::new();
        genotypes.insert("rs3892097".to_string(), "GA".to_string());

        let result = caller.call_gene("CYP2D6", &genotypes).expect("call");
        assert_eq!(result.diplotype, "*1/*4");
        assert!((result.activity_score - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn call_gene_homozygous_ref() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");

        let mut genotypes = HashMap::new();
        genotypes.insert("rs3892097".to_string(), "GG".to_string());

        let result = caller.call_gene("CYP2D6", &genotypes).expect("call");
        assert_eq!(result.diplotype, "*1/*1");
        assert!((result.activity_score - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn call_gene_no_data_returns_none() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");

        let genotypes = HashMap::new();
        assert!(caller.call_gene("CYP2D6", &genotypes).is_none());
    }

    #[test]
    fn call_gene_unknown_gene_returns_none() {
        let conn = setup_test_db();
        let caller = StarAlleleCaller::from_db(&conn).expect("load");

        let genotypes = HashMap::new();
        assert!(caller.call_gene("UNKNOWN_GENE", &genotypes).is_none());
    }
}
