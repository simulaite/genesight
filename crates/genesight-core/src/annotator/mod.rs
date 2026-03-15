//! Variant annotation engine.
//!
//! Matches parsed variants against local database entries to produce
//! annotated findings. The annotation step collects raw data from all
//! databases; scoring and interpretation happen in the [`scorer`](crate::scorer) module.

pub mod clinical;
pub mod frequency;
pub mod pharmacogenomics;
pub mod traits;

use rusqlite::Connection;

use crate::db::{self, DbError};
use crate::models::annotation::AnnotatedVariant;
use crate::models::config::AnnotationConfig;
use crate::models::variant::Variant;

/// Errors that can occur during variant annotation.
#[derive(Debug, thiserror::Error)]
pub enum AnnotateError {
    /// A database query failed.
    #[error("database error during annotation: {0}")]
    Database(#[from] DbError),

    /// A SQLite error occurred directly (e.g., from temp table operations).
    #[error("SQLite error: {0}")]
    Sqlite(#[from] rusqlite::Error),
}

/// Annotate a batch of variants against all available databases.
///
/// Queries ClinVar, GWAS Catalog, allele frequencies, and PharmGKB from the
/// main database connection. Optionally queries SNPedia from a separate connection.
///
/// Returns only variants that have at least one annotation (ClinVar, SNPedia,
/// GWAS hits, frequency data, or pharmacogenomics). Variants with no database
/// matches are silently dropped.
///
/// # Arguments
///
/// * `variants` - Parsed variants from a DNA raw data file
/// * `main_db` - Connection to `genesight.db` (ClinVar, GWAS, frequencies, PharmGKB)
/// * `snpedia_db` - Optional connection to `snpedia.db`
pub fn annotate_variants(
    variants: &[Variant],
    main_db: &Connection,
    snpedia_db: Option<&Connection>,
) -> Result<Vec<AnnotatedVariant>, AnnotateError> {
    annotate_variants_with_config(variants, main_db, snpedia_db, &AnnotationConfig::default())
}

/// Annotate a batch of variants with selective database queries.
///
/// Like [`annotate_variants`], but accepts an [`AnnotationConfig`] to control
/// which databases are queried. Disabled databases produce no annotations of
/// that type.
pub fn annotate_variants_with_config(
    variants: &[Variant],
    main_db: &Connection,
    snpedia_db: Option<&Connection>,
    config: &AnnotationConfig,
) -> Result<Vec<AnnotatedVariant>, AnnotateError> {
    // Collect rsIDs from variants that have them
    let rsid_variants: Vec<(&Variant, &str)> = variants
        .iter()
        .filter_map(|v| v.rsid.as_deref().map(|rsid| (v, rsid)))
        .collect();

    if rsid_variants.is_empty() {
        tracing::info!("no variants with rsIDs to annotate");
        return Ok(Vec::new());
    }

    let rsids: Vec<&str> = rsid_variants.iter().map(|(_, rsid)| *rsid).collect();

    tracing::info!(
        total_variants = variants.len(),
        with_rsid = rsids.len(),
        ?config,
        "starting batch annotation"
    );

    // Batch-query only enabled databases
    let clinvar_map = if config.clinvar {
        db::clinvar::query_batch(main_db, &rsids)?
    } else {
        std::collections::HashMap::new()
    };
    let gwas_map = if config.gwas {
        db::gwas::query_batch(main_db, &rsids)?
    } else {
        std::collections::HashMap::new()
    };
    let freq_map = if config.frequencies {
        db::dbsnp::query_batch(main_db, &rsids)?
    } else {
        std::collections::HashMap::new()
    };
    let pharma_map = if config.pharmacogenomics {
        db::pharmgkb::query_batch(main_db, &rsids)?
    } else {
        std::collections::HashMap::new()
    };

    let snpedia_map = match snpedia_db {
        Some(conn) => db::snpedia::query_batch(conn, &rsids)?,
        None => {
            tracing::debug!("SNPedia database not available, skipping");
            std::collections::HashMap::new()
        }
    };

    tracing::info!(
        clinvar = clinvar_map.len(),
        gwas = gwas_map.len(),
        frequency = freq_map.len(),
        pharmgkb = pharma_map.len(),
        snpedia = snpedia_map.len(),
        "database queries complete"
    );

    // Merge results into AnnotatedVariant structs
    let mut annotated = Vec::new();

    for (variant, rsid) in &rsid_variants {
        let clinvar = clinvar_map.get(*rsid).cloned();
        let snpedia = snpedia_map.get(*rsid).cloned();
        let gwas_hits = gwas_map.get(*rsid).cloned().unwrap_or_default();
        let frequency = freq_map.get(*rsid).cloned();
        let pharmacogenomics = pharma_map.get(*rsid).cloned();

        // Only include variants with at least one annotation
        let has_annotation = clinvar.is_some()
            || snpedia.is_some()
            || !gwas_hits.is_empty()
            || frequency.is_some()
            || pharmacogenomics.is_some();

        if has_annotation {
            annotated.push(AnnotatedVariant {
                variant: (*variant).clone(),
                clinvar,
                snpedia,
                gwas_hits,
                frequency,
                pharmacogenomics,
            });
        }
    }

    tracing::info!(annotated = annotated.len(), "annotation complete");

    Ok(annotated)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::variant::{Genotype, SourceFormat};

    fn setup_main_db() -> Connection {
        let conn = Connection::open_in_memory().expect("open in-memory db");
        conn.execute_batch(
            r#"
            CREATE TABLE clinvar (
                rsid TEXT, clinical_significance TEXT, review_status INTEGER,
                conditions TEXT, gene_symbol TEXT
            );
            CREATE TABLE gwas (
                rsid TEXT, trait TEXT NOT NULL, p_value REAL, odds_ratio REAL,
                beta REAL, risk_allele TEXT, risk_allele_frequency REAL,
                pubmed_id TEXT, mapped_gene TEXT
            );
            CREATE TABLE frequencies (
                rsid TEXT, af_total REAL, af_afr REAL, af_amr REAL,
                af_eas REAL, af_eur REAL, af_sas REAL, source TEXT
            );
            CREATE TABLE pharmacogenomics (
                rsid TEXT, drug TEXT NOT NULL, phenotype_category TEXT,
                evidence_level TEXT, clinical_recommendation TEXT, gene_symbol TEXT
            );

            INSERT INTO clinvar VALUES ('rs123', 'Pathogenic', 3, '["Breast cancer"]', 'BRCA1');
            INSERT INTO frequencies VALUES ('rs123', 0.001, NULL, NULL, NULL, 0.002, NULL, 'gnomad');
            INSERT INTO frequencies VALUES ('rs456', 0.50, NULL, NULL, NULL, NULL, NULL, 'dbsnp');
            "#,
        )
        .expect("setup");
        conn
    }

    fn make_variant(rsid: Option<&str>) -> Variant {
        Variant {
            rsid: rsid.map(String::from),
            chromosome: "1".to_string(),
            position: 100000,
            genotype: Genotype::Heterozygous('A', 'G'),
            source_format: SourceFormat::TwentyThreeAndMe,
        }
    }

    #[test]
    fn annotate_returns_only_annotated_variants() {
        let db = setup_main_db();
        let variants = vec![
            make_variant(Some("rs123")),
            make_variant(Some("rs999")), // not in any database
            make_variant(None),          // no rsid
        ];

        let result = annotate_variants(&variants, &db, None).expect("annotate");

        // rs123 has clinvar + frequency, rs456 has frequency only (not in input),
        // rs999 has nothing, None-rsid is skipped
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].variant.rsid.as_deref(), Some("rs123"));
        assert!(result[0].clinvar.is_some());
        assert!(result[0].frequency.is_some());
    }

    #[test]
    fn annotate_empty_variants() {
        let db = setup_main_db();
        let result = annotate_variants(&[], &db, None).expect("annotate");
        assert!(result.is_empty());
    }

    #[test]
    fn annotate_with_snpedia() {
        let main_db = setup_main_db();
        let snpedia_db = Connection::open_in_memory().expect("open snpedia");
        snpedia_db
            .execute_batch(
                r#"
                CREATE TABLE snpedia (
                    rsid TEXT PRIMARY KEY, magnitude REAL, repute TEXT,
                    summary TEXT, genotype_descriptions TEXT
                );
                INSERT INTO snpedia VALUES ('rs123', 3.0, 'bad', 'Associated with cancer risk', NULL);
                "#,
            )
            .expect("setup snpedia");

        let variants = vec![make_variant(Some("rs123"))];
        let result = annotate_variants(&variants, &main_db, Some(&snpedia_db)).expect("annotate");

        assert_eq!(result.len(), 1);
        assert!(result[0].snpedia.is_some());
        assert!(result[0].clinvar.is_some());
    }
}
