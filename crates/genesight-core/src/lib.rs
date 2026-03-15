//! GeneSight Core Library
//!
//! Privacy-first DNA analysis engine that annotates personal genetic data
//! against local copies of public genome databases.
//!
//! # Architecture
//!
//! - **parser** — Read DNA raw data files (23andMe, AncestryDNA, VCF)
//! - **db** — Query local SQLite databases (ClinVar, SNPedia, GWAS, dbSNP, PharmGKB)
//! - **annotator** — Match variants against database entries
//! - **scorer** — Assign confidence tiers and risk scores
//! - **report** — Generate human-readable reports (Markdown, JSON, HTML)
//! - **models** — Shared types and data structures
//!
//! # Usage
//!
//! The primary entry point is the [`analyze`] function, which runs the full
//! pipeline: annotate variants, score them, filter by requested tiers, and
//! build a report.
//!
//! ```rust,no_run
//! use genesight_core::models::{Variant, ConfidenceTier};
//! use rusqlite::Connection;
//!
//! let main_db = Connection::open("genesight.db").unwrap();
//! let variants: Vec<Variant> = vec![]; // parsed from DNA file
//! let tiers = [ConfidenceTier::Tier1Reliable, ConfidenceTier::Tier2Probable];
//!
//! let report = genesight_core::analyze(&variants, &main_db, None, &tiers).unwrap();
//! ```

pub mod annotator;
pub mod db;
pub mod models;
pub mod parser;
pub mod report;
pub mod scorer;

use models::{AnnotationConfig, ConfidenceTier, Report, Variant};
use rusqlite::Connection;

/// Errors that can occur during the full analysis pipeline.
#[derive(Debug, thiserror::Error)]
pub enum AnalyzeError {
    /// Annotation phase failed.
    #[error("annotation failed: {0}")]
    Annotate(#[from] annotator::AnnotateError),

    /// Report rendering failed.
    #[error("report generation failed: {0}")]
    Report(#[from] report::ReportError),
}

/// Medical disclaimer text included in every report.
const DISCLAIMER: &str = "\
This report is for informational and educational purposes only. It is NOT \
medical advice and should NOT be used for clinical decision-making. Genetic \
variants are interpreted using publicly available databases, which may contain \
errors or incomplete information. Many genetic findings have low predictive \
value for individual health outcomes. Always consult a qualified healthcare \
professional or certified genetic counselor before making any health decisions \
based on genetic data.";

/// Standard data source attributions.
const ATTRIBUTIONS: &[&str] = &[
    "ClinVar: National Center for Biotechnology Information (NCBI), National Library of Medicine (NLM) — public domain",
    "GWAS Catalog: NHGRI-EBI Catalog of human genome-wide association studies — open access (EBI Terms of Use)",
    "gnomAD / dbSNP: Genome Aggregation Database / NCBI dbSNP — open access",
    "PharmGKB: Pharmacogenomics Knowledge Base — Creative Commons Attribution-ShareAlike 4.0",
];

/// SNPedia attribution, added only when SNPedia data is used.
const SNPEDIA_ATTRIBUTION: &str =
    "SNPedia: SNPedia.com — Creative Commons Attribution-NonCommercial-ShareAlike 3.0";

/// Run the full analysis pipeline: annotate, score, filter, and build report.
///
/// This is the primary entry point for the GeneSight core library. It takes
/// parsed variants and database connections, runs annotation and scoring,
/// filters results to the requested confidence tiers, and returns a complete
/// `Report` struct.
///
/// # Arguments
///
/// * `variants` - Parsed variants from a DNA raw data file
/// * `main_db` - Connection to `genesight.db` (ClinVar, GWAS, frequencies, PharmGKB)
/// * `snpedia_db` - Optional connection to `snpedia.db`
/// * `tiers` - Which confidence tiers to include in the report (empty = all tiers)
///
/// # Errors
///
/// Returns `AnalyzeError::Annotate` if database queries fail, or
/// `AnalyzeError::Report` if report generation fails.
pub fn analyze(
    variants: &[Variant],
    main_db: &Connection,
    snpedia_db: Option<&Connection>,
    tiers: &[ConfidenceTier],
) -> Result<Report, AnalyzeError> {
    analyze_with_config(
        variants,
        main_db,
        snpedia_db,
        tiers,
        &AnnotationConfig::default(),
    )
}

/// Run the full analysis pipeline with selective database configuration.
///
/// Like [`analyze`], but accepts an [`AnnotationConfig`] to control which
/// databases are queried during annotation.
pub fn analyze_with_config(
    variants: &[Variant],
    main_db: &Connection,
    snpedia_db: Option<&Connection>,
    tiers: &[ConfidenceTier],
    config: &AnnotationConfig,
) -> Result<Report, AnalyzeError> {
    tracing::info!(
        total_variants = variants.len(),
        tier_filter = ?tiers,
        snpedia = snpedia_db.is_some(),
        "starting analysis pipeline"
    );

    // Step 1: Annotate variants against selected databases
    let annotated =
        annotator::annotate_variants_with_config(variants, main_db, snpedia_db, config)?;
    let annotated_count = annotated.len();

    // Step 2: Score annotated variants
    let scored = scorer::score_variants(&annotated);

    // Step 3: Filter by requested tiers (empty slice = include all)
    let filtered: Vec<_> = if tiers.is_empty() {
        scored
    } else {
        scored
            .into_iter()
            .filter(|r| tiers.contains(&r.tier))
            .collect()
    };

    tracing::info!(
        annotated = annotated_count,
        scored_total = filtered.len(),
        "pipeline complete, building report"
    );

    // Step 4: Build attributions list
    let mut attributions: Vec<String> = ATTRIBUTIONS.iter().map(|s| (*s).to_string()).collect();
    if snpedia_db.is_some() {
        attributions.push(SNPEDIA_ATTRIBUTION.to_string());
    }

    // Step 5: Build the Report
    let report = Report {
        total_variants: variants.len(),
        annotated_variants: annotated_count,
        results: filtered,
        attributions,
        disclaimer: DISCLAIMER.to_string(),
    };

    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use models::variant::{Genotype, SourceFormat};

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
            INSERT INTO pharmacogenomics VALUES (
                'rs456', 'Codeine', 'Poor Metabolizer', '1A',
                'Consider alternative analgesic', 'CYP2D6'
            );
            "#,
        )
        .expect("setup");
        conn
    }

    fn make_variant(rsid: &str) -> Variant {
        Variant {
            rsid: Some(rsid.to_string()),
            chromosome: "1".to_string(),
            position: 100000,
            genotype: Genotype::Heterozygous('A', 'G'),
            source_format: SourceFormat::TwentyThreeAndMe,
        }
    }

    #[test]
    fn analyze_full_pipeline() {
        let db = setup_main_db();
        let variants = vec![
            make_variant("rs123"),
            make_variant("rs456"),
            make_variant("rs999"),
        ];

        let report = analyze(&variants, &db, None, &[]).expect("analyze");

        assert_eq!(report.total_variants, 3);
        assert_eq!(report.annotated_variants, 2);
        assert!(!report.results.is_empty());
        assert!(!report.disclaimer.is_empty());
        assert!(!report.attributions.is_empty());
    }

    #[test]
    fn analyze_with_tier_filter() {
        let db = setup_main_db();
        let variants = vec![make_variant("rs123"), make_variant("rs456")];

        let report =
            analyze(&variants, &db, None, &[ConfidenceTier::Tier1Reliable]).expect("analyze");

        // Both rs123 (ClinVar pathogenic 3-star) and rs456 (PharmGKB 1A) are Tier1
        for result in &report.results {
            assert_eq!(result.tier, ConfidenceTier::Tier1Reliable);
        }
    }

    #[test]
    fn analyze_empty_variants() {
        let db = setup_main_db();
        let report = analyze(&[], &db, None, &[]).expect("analyze");
        assert_eq!(report.total_variants, 0);
        assert_eq!(report.annotated_variants, 0);
        assert!(report.results.is_empty());
    }

    #[test]
    fn analyze_with_snpedia_adds_attribution() {
        let db = setup_main_db();
        let snpedia = Connection::open_in_memory().expect("open snpedia");
        snpedia
            .execute_batch(
                "CREATE TABLE snpedia (
                    rsid TEXT PRIMARY KEY, magnitude REAL, repute TEXT,
                    summary TEXT, genotype_descriptions TEXT
                );",
            )
            .expect("setup snpedia");

        let report = analyze(&[], &db, Some(&snpedia), &[]).expect("analyze");
        assert!(report.attributions.iter().any(|a| a.contains("SNPedia")));
    }

    #[test]
    fn analyze_without_snpedia_omits_attribution() {
        let db = setup_main_db();
        let report = analyze(&[], &db, None, &[]).expect("analyze");
        assert!(!report.attributions.iter().any(|a| a.contains("SNPedia")));
    }

    #[test]
    fn report_renders_to_all_formats() {
        let db = setup_main_db();
        let variants = vec![make_variant("rs123")];
        let report_data = analyze(&variants, &db, None, &[]).expect("analyze");

        let md = report::render(&report_data, report::OutputFormat::Markdown).expect("md");
        assert!(md.contains("GeneSight"));

        let json = report::render(&report_data, report::OutputFormat::Json).expect("json");
        assert!(json.contains("rs123"));

        let html = report::render(&report_data, report::OutputFormat::Html).expect("html");
        assert!(html.contains("<!DOCTYPE html>"));
    }
}
