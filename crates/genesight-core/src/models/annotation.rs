use serde::{Deserialize, Serialize};

use super::variant::Variant;

/// A variant with all database annotations attached.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedVariant {
    /// The original parsed variant
    pub variant: Variant,
    /// ClinVar annotation (if found)
    pub clinvar: Option<ClinVarAnnotation>,
    /// SNPedia annotation (if found, from optional snpedia.db)
    pub snpedia: Option<SnpediaAnnotation>,
    /// GWAS Catalog hits (zero or more)
    pub gwas_hits: Vec<GwasHit>,
    /// Allele frequency data (gnomAD/dbSNP)
    pub frequency: Option<AlleleFrequency>,
    /// Pharmacogenomic annotation (if found)
    pub pharmacogenomics: Option<PharmaAnnotation>,
    /// Reference allele at this position (from the `variants` table).
    /// Used for allele matching to determine if the user carries the variant.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ref_allele: Option<String>,
    /// Alternate allele at this position (from the `variants` table).
    /// Used for allele matching to determine if the user carries the variant.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alt_allele: Option<String>,
}

/// Classification context for a ClinVar entry.
///
/// Since 2024, ClinVar separates germline, somatic, and oncogenicity
/// classifications. In a consumer DNA (germline) context, somatic and
/// oncogenicity classifications should be treated as informational only.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ClinVarClassificationType {
    /// Germline (inherited) classification — most relevant for consumer DNA.
    Germline,
    /// Somatic (tumor tissue) classification — informational in germline context.
    Somatic,
    /// Oncogenicity classification — informational in germline context.
    Oncogenicity,
}

impl ClinVarClassificationType {
    /// Parse a classification type from its database string representation.
    ///
    /// Unrecognized values default to `Germline` for backward compatibility
    /// with databases that predate the classification type column.
    pub fn from_db_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "germline" => Self::Germline,
            "somatic" => Self::Somatic,
            "oncogenicity" => Self::Oncogenicity,
            _ => Self::Germline,
        }
    }

    /// Return a human-readable label for display in reports.
    pub fn label(self) -> &'static str {
        match self {
            Self::Germline => "Germline",
            Self::Somatic => "Somatic",
            Self::Oncogenicity => "Oncogenicity",
        }
    }
}

impl std::fmt::Display for ClinVarClassificationType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.label())
    }
}

/// ClinVar clinical classification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClinVarAnnotation {
    /// Clinical significance (e.g., "Pathogenic", "Benign")
    pub significance: String,
    /// Review status (0-4 stars)
    pub review_stars: u8,
    /// Associated conditions/diseases
    pub conditions: Vec<String>,
    /// Gene symbol (e.g., "BRCA1")
    pub gene_symbol: Option<String>,
    /// Classification context (germline, somatic, or oncogenicity).
    ///
    /// Defaults to `Germline` for databases that predate the 2024
    /// ClinVar classification split.
    pub classification_type: ClinVarClassificationType,
}

/// SNPedia wiki annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpediaAnnotation {
    /// Importance magnitude (0-10, higher = more important)
    pub magnitude: f64,
    /// Repute: "good", "bad", or None
    pub repute: Option<String>,
    /// Human-readable summary
    pub summary: String,
    /// Genotype-specific descriptions (e.g., {"AA": "...", "AG": "..."})
    pub genotype_descriptions: Option<std::collections::HashMap<String, String>>,
}

/// A hit from the GWAS Catalog.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GwasHit {
    /// Associated trait or disease
    pub trait_name: String,
    /// p-value of the association
    pub p_value: f64,
    /// Odds ratio (if available)
    pub odds_ratio: Option<f64>,
    /// Beta coefficient (if available)
    pub beta: Option<f64>,
    /// Risk allele
    pub risk_allele: Option<String>,
    /// Risk allele frequency
    pub risk_allele_frequency: Option<f64>,
    /// PubMed ID of the study
    pub pubmed_id: Option<String>,
    /// Mapped gene
    pub mapped_gene: Option<String>,
}

/// Allele frequency from gnomAD or dbSNP.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleFrequency {
    /// Overall allele frequency
    pub af_total: f64,
    /// African population
    pub af_afr: Option<f64>,
    /// American population
    pub af_amr: Option<f64>,
    /// East Asian population
    pub af_eas: Option<f64>,
    /// European (non-Finnish) population
    pub af_eur: Option<f64>,
    /// South Asian population
    pub af_sas: Option<f64>,
    /// Data source ("gnomad" or "dbsnp")
    pub source: String,
}

/// Pharmacogenomic annotation from PharmGKB.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmaAnnotation {
    /// Affected gene (e.g., "CYP2D6")
    pub gene: String,
    /// Affected drug
    pub drug: String,
    /// Metabolizer phenotype (e.g., "Poor Metabolizer")
    pub phenotype_category: Option<String>,
    /// PharmGKB evidence level (1A, 1B, 2A, 2B, 3, 4)
    pub evidence_level: String,
    /// Clinical recommendation
    pub clinical_recommendation: Option<String>,
}
