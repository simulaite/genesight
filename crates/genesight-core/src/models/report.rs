use serde::{Deserialize, Serialize};

use super::{AnnotatedVariant, ConfidenceTier, GenomeAssembly};

/// A scored result with confidence tier and human-readable summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoredResult {
    /// The annotated variant
    pub variant: AnnotatedVariant,
    /// Assigned confidence tier
    pub tier: ConfidenceTier,
    /// Category of this result
    pub category: ResultCategory,
    /// Human-readable summary of the finding
    pub summary: String,
    /// More detailed explanation
    pub details: String,
    /// Caveats or limitations that affect interpretation of this result.
    ///
    /// Examples: palindromic SNP strand ambiguity, low-confidence allele
    /// frequency data, conflicting database entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub limitations: Vec<String>,
}

/// Categories for organizing results in the report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ResultCategory {
    /// Monogenic disease risk (e.g., BRCA1/2, CFTR)
    MonogenicDisease,
    /// Carrier status for recessive conditions
    CarrierStatus,
    /// Drug metabolism and interactions
    Pharmacogenomics,
    /// Polygenic risk scores (diabetes, heart disease)
    PolygenicRiskScore,
    /// Physical traits (hair color, lactose tolerance)
    PhysicalTrait,
    /// Complex traits (speculative)
    ComplexTrait,
    /// Ancestry markers
    Ancestry,
}

/// A complete analysis report for a DNA sample.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    /// Total number of variants in the input file
    pub total_variants: usize,
    /// Number of variants with at least one annotation
    pub annotated_variants: usize,
    /// Scored results grouped by tier
    pub results: Vec<ScoredResult>,
    /// Data source attributions (required by licenses)
    pub attributions: Vec<String>,
    /// Medical disclaimer (mandatory)
    pub disclaimer: String,
    /// Genome assembly detected from the input file
    pub input_assembly: GenomeAssembly,
    /// Genome assembly of the reference database
    pub db_assembly: GenomeAssembly,
    /// Warnings about assembly mismatches or unknown assemblies
    pub assembly_warnings: Vec<String>,
}
