use serde::{Deserialize, Serialize};

/// Confidence tier for analysis results.
///
/// Every result MUST be assigned a confidence tier. The report displays
/// this prominently to help users understand the reliability of each finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum ConfidenceTier {
    /// Tier 1: Reliable (>95% predictive value)
    ///
    /// Monogenic diseases, carrier status, pharmacogenetics, simple traits.
    /// Sources: ClinVar (review status >= 2 stars), PharmGKB (Level 1A/1B).
    Tier1Reliable,

    /// Tier 2: Probable (60-85% predictive value)
    ///
    /// Polygenic risk scores, physical traits.
    /// Sources: GWAS Catalog (genome-wide significant, moderate effect),
    /// SNPedia (magnitude 2-3.9).
    Tier2Probable,

    /// Tier 3: Speculative (50-65% predictive value)
    ///
    /// Complex diseases, personality traits, athletic aptitude.
    /// Sources: GWAS Catalog (low effect size), SNPedia (magnitude < 2).
    Tier3Speculative,
}

impl std::fmt::Display for ConfidenceTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfidenceTier::Tier1Reliable => write!(f, "Tier 1: Reliable (>95%)"),
            ConfidenceTier::Tier2Probable => write!(f, "Tier 2: Probable (60-85%)"),
            ConfidenceTier::Tier3Speculative => write!(f, "Tier 3: Speculative (50-65%)"),
        }
    }
}
