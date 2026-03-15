//! Risk scoring and confidence tier assignment.
//!
//! Every finding MUST be assigned a [`ConfidenceTier`] based on the
//! quality and type of evidence. Variants with only frequency data
//! (no clinical, GWAS, pharmacogenomic, or SNPedia annotations) are
//! not scored and are excluded from results.

pub mod monogenic;
pub mod pharmaco;
pub mod polygenic;
pub mod traits;

use crate::models::annotation::{AnnotatedVariant, GwasHit};
use crate::models::confidence::ConfidenceTier;
use crate::models::report::{ResultCategory, ScoredResult};

/// Score annotated variants and assign confidence tiers.
///
/// Evaluates each annotated variant against scoring rules for different
/// evidence types (ClinVar, PharmGKB, GWAS, SNPedia) and produces
/// `ScoredResult` entries with tier, category, and human-readable summaries.
///
/// Variants with only allele frequency data are skipped since frequency
/// alone is not a clinical finding.
///
/// A single annotated variant may produce multiple scored results if it
/// has annotations from multiple databases (e.g., both ClinVar and GWAS).
pub fn score_variants(annotated: &[AnnotatedVariant]) -> Vec<ScoredResult> {
    let mut results = Vec::new();

    for av in annotated {
        let rsid = av.variant.rsid.as_deref().unwrap_or("unknown");
        let genotype = av.variant.genotype.to_string();

        // Score ClinVar annotations
        if let Some(clinvar) = &av.clinvar {
            if let Some(scored) = score_clinvar(av, rsid, &genotype, clinvar) {
                results.push(scored);
            }
        }

        // Score pharmacogenomic annotations
        if let Some(pharma) = &av.pharmacogenomics {
            if let Some(scored) = score_pharma(av, rsid, &genotype, pharma) {
                results.push(scored);
            }
        }

        // Score GWAS hits (each hit scored independently)
        for hit in &av.gwas_hits {
            if let Some(scored) = score_gwas_hit(av, rsid, &genotype, hit) {
                results.push(scored);
            }
        }

        // Score SNPedia annotations
        if let Some(snpedia) = &av.snpedia {
            if let Some(scored) = score_snpedia(av, rsid, &genotype, snpedia) {
                results.push(scored);
            }
        }
    }

    // Sort by tier (Tier1 first) then by category
    results.sort_by(|a, b| a.tier.cmp(&b.tier).then(a.category.cmp(&b.category)));

    tracing::info!(scored = results.len(), "scoring complete");
    results
}

/// Score a ClinVar annotation.
///
/// - review_stars >= 2 + pathogenic/likely pathogenic => Tier1, MonogenicDisease
/// - significance containing "carrier" or benign => CarrierStatus (Tier1 if stars >= 2, else Tier2)
/// - Other significance values with low stars => Tier2
fn score_clinvar(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    clinvar: &crate::models::annotation::ClinVarAnnotation,
) -> Option<ScoredResult> {
    let sig_lower = clinvar.significance.to_lowercase();
    let gene = clinvar.gene_symbol.as_deref().unwrap_or("unknown gene");
    let conditions_str = if clinvar.conditions.is_empty() {
        "unspecified condition".to_string()
    } else {
        clinvar.conditions.join(", ")
    };

    let is_conflicting = sig_lower.contains("conflicting");
    let is_uncertain = sig_lower.contains("uncertain significance");
    let is_pathogenic = sig_lower.contains("pathogenic")
        && !sig_lower.contains("benign")
        && !is_conflicting
        && !is_uncertain;
    let is_carrier_or_benign = sig_lower.contains("carrier") || sig_lower.contains("benign");

    // Skip conflicting and uncertain classifications — not actionable
    if is_conflicting || is_uncertain {
        return None;
    }

    if is_pathogenic && clinvar.review_stars >= 2 {
        Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier1Reliable,
            category: ResultCategory::MonogenicDisease,
            summary: format!(
                "{gene} ({rsid}) — {sig} ({stars}-star review)",
                sig = clinvar.significance,
                stars = clinvar.review_stars,
            ),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star. This variant has strong clinical evidence.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
        })
    } else if is_carrier_or_benign {
        let tier = if clinvar.review_stars >= 2 {
            ConfidenceTier::Tier1Reliable
        } else {
            ConfidenceTier::Tier2Probable
        };
        Some(ScoredResult {
            variant: av.clone(),
            tier,
            category: ResultCategory::CarrierStatus,
            summary: format!("{gene} ({rsid}) — {sig}", sig = clinvar.significance,),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
        })
    } else if is_pathogenic {
        // Pathogenic but low review stars
        Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier2Probable,
            category: ResultCategory::MonogenicDisease,
            summary: format!(
                "{gene} ({rsid}) — {sig} ({stars}-star review, limited evidence)",
                sig = clinvar.significance,
                stars = clinvar.review_stars,
            ),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star. Lower review status indicates limited evidence.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
        })
    } else {
        // Other significance values (VUS, conflicting, etc.) — skip
        None
    }
}

/// Score a PharmGKB annotation.
///
/// - evidence_level 1A or 1B => Tier1, Pharmacogenomics
/// - evidence_level 2A or 2B => Tier2, Pharmacogenomics
/// - Lower evidence levels => Tier3
fn score_pharma(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    pharma: &crate::models::annotation::PharmaAnnotation,
) -> Option<ScoredResult> {
    let level = pharma.evidence_level.trim();

    let tier = match level {
        "1A" | "1B" => ConfidenceTier::Tier1Reliable,
        "2A" | "2B" => ConfidenceTier::Tier2Probable,
        _ => ConfidenceTier::Tier3Speculative,
    };

    let phenotype = pharma
        .phenotype_category
        .as_deref()
        .unwrap_or("Unknown phenotype");

    let recommendation = pharma
        .clinical_recommendation
        .as_deref()
        .unwrap_or("No specific recommendation available");

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category: ResultCategory::Pharmacogenomics,
        summary: format!(
            "{gene} ({rsid}) — {drug}: {phenotype}",
            gene = pharma.gene,
            drug = pharma.drug,
        ),
        details: format!(
            "Genotype: {genotype}. Gene: {}. Drug: {}. Phenotype: {phenotype}. \
             Evidence level: {level}. Recommendation: {recommendation}.",
            pharma.gene, pharma.drug,
        ),
    })
}

/// Score a single GWAS hit.
///
/// - p_value < 5e-8 and odds_ratio > 1.5 => Tier2, PolygenicRiskScore
/// - Otherwise => Tier3, ComplexTrait
fn score_gwas_hit(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    hit: &GwasHit,
) -> Option<ScoredResult> {
    let genome_wide_significant = hit.p_value < 5e-8;
    let moderate_effect = hit.odds_ratio.is_some_and(|or| or > 1.5);
    let gene = hit.mapped_gene.as_deref().unwrap_or("intergenic");

    let (tier, category) = if genome_wide_significant && moderate_effect {
        (
            ConfidenceTier::Tier2Probable,
            ResultCategory::PolygenicRiskScore,
        )
    } else {
        (
            ConfidenceTier::Tier3Speculative,
            ResultCategory::ComplexTrait,
        )
    };

    let effect_desc = if let Some(or) = hit.odds_ratio {
        format!("odds ratio {or:.2}")
    } else if let Some(beta) = hit.beta {
        format!("beta {beta:.3}")
    } else {
        "effect size not reported".to_string()
    };

    let pubmed = hit
        .pubmed_id
        .as_deref()
        .map(|id| format!(" (PMID: {id})"))
        .unwrap_or_default();

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category,
        summary: format!(
            "{gene} ({rsid}) — {trait_name}: {effect_desc}",
            trait_name = hit.trait_name,
        ),
        details: format!(
            "Genotype: {genotype}. Trait: {}. p-value: {:.2e}, {effect_desc}. \
             Mapped gene: {gene}.{pubmed}",
            hit.trait_name, hit.p_value,
        ),
    })
}

/// Score a SNPedia annotation.
///
/// - magnitude >= 3 => Tier2, PhysicalTrait
/// - magnitude < 3 => Tier3, ComplexTrait
fn score_snpedia(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    snpedia: &crate::models::annotation::SnpediaAnnotation,
) -> Option<ScoredResult> {
    // Skip very low magnitude entries (essentially noise)
    if snpedia.magnitude < 0.5 {
        return None;
    }

    let (tier, category) = if snpedia.magnitude >= 3.0 {
        (ConfidenceTier::Tier2Probable, ResultCategory::PhysicalTrait)
    } else {
        (
            ConfidenceTier::Tier3Speculative,
            ResultCategory::ComplexTrait,
        )
    };

    let repute_desc = match snpedia.repute.as_deref() {
        Some("good") => " (positive)",
        Some("bad") => " (negative)",
        _ => "",
    };

    // Look up genotype-specific description if available
    let geno_desc = snpedia
        .genotype_descriptions
        .as_ref()
        .and_then(|descs| descs.get(genotype))
        .map(|d| format!(" Your genotype ({genotype}): {d}."))
        .unwrap_or_default();

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category,
        summary: format!(
            "{rsid} — {summary}{repute_desc} (magnitude {mag:.1})",
            summary = snpedia.summary,
            mag = snpedia.magnitude,
        ),
        details: format!(
            "Genotype: {genotype}. SNPedia magnitude: {mag:.1}/10. {summary}.{geno_desc} \
             Source: SNPedia (CC-BY-NC-SA 3.0).",
            mag = snpedia.magnitude,
            summary = snpedia.summary,
        ),
    })
}

/// Display implementation for `ResultCategory` used in reports.
impl std::fmt::Display for ResultCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ResultCategory::MonogenicDisease => write!(f, "Monogenic Disease Risk"),
            ResultCategory::CarrierStatus => write!(f, "Carrier Status"),
            ResultCategory::Pharmacogenomics => write!(f, "Pharmacogenomics"),
            ResultCategory::PolygenicRiskScore => write!(f, "Polygenic Risk Score"),
            ResultCategory::PhysicalTrait => write!(f, "Physical Traits"),
            ResultCategory::ComplexTrait => write!(f, "Complex Traits"),
            ResultCategory::Ancestry => write!(f, "Ancestry"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::annotation::*;
    use crate::models::variant::{Genotype, SourceFormat, Variant};

    fn make_variant(rsid: &str) -> Variant {
        Variant {
            rsid: Some(rsid.to_string()),
            chromosome: "1".to_string(),
            position: 100000,
            genotype: Genotype::Heterozygous('A', 'G'),
            source_format: SourceFormat::TwentyThreeAndMe,
        }
    }

    fn make_annotated(rsid: &str) -> AnnotatedVariant {
        AnnotatedVariant {
            variant: make_variant(rsid),
            clinvar: None,
            snpedia: None,
            gwas_hits: Vec::new(),
            frequency: None,
            pharmacogenomics: None,
        }
    }

    #[test]
    fn clinvar_pathogenic_high_stars_is_tier1() {
        let mut av = make_annotated("rs123");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Breast cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
    }

    #[test]
    fn clinvar_pathogenic_low_stars_is_tier2() {
        let mut av = make_annotated("rs123");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Likely pathogenic".to_string(),
            review_stars: 1,
            conditions: vec!["Some disease".to_string()],
            gene_symbol: Some("GENE1".to_string()),
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
    }

    #[test]
    fn clinvar_benign_is_carrier_status() {
        let mut av = make_annotated("rs456");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Benign".to_string(),
            review_stars: 2,
            conditions: vec![],
            gene_symbol: Some("TP53".to_string()),
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].category, ResultCategory::CarrierStatus);
    }

    #[test]
    fn pharmgkb_1a_is_tier1() {
        let mut av = make_annotated("rs1065852");
        av.pharmacogenomics = Some(PharmaAnnotation {
            gene: "CYP2D6".to_string(),
            drug: "Codeine".to_string(),
            phenotype_category: Some("Poor Metabolizer".to_string()),
            evidence_level: "1A".to_string(),
            clinical_recommendation: Some("Consider alternative".to_string()),
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert_eq!(results[0].category, ResultCategory::Pharmacogenomics);
    }

    #[test]
    fn pharmgkb_2a_is_tier2() {
        let mut av = make_annotated("rs999");
        av.pharmacogenomics = Some(PharmaAnnotation {
            gene: "DPYD".to_string(),
            drug: "Fluorouracil".to_string(),
            phenotype_category: None,
            evidence_level: "2A".to_string(),
            clinical_recommendation: None,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
    }

    #[test]
    fn gwas_significant_strong_effect_is_tier2_polygenic() {
        let mut av = make_annotated("rs100");
        av.gwas_hits = vec![GwasHit {
            trait_name: "Type 2 Diabetes".to_string(),
            p_value: 1e-12,
            odds_ratio: Some(2.0),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: Some("12345".to_string()),
            mapped_gene: Some("TCF7L2".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
        assert_eq!(results[0].category, ResultCategory::PolygenicRiskScore);
    }

    #[test]
    fn gwas_weak_effect_is_tier3_complex() {
        let mut av = make_annotated("rs200");
        av.gwas_hits = vec![GwasHit {
            trait_name: "Height".to_string(),
            p_value: 1e-9,
            odds_ratio: Some(1.1),
            beta: None,
            risk_allele: None,
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: Some("HMGA2".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert_eq!(results[0].category, ResultCategory::ComplexTrait);
    }

    #[test]
    fn snpedia_high_magnitude_is_tier2() {
        let mut av = make_annotated("rs300");
        av.snpedia = Some(SnpediaAnnotation {
            magnitude: 4.0,
            repute: Some("bad".to_string()),
            summary: "Important finding".to_string(),
            genotype_descriptions: None,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
        assert_eq!(results[0].category, ResultCategory::PhysicalTrait);
    }

    #[test]
    fn snpedia_low_magnitude_is_tier3() {
        let mut av = make_annotated("rs400");
        av.snpedia = Some(SnpediaAnnotation {
            magnitude: 2.0,
            repute: None,
            summary: "Minor trait association".to_string(),
            genotype_descriptions: None,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert_eq!(results[0].category, ResultCategory::ComplexTrait);
    }

    #[test]
    fn snpedia_very_low_magnitude_is_skipped() {
        let mut av = make_annotated("rs500");
        av.snpedia = Some(SnpediaAnnotation {
            magnitude: 0.1,
            repute: None,
            summary: "Noise".to_string(),
            genotype_descriptions: None,
        });

        let results = score_variants(&[av]);
        assert!(results.is_empty());
    }

    #[test]
    fn frequency_only_is_not_scored() {
        let mut av = make_annotated("rs600");
        av.frequency = Some(AlleleFrequency {
            af_total: 0.25,
            af_afr: None,
            af_amr: None,
            af_eas: None,
            af_eur: None,
            af_sas: None,
            source: "gnomad".to_string(),
        });

        let results = score_variants(&[av]);
        assert!(results.is_empty());
    }

    #[test]
    fn multiple_annotations_produce_multiple_scores() {
        let mut av = make_annotated("rs700");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Disease X".to_string()],
            gene_symbol: Some("GENEX".to_string()),
        });
        av.gwas_hits = vec![GwasHit {
            trait_name: "Disease X risk".to_string(),
            p_value: 1e-15,
            odds_ratio: Some(3.0),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.1),
            pubmed_id: None,
            mapped_gene: Some("GENEX".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn results_are_sorted_by_tier() {
        let mut av1 = make_annotated("rs10");
        av1.snpedia = Some(SnpediaAnnotation {
            magnitude: 1.5,
            repute: None,
            summary: "Low".to_string(),
            genotype_descriptions: None,
        });

        let mut av2 = make_annotated("rs20");
        av2.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Serious".to_string()],
            gene_symbol: Some("GENE".to_string()),
        });

        let results = score_variants(&[av1, av2]);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert_eq!(results[1].tier, ConfidenceTier::Tier3Speculative);
    }
}
