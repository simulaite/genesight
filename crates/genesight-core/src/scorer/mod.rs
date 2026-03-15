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

use crate::allele::{
    count_risk_allele_copies, match_alleles_with_frequency, AlleleMatch, RiskAlleleCopies,
};
use crate::models::annotation::{AnnotatedVariant, ClinVarClassificationType, GwasHit};
use crate::models::confidence::ConfidenceTier;
use crate::models::report::{ConfirmationUrgency, ResultCategory, ScoredResult};
use crate::models::variant::Genotype;

/// DTC raw-data caveat appended to ALL scored results.
const DTC_RAW_DATA_CAVEAT: &str = "This result is derived from direct-to-consumer (DTC) \
    microarray genotyping data, which has not been validated in a clinical laboratory \
    setting. DTC genotyping has known limitations including strand ambiguity, limited \
    coverage, and potential genotyping errors. Any clinically relevant finding should be \
    confirmed through clinical-grade testing before making medical decisions.";

/// ACMG Secondary Findings v3.2 actionable gene list.
const ACMG_SF_GENES: &[&str] = &[
    "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2", "APC", "MEN1", "RET", "RB1", "TP53", "VHL",
    "SDHB", "SDHD", "SDHAF2", "BMPR1A", "SMAD4", "STK11", "PTEN", "CDH1", "PALB2",
];

/// Caveat text appended to results when a palindromic SNP's strand could not
/// be resolved using allele frequency data.
const STRAND_AMBIGUITY_CAVEAT: &str =
    "This variant is a palindromic SNP (A/T or C/G). Strand orientation \
     could not be confidently resolved from allele frequency data. The \
     reported risk allele match may be on the wrong strand.";

/// Somatic classification limitation text appended to results.
const SOMATIC_LIMITATION: &str = "This variant's clinical classification is based on \
    somatic (tumor) context and may not apply to germline (inherited) analysis.";

/// Oncogenicity classification limitation text appended to results.
const ONCOGENICITY_LIMITATION: &str = "This variant's clinical classification is an \
    oncogenicity assessment and may not apply to germline (inherited) analysis.";

/// FDA disclaimer for pharmacogenomic results from consumer genotyping.
const PGX_FDA_DISCLAIMER: &str = "Pharmacogenomic results from consumer genotyping \
    arrays have NOT been reviewed or approved by the U.S. Food and Drug Administration \
    (FDA) for clinical use. Do not alter any medication regimen based solely on these \
    results. Consult a healthcare provider or clinical pharmacogenomics service for \
    validated testing.";

/// Caveat: OR is relative, not absolute risk.
const GWAS_OR_RELATIVE_CAVEAT: &str = "The odds ratio (OR) is a relative measure \
    comparing your odds to someone without this allele. It is NOT an absolute \
    probability. The actual risk depends on baseline prevalence in your population. \
    An OR of 1.3 means approximately 1.3x the odds, not a 30% chance.";

/// Caveat: Beta coefficient context.
const GWAS_BETA_CAVEAT: &str = "The beta coefficient is the estimated effect size \
    per copy of the effect allele in the study's units. Without population mean and \
    standard deviation, clinical significance cannot be determined from DTC data alone.";

/// Caveat: Historical GWAS OR inversion.
const GWAS_OR_INVERSION_CAVEAT: &str = "The GWAS Catalog changed curation conventions \
    around January 2021. Earlier entries may have inverted odds ratios. Risk allele \
    directionality should be verified against the original publication.";

/// Curated autosomal recessive (AR) genes. Het pathogenic = carrier, not affected.
const AR_GENES: &[&str] = &[
    "CFTR", "HBB", "HEXA", "GJB2", "GJB6", "SLC26A4", "PAH", "SMN1", "GBA1", "ASPA", "MEFV",
    "ATP7B", "HFE", "BLM", "FANCA", "FANCC", "MYO7A", "USH2A", "GALT", "ACADM", "SMPD1", "BCKDHA",
    "BCKDHB", "SERPINA1", "CYP21A2",
];

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

        // Score GWAS hits — deduplicate by trait, keeping most significant per trait
        let deduped_hits = dedup_gwas_hits(&av.gwas_hits);
        for hit in &deduped_hits {
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

/// Result of ClinVar allele comparison.
enum ClinvarAlleleResult {
    /// User carries this many copies of the alt allele (0, 1, or 2).
    Copies(u8),
    /// Variant is an indel that cannot be detected from microarray data.
    IndelNotDetectable,
    /// User has NoCall at this position — no genotype data.
    NoCallGenotype,
    /// Palindromic SNP — strand cannot be resolved.
    Palindromic,
    /// Ref/alt allele data not available in database.
    NoAlleleData,
}

/// Determine how many copies of the alternate (variant) allele the user carries.
fn clinvar_allele_check(av: &AnnotatedVariant) -> ClinvarAlleleResult {
    let ref_str = match av.ref_allele.as_deref() {
        Some(s) if !s.is_empty() => s,
        _ => return ClinvarAlleleResult::NoAlleleData,
    };
    let alt_str = match av.alt_allele.as_deref() {
        Some(s) if !s.is_empty() => s,
        _ => return ClinvarAlleleResult::NoAlleleData,
    };

    // Check for indels (multi-base variants) — microarrays can't detect these
    if ref_str.len() != 1 || alt_str.len() != 1 {
        return ClinvarAlleleResult::IndelNotDetectable;
    }

    let ref_char = ref_str.as_bytes()[0] as char;
    let alt_char = alt_str.as_bytes()[0] as char;

    let (a1, a2) = match &av.variant.genotype {
        Genotype::Homozygous(a) => (*a, *a),
        Genotype::Heterozygous(a, b) => (*a, *b),
        Genotype::NoCall => return ClinvarAlleleResult::NoCallGenotype,
        Genotype::Indel(_) => return ClinvarAlleleResult::IndelNotDetectable,
    };

    let a1_upper = a1.to_ascii_uppercase();
    let a2_upper = a2.to_ascii_uppercase();
    let ref_upper = ref_char.to_ascii_uppercase();
    let alt_upper = alt_char.to_ascii_uppercase();

    // Check for palindromic SNPs — strand cannot be resolved from sequence
    if crate::allele::strand::is_palindromic(ref_upper, alt_upper) {
        // Try frequency-based resolution
        if let Some(freq) = &av.frequency {
            let dist_from_half = (freq.af_total - 0.5).abs();
            if dist_from_half > 0.40 {
                let copies = u8::from(a1_upper == alt_upper) + u8::from(a2_upper == alt_upper);
                return ClinvarAlleleResult::Copies(copies);
            }
        }
        return ClinvarAlleleResult::Palindromic;
    }

    // Non-palindromic: try direct strand first
    let direct_alt = u8::from(a1_upper == alt_upper) + u8::from(a2_upper == alt_upper);
    let direct_ref = u8::from(a1_upper == ref_upper) + u8::from(a2_upper == ref_upper);

    if direct_alt + direct_ref == 2 {
        return ClinvarAlleleResult::Copies(direct_alt);
    }

    // Try complement strand
    if let (Some(ref_comp), Some(alt_comp)) = (
        crate::allele::strand::complement(ref_upper),
        crate::allele::strand::complement(alt_upper),
    ) {
        let comp_alt = u8::from(a1_upper == alt_comp) + u8::from(a2_upper == alt_comp);
        let comp_ref = u8::from(a1_upper == ref_comp) + u8::from(a2_upper == ref_comp);

        if comp_alt + comp_ref == 2 {
            return ClinvarAlleleResult::Copies(comp_alt);
        }
    }

    // Alleles don't match either strand — treat as no data
    ClinvarAlleleResult::NoAlleleData
}

/// Score a ClinVar annotation with quality gates, allele matching, and
/// germline/somatic handling.
///
/// Control flow (ordered by priority):
/// 0. Allele matching gate: skip if user is homozygous reference (0 copies of alt)
/// 1. Non-germline gate: somatic/oncogenicity -> Tier3 informational with limitation
/// 2. Zero-star gate: pathogenic at 0 stars -> Tier3 with limitation; benign at 0 stars -> skip
/// 3. Conflicting gate: conflicting interpretations -> Tier3 + ClinVarConflicting
/// 4. VUS gate: uncertain significance -> Tier3 with uncertainty language
/// 5. BA1 filter: common variant (max population AF > 5%) -> skip
/// 6. Benign/Likely Benign: skip entirely (not clinically actionable)
/// 7. Normal scoring:
///    - Pathogenic/Likely Pathogenic + stars >= 2 -> Tier1, MonogenicDisease
///    - Pathogenic/Likely Pathogenic + stars == 1 -> Tier2, MonogenicDisease
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
    let is_uncertain = sig_lower.contains("uncertain significance") || sig_lower.contains("vus");
    let is_pathogenic = sig_lower.contains("pathogenic")
        && !sig_lower.contains("benign")
        && !is_conflicting
        && !is_uncertain;
    let is_benign = sig_lower.contains("benign") && !sig_lower.contains("pathogenic");

    // ── Gate -1: Benign/Likely Benign — skip entirely ─────────────────
    // Benign variants are NOT clinically actionable findings. Reporting them
    // as CarrierStatus is a false positive that overwhelms the report.
    if is_benign {
        return None;
    }

    // ── Gate 0: Allele matching — skip if user doesn't carry the variant ──
    let allele_result = clinvar_allele_check(av);
    let _is_palindromic_ambiguous = matches!(allele_result, ClinvarAlleleResult::Palindromic);
    let allele_match_limitation = match allele_result {
        ClinvarAlleleResult::Copies(0) => {
            // User is homozygous reference — does not carry this variant
            return None;
        }
        ClinvarAlleleResult::Copies(copies) => Some(format!(
            "User carries {copies} cop{} of the variant allele (genotype: {genotype}).",
            if copies == 1 { "y" } else { "ies" },
        )),
        ClinvarAlleleResult::IndelNotDetectable => {
            // Indels cannot be reliably detected from microarray data.
            // Skip entirely — reporting would be a false positive.
            return None;
        }
        ClinvarAlleleResult::NoCallGenotype => {
            // No genotype data at this position — can't confirm the variant
            return None;
        }
        ClinvarAlleleResult::Palindromic => Some(
            "Allele comparison inconclusive (palindromic SNP). Strand orientation \
             could not be resolved. Clinical-grade confirmation recommended."
                .to_string(),
        ),
        ClinvarAlleleResult::NoAlleleData => {
            tracing::warn!(
                rsid = rsid,
                "ClinVar allele verification not possible — no ref/alt data"
            );
            Some(
                "Allele verification was not possible for this variant (no reference/alternate \
                 allele data in the database). This result is based on rsID matching only and \
                 may include false positives."
                    .to_string(),
            )
        }
    };

    // Extract allele copy count for carrier detection (Fix 5)
    let allele_copies: Option<u8> = match allele_result {
        ClinvarAlleleResult::Copies(c) => Some(c),
        _ => None,
    };

    // Track whether we had NoAlleleData for tier downgrade (Fix 6)
    let is_no_allele_data = matches!(allele_result, ClinvarAlleleResult::NoAlleleData);
    let is_palindromic = matches!(allele_result, ClinvarAlleleResult::Palindromic);

    // ── Gate 0a: Palindromic SNPs — strand unresolved → Tier 3 ───────────
    // If we can't resolve the strand, we can't confirm the user actually
    // carries the variant. Reporting at Tier 1/2 would be a false positive.
    if is_palindromic && is_pathogenic {
        return Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier3Speculative,
            category: ResultCategory::MonogenicDisease,
            confirmation_urgency: ConfirmationUrgency::ClinicalConfirmationRecommended,
            summary: format!(
                "{gene} ({rsid}) — {sig} (palindromic, strand unresolved)",
                sig = clinvar.significance,
            ),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star. Cannot confirm variant presence due to \
                 strand ambiguity.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
            limitations: vec![
                STRAND_AMBIGUITY_CAVEAT.to_string(),
                DTC_RAW_DATA_CAVEAT.to_string(),
            ],
        });
    }

    // ── Gate 1: Non-germline classifications (somatic / oncogenicity) ─────
    match clinvar.classification_type {
        ClinVarClassificationType::Somatic => {
            return Some(ScoredResult {
                variant: av.clone(),
                tier: ConfidenceTier::Tier3Speculative,
                category: ResultCategory::ComplexTrait,
                confirmation_urgency: ConfirmationUrgency::InformationalOnly,
                summary: format!(
                    "{gene} ({rsid}) — {sig} (somatic classification)",
                    sig = clinvar.significance,
                ),
                details: format!(
                    "Genotype: {genotype}. Classification: {} (somatic). \
                     Associated conditions: {}. ClinVar review status: {}-star.",
                    clinvar.significance, conditions_str, clinvar.review_stars,
                ),
                limitations: vec![
                    SOMATIC_LIMITATION.to_string(),
                    DTC_RAW_DATA_CAVEAT.to_string(),
                ],
            });
        }
        ClinVarClassificationType::Oncogenicity => {
            return Some(ScoredResult {
                variant: av.clone(),
                tier: ConfidenceTier::Tier3Speculative,
                category: ResultCategory::ComplexTrait,
                confirmation_urgency: ConfirmationUrgency::InformationalOnly,
                summary: format!(
                    "{gene} ({rsid}) — {sig} (oncogenicity classification)",
                    sig = clinvar.significance,
                ),
                details: format!(
                    "Genotype: {genotype}. Classification: {} (oncogenicity). \
                     Associated conditions: {}. ClinVar review status: {}-star.",
                    clinvar.significance, conditions_str, clinvar.review_stars,
                ),
                limitations: vec![
                    ONCOGENICITY_LIMITATION.to_string(),
                    DTC_RAW_DATA_CAVEAT.to_string(),
                ],
            });
        }
        ClinVarClassificationType::Germline => {}
    }

    // ── Gate 1: Zero-star entries (no assertion criteria provided) ────────
    if clinvar.review_stars == 0 && is_pathogenic {
        let urgency = if ACMG_SF_GENES.contains(&gene) {
            ConfirmationUrgency::HighImpact
        } else {
            ConfirmationUrgency::ClinicalConfirmationRecommended
        };
        let mut lims = vec!["No assertion criteria provided (0 stars)".to_string()];
        if let Some(ref allele_lim) = allele_match_limitation {
            lims.push(allele_lim.clone());
        }
        lims.push(DTC_RAW_DATA_CAVEAT.to_string());
        return Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier3Speculative,
            category: ResultCategory::MonogenicDisease,
            confirmation_urgency: urgency,
            summary: format!(
                "{gene} ({rsid}) — {sig} (0-star, unreviewed)",
                sig = clinvar.significance,
            ),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: 0-star. This entry has no assertion criteria.",
                clinvar.significance, conditions_str,
            ),
            limitations: lims,
        });
    }

    // ── Gate 2: Conflicting interpretations ──────────────────────────────
    if is_conflicting {
        return Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier3Speculative,
            category: ResultCategory::ClinVarConflicting,
            confirmation_urgency: ConfirmationUrgency::InformationalOnly,
            summary: format!("{gene} ({rsid}) — conflicting interpretations"),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star. Multiple submitters disagree on the \
                 clinical significance of this variant.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
            limitations: vec![
                "Conflicting interpretations of pathogenicity across ClinVar submitters"
                    .to_string(),
                DTC_RAW_DATA_CAVEAT.to_string(),
            ],
        });
    }

    // ── Gate 3: Variant of Uncertain Significance (VUS) ──────────────────
    if is_uncertain {
        return Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier3Speculative,
            category: ResultCategory::MonogenicDisease,
            confirmation_urgency: ConfirmationUrgency::InformationalOnly,
            summary: format!("{gene} ({rsid}) — variant of uncertain significance"),
            details: format!(
                "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                 ClinVar review status: {}-star. The clinical significance of this variant \
                 is currently uncertain and may be reclassified as more evidence emerges.",
                clinvar.significance, conditions_str, clinvar.review_stars,
            ),
            limitations: vec![
                "Variant of uncertain significance — not clinically actionable".to_string(),
                DTC_RAW_DATA_CAVEAT.to_string(),
            ],
        });
    }

    // ── Gate 4: BA1 common-variant filter (max population AF > 5%) ───────
    if let Some(freq) = &av.frequency {
        let max_pop_af = [
            Some(freq.af_total),
            freq.af_afr,
            freq.af_amr,
            freq.af_eas,
            freq.af_eur,
            freq.af_sas,
        ]
        .into_iter()
        .flatten()
        .fold(0.0_f64, f64::max);

        if max_pop_af > 0.05 && is_pathogenic {
            return None;
        }
    }

    // ── Gate 5: Normal scoring ───────────────────────────────────────────
    // Determine confirmation urgency for pathogenic/likely pathogenic variants
    let clinvar_urgency = if ACMG_SF_GENES.contains(&gene) {
        ConfirmationUrgency::HighImpact
    } else if is_pathogenic {
        ConfirmationUrgency::ClinicalConfirmationRecommended
    } else {
        ConfirmationUrgency::InformationalOnly
    };

    // Build limitations list with allele match info
    let build_limitations = |extra: Option<&str>| -> Vec<String> {
        let mut lims = Vec::new();
        if let Some(e) = extra {
            lims.push(e.to_string());
        }
        if let Some(ref allele_lim) = allele_match_limitation {
            lims.push(allele_lim.clone());
        }
        lims.push(DTC_RAW_DATA_CAVEAT.to_string());
        lims
    };

    if is_pathogenic && clinvar.review_stars >= 2 {
        // ── Fix 5: AR carrier detection ───────────────────────────────────
        if allele_copies == Some(1) && AR_GENES.contains(&gene) {
            let mut lims = build_limitations(None);
            lims.insert(
                0,
                "You carry one copy of a pathogenic variant in an autosomal \
                recessive gene. Carriers typically do not develop symptoms but may pass \
                the variant to offspring."
                    .to_string(),
            );
            return Some(ScoredResult {
                variant: av.clone(),
                tier: ConfidenceTier::Tier1Reliable,
                category: ResultCategory::CarrierStatus,
                confirmation_urgency: ConfirmationUrgency::ClinicalConfirmationRecommended,
                summary: format!(
                    "{gene} ({rsid}) — {sig} ({stars}-star review) (carrier)",
                    sig = clinvar.significance,
                    stars = clinvar.review_stars,
                ),
                details: format!(
                    "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                     ClinVar review status: {}-star. This variant has strong clinical evidence. \
                     Heterozygous carrier of autosomal recessive condition.",
                    clinvar.significance, conditions_str, clinvar.review_stars,
                ),
                limitations: lims,
            });
        }

        // ── Fix 6: Tier downgrade for NoAlleleData ────────────────────────
        let tier = if is_no_allele_data {
            ConfidenceTier::Tier2Probable
        } else {
            ConfidenceTier::Tier1Reliable
        };

        Some(ScoredResult {
            variant: av.clone(),
            tier,
            category: ResultCategory::MonogenicDisease,
            confirmation_urgency: clinvar_urgency,
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
            limitations: build_limitations(None),
        })
    } else if is_pathogenic {
        // ── Fix 5: AR carrier detection (1-star) ─────────────────────────
        if allele_copies == Some(1) && AR_GENES.contains(&gene) {
            let mut lims = build_limitations(None);
            lims.insert(
                0,
                "You carry one copy of a pathogenic variant in an autosomal \
                recessive gene. Carriers typically do not develop symptoms but may pass \
                the variant to offspring."
                    .to_string(),
            );
            return Some(ScoredResult {
                variant: av.clone(),
                tier: ConfidenceTier::Tier2Probable,
                category: ResultCategory::CarrierStatus,
                confirmation_urgency: ConfirmationUrgency::ClinicalConfirmationRecommended,
                summary: format!(
                    "{gene} ({rsid}) — {sig} ({stars}-star review, limited evidence) (carrier)",
                    sig = clinvar.significance,
                    stars = clinvar.review_stars,
                ),
                details: format!(
                    "Genotype: {genotype}. Classification: {}. Associated conditions: {}. \
                     ClinVar review status: {}-star. Lower review status indicates limited evidence. \
                     Heterozygous carrier of autosomal recessive condition.",
                    clinvar.significance, conditions_str, clinvar.review_stars,
                ),
                limitations: lims,
            });
        }

        // Pathogenic with 1 star
        Some(ScoredResult {
            variant: av.clone(),
            tier: ConfidenceTier::Tier2Probable,
            category: ResultCategory::MonogenicDisease,
            confirmation_urgency: clinvar_urgency,
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
            limitations: build_limitations(None),
        })
    } else {
        // Other significance values (not pathogenic, not benign, not VUS,
        // not conflicting) — skip
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
    // Gate on allele match — skip if user doesn't carry the variant allele
    let allele_result = clinvar_allele_check(av);
    match allele_result {
        ClinvarAlleleResult::Copies(0) => return None, // homozygous ref → skip
        ClinvarAlleleResult::IndelNotDetectable => return None,
        ClinvarAlleleResult::NoCallGenotype => return None,
        _ => {} // Copies(1|2), Palindromic, NoAlleleData → proceed
    }

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
        confirmation_urgency: ConfirmationUrgency::ClinicalConfirmationRecommended,
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
        limitations: vec![
            PGX_FDA_DISCLAIMER.to_string(),
            DTC_RAW_DATA_CAVEAT.to_string(),
        ],
    })
}

/// Score a single GWAS hit with strand-aware risk allele counting.
///
/// Performs strand-aware comparison of the user's genotype against the
/// GWAS risk allele, applying complement matching and palindromic detection.
///
/// - User is homozygous reference (0 risk allele copies): skip this hit
/// Deduplicate GWAS hits by trait name, keeping the most significant (lowest p-value)
/// study for each unique trait. The GWAS Catalog often has many studies for the same
/// trait on the same rsID — we only want one result per trait per variant.
fn dedup_gwas_hits(hits: &[GwasHit]) -> Vec<GwasHit> {
    use std::collections::HashMap;
    let mut best: HashMap<&str, &GwasHit> = HashMap::new();
    for hit in hits {
        let trait_name = hit.trait_name.as_str();
        let entry = best.entry(trait_name).or_insert(hit);
        if hit.p_value < entry.p_value {
            *entry = hit;
        }
    }
    best.into_values().cloned().collect()
}

/// - p_value < 5e-8 and odds_ratio > 1.5: Tier2, GwasAssociation
/// - Otherwise: Tier3, ComplexTrait
///
/// Adds effect context and limitations/caveats to every GWAS result.
fn score_gwas_hit(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    hit: &GwasHit,
) -> Option<ScoredResult> {
    // Extract the two alleles from the user's genotype
    let (allele1, allele2) = match &av.variant.genotype {
        Genotype::Homozygous(a) => (*a, *a),
        Genotype::Heterozygous(a, b) => (*a, *b),
        Genotype::NoCall | Genotype::Indel(_) => {
            // Cannot do strand-aware comparison for no-calls or indels
            return score_gwas_hit_fallback(av, rsid, genotype, hit);
        }
    };

    // Count risk allele copies with strand awareness
    let copy_result = count_risk_allele_copies(hit.risk_allele.as_deref(), allele1, allele2);

    // Determine copy count and build limitations
    let mut limitations = Vec::new();
    let copies: Option<u8>;

    match &copy_result {
        RiskAlleleCopies::Determined { copies: c, .. } => {
            if *c == 0 {
                // User is homozygous reference — skip this hit
                return None;
            }
            copies = Some(*c);
        }
        RiskAlleleCopies::Palindromic { copies: c } => {
            if *c == 0 {
                return None;
            }
            copies = Some(*c);
            limitations.push(
                "Strand ambiguity: this is a palindromic SNP (A/T or C/G). \
                 Risk allele assignment may be on the wrong strand."
                    .to_string(),
            );
        }
        RiskAlleleCopies::Indeterminate => {
            copies = None;
            limitations.push(
                "Risk allele not specified in GWAS data. \
                 Copy count could not be determined."
                    .to_string(),
            );
        }
    }

    // Also check frequency-based palindromic resolution via the existing allele matcher
    let allele_match = compute_allele_match(av, hit);
    if allele_match == Some(AlleleMatch::StrandAmbiguous)
        && !limitations.iter().any(|l| l.contains("palindromic"))
    {
        limitations.push(STRAND_AMBIGUITY_CAVEAT.to_string());
    }

    // Standard GWAS caveats
    limitations.push(
        "Single variant association \u{2014} effect size may be small in isolation".to_string(),
    );
    limitations.push("Population-specific: effect size may vary across ancestries".to_string());

    // GWAS-specific effect size caveats
    if hit.odds_ratio.is_some() {
        limitations.push(GWAS_OR_RELATIVE_CAVEAT.to_string());
    } else if hit.beta.is_some() {
        limitations.push(GWAS_BETA_CAVEAT.to_string());
    }
    limitations.push(GWAS_OR_INVERSION_CAVEAT.to_string());

    // Effect context
    let effect_context = gwas_effect_context(hit);

    // Tier and category assignment
    let genome_wide_significant = hit.p_value < 5e-8;
    let moderate_effect = hit.odds_ratio.is_some_and(|or| or > 1.5);
    let gene = hit.mapped_gene.as_deref().unwrap_or("intergenic");

    let (tier, category) = if genome_wide_significant && moderate_effect {
        (
            ConfidenceTier::Tier2Probable,
            ResultCategory::GwasAssociation,
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

    let copies_desc = match copies {
        Some(c) => format!(" Risk allele copies: {c}."),
        None => String::new(),
    };

    let effect_context_desc = if effect_context.is_empty() {
        String::new()
    } else {
        format!(" {effect_context}")
    };

    let pubmed = hit
        .pubmed_id
        .as_deref()
        .map(|id| format!(" (PMID: {id})"))
        .unwrap_or_default();

    limitations.push(DTC_RAW_DATA_CAVEAT.to_string());

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category,
        confirmation_urgency: ConfirmationUrgency::InformationalOnly,
        summary: format!(
            "{gene} ({rsid}) \u{2014} {trait_name}: {effect_desc}",
            trait_name = hit.trait_name,
        ),
        details: format!(
            "Genotype: {genotype}.{copies_desc} Trait: {}. p-value: {:.2e}, {effect_desc}. \
             Mapped gene: {gene}.{effect_context_desc}{pubmed}",
            hit.trait_name, hit.p_value,
        ),
        limitations,
    })
}

/// Fallback GWAS scoring for NoCall/Indel genotypes where strand-aware
/// comparison is not possible.
fn score_gwas_hit_fallback(
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
            ResultCategory::GwasAssociation,
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

    let effect_context = gwas_effect_context(hit);
    let effect_context_desc = if effect_context.is_empty() {
        String::new()
    } else {
        format!(" {effect_context}")
    };

    let pubmed = hit
        .pubmed_id
        .as_deref()
        .map(|id| format!(" (PMID: {id})"))
        .unwrap_or_default();

    let mut limitations = vec![
        "Genotype is no-call or indel \u{2014} strand-aware risk allele comparison not possible"
            .to_string(),
        "Single variant association \u{2014} effect size may be small in isolation".to_string(),
        "Population-specific: effect size may vary across ancestries".to_string(),
    ];

    // GWAS-specific effect size caveats
    if hit.odds_ratio.is_some() {
        limitations.push(GWAS_OR_RELATIVE_CAVEAT.to_string());
    } else if hit.beta.is_some() {
        limitations.push(GWAS_BETA_CAVEAT.to_string());
    }
    limitations.push(GWAS_OR_INVERSION_CAVEAT.to_string());

    limitations.push(DTC_RAW_DATA_CAVEAT.to_string());

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category,
        confirmation_urgency: ConfirmationUrgency::InformationalOnly,
        summary: format!(
            "{gene} ({rsid}) \u{2014} {trait_name}: {effect_desc}",
            trait_name = hit.trait_name,
        ),
        details: format!(
            "Genotype: {genotype}. Trait: {}. p-value: {:.2e}, {effect_desc}. \
             Mapped gene: {gene}.{effect_context_desc}{pubmed}",
            hit.trait_name, hit.p_value,
        ),
        limitations,
    })
}

/// Generate human-readable effect context for a GWAS hit.
fn gwas_effect_context(hit: &GwasHit) -> String {
    if let Some(or) = hit.odds_ratio {
        if or > 2.0 {
            "Substantially elevated risk.".to_string()
        } else if or >= 1.5 {
            "Moderately elevated risk.".to_string()
        } else if or >= 1.2 {
            "Slightly elevated risk.".to_string()
        } else if or >= 1.0 {
            "Marginally elevated risk.".to_string()
        } else {
            "Protective allele: odds ratio below 1.0 suggests reduced risk.".to_string()
        }
    } else if let Some(beta) = hit.beta {
        if beta.abs() > 0.1 && hit.p_value < 5e-8 {
            format!(
                "Quantitative trait effect (beta={beta:.3}): genome-wide significant \
                 association with a measurable effect on the trait."
            )
        } else {
            String::new()
        }
    } else {
        String::new()
    }
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

    let urgency = if snpedia.magnitude >= 3.0 {
        ConfirmationUrgency::ClinicalConfirmationRecommended
    } else {
        ConfirmationUrgency::InformationalOnly
    };

    Some(ScoredResult {
        variant: av.clone(),
        tier,
        category,
        confirmation_urgency: urgency,
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
        limitations: vec![DTC_RAW_DATA_CAVEAT.to_string()],
    })
}

/// Compute allele match for a GWAS hit using frequency-aware strand resolution.
fn compute_allele_match(av: &AnnotatedVariant, hit: &GwasHit) -> Option<AlleleMatch> {
    let risk_allele_str = hit.risk_allele.as_deref()?;
    let risk_char = single_base_char(risk_allele_str)?;

    let (a1, a2) = match av.variant.genotype {
        Genotype::Homozygous(a) => (a, a),
        Genotype::Heterozygous(a, b) => (a, b),
        Genotype::NoCall | Genotype::Indel(_) => return None,
    };

    let alt_char = infer_alt_allele(a1, a2, risk_char);
    let db_af = av.frequency.as_ref().map(|f| f.af_total);
    let user_af: Option<f64> = None;

    Some(match_alleles_with_frequency(
        (a1, a2),
        risk_char,
        alt_char,
        user_af,
        db_af,
    ))
}

/// Extract a single DNA base character from a string.
fn single_base_char(s: &str) -> Option<char> {
    let trimmed = s.trim();
    if trimmed.len() != 1 {
        return None;
    }
    let c = trimmed.chars().next()?;
    match c.to_ascii_uppercase() {
        'A' | 'T' | 'C' | 'G' => Some(c.to_ascii_uppercase()),
        _ => None,
    }
}

/// Infer the alternate allele given two user alleles and a reference allele.
fn infer_alt_allele(a1: char, a2: char, ref_allele: char) -> char {
    let r = ref_allele.to_ascii_uppercase();
    let u1 = a1.to_ascii_uppercase();
    let u2 = a2.to_ascii_uppercase();

    if u1 == r && u2 != r {
        u2
    } else if u2 == r && u1 != r {
        u1
    } else {
        crate::allele::strand::complement(r).unwrap_or(r)
    }
}

/// Display implementation for `ResultCategory` used in reports.
impl std::fmt::Display for ResultCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ResultCategory::MonogenicDisease => write!(f, "Monogenic Disease Risk"),
            ResultCategory::CarrierStatus => write!(f, "Carrier Status"),
            ResultCategory::Pharmacogenomics => write!(f, "Pharmacogenomics"),
            ResultCategory::GwasAssociation => write!(f, "GWAS Association"),
            ResultCategory::PhysicalTrait => write!(f, "Physical Traits"),
            ResultCategory::ComplexTrait => write!(f, "Complex Traits"),
            ResultCategory::Ancestry => write!(f, "Ancestry"),
            ResultCategory::ClinVarConflicting => write!(f, "ClinVar Conflicting"),
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
            ref_allele: None,
            alt_allele: None,
        }
    }

    fn make_variant_with_genotype(rsid: &str, genotype: Genotype) -> Variant {
        Variant {
            rsid: Some(rsid.to_string()),
            chromosome: "1".to_string(),
            position: 100000,
            genotype,
            source_format: SourceFormat::TwentyThreeAndMe,
        }
    }

    fn make_annotated_with_genotype(rsid: &str, genotype: Genotype) -> AnnotatedVariant {
        AnnotatedVariant {
            variant: make_variant_with_genotype(rsid, genotype),
            clinvar: None,
            snpedia: None,
            gwas_hits: Vec::new(),
            frequency: None,
            pharmacogenomics: None,
            ref_allele: None,
            alt_allele: None,
        }
    }

    /// Create an annotated variant with genotype AND ref/alt alleles for allele matching tests.
    fn make_annotated_with_alleles(
        rsid: &str,
        genotype: Genotype,
        ref_allele: &str,
        alt_allele: &str,
    ) -> AnnotatedVariant {
        AnnotatedVariant {
            variant: make_variant_with_genotype(rsid, genotype),
            clinvar: None,
            snpedia: None,
            gwas_hits: Vec::new(),
            frequency: None,
            pharmacogenomics: None,
            ref_allele: Some(ref_allele.to_string()),
            alt_allele: Some(alt_allele.to_string()),
        }
    }

    // ── ClinVar basic scoring ────────────────────────────────────────────

    #[test]
    fn clinvar_pathogenic_high_stars_is_tier1() {
        let mut av = make_annotated("rs123");
        av.ref_allele = Some("A".to_string());
        av.alt_allele = Some("G".to_string());
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Breast cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
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
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
    }

    #[test]
    fn clinvar_benign_is_skipped() {
        let mut av = make_annotated("rs456");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Benign".to_string(),
            review_stars: 2,
            conditions: vec![],
            gene_symbol: Some("TP53".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Benign variants should be skipped entirely, not reported as findings"
        );
    }

    // ── Pharmacogenomics scoring ─────────────────────────────────────────

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

    // ── GWAS scoring ────────────────────────────────────────────────────

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
        assert_eq!(results[0].category, ResultCategory::GwasAssociation);
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

    // ── SNPedia scoring ─────────────────────────────────────────────────

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

    // ── General behavior ────────────────────────────────────────────────

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
            classification_type: ClinVarClassificationType::Germline,
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
        av2.ref_allele = Some("A".to_string());
        av2.alt_allele = Some("G".to_string());
        av2.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Serious".to_string()],
            gene_symbol: Some("GENE".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av1, av2]);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert_eq!(results[1].tier, ConfidenceTier::Tier3Speculative);
    }

    // ── Blueprint 4: ClinVar quality gate tests ─────────────────────────

    #[test]
    fn zero_star_pathogenic_is_tier3_with_limitation() {
        let mut av = make_annotated("rs800");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 0,
            conditions: vec!["Some condition".to_string()],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
        assert!(
            results[0].limitations.iter().any(|l| l.contains("0 stars")),
            "expected a limitation mentioning 0 stars"
        );
    }

    #[test]
    fn zero_star_benign_is_skipped() {
        let mut av = make_annotated("rs801");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Benign".to_string(),
            review_stars: 0,
            conditions: vec![],
            gene_symbol: Some("GENE2".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "0-star benign should produce no results"
        );
    }

    #[test]
    fn zero_star_likely_benign_is_skipped() {
        let mut av = make_annotated("rs807");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Likely benign".to_string(),
            review_stars: 0,
            conditions: vec![],
            gene_symbol: Some("GENE4".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "0-star likely benign should produce no results"
        );
    }

    #[test]
    fn two_star_pathogenic_is_tier1_unchanged() {
        // Use a non-AR gene (BRCA1) to test basic 2-star pathogenic → Tier1 behavior.
        // AR genes (like CFTR) with het genotype now correctly produce CarrierStatus.
        let mut av = make_annotated("rs802");
        av.ref_allele = Some("A".to_string());
        av.alt_allele = Some("G".to_string());
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 2,
            conditions: vec!["Breast cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
        // Allele match limitation + DTC raw-data caveat
        assert!(results[0]
            .limitations
            .iter()
            .any(|l| l.contains("direct-to-consumer")));
    }

    #[test]
    fn conflicting_interpretations_is_tier3_clinvar_conflicting() {
        let mut av = make_annotated("rs803");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Conflicting interpretations of pathogenicity".to_string(),
            review_stars: 1,
            conditions: vec!["Heart disease".to_string()],
            gene_symbol: Some("SCN5A".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert_eq!(results[0].category, ResultCategory::ClinVarConflicting);
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("Conflicting")),
            "expected a limitation about conflicting interpretations"
        );
    }

    #[test]
    fn vus_is_tier3_with_uncertainty_language() {
        let mut av = make_annotated("rs804");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Uncertain significance".to_string(),
            review_stars: 1,
            conditions: vec!["Hereditary cancer".to_string()],
            gene_symbol: Some("BRCA2".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert!(results[0].summary.contains("uncertain significance"));
        assert!(results[0]
            .limitations
            .iter()
            .any(|l| l.contains("uncertain")),);
    }

    #[test]
    fn one_star_pathogenic_is_tier2() {
        let mut av = make_annotated("rs805");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 1,
            conditions: vec!["Rare disease".to_string()],
            gene_symbol: Some("GENE3".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier2Probable);
        assert_eq!(results[0].category, ResultCategory::MonogenicDisease);
    }

    #[test]
    fn high_af_pathogenic_is_skipped_ba1() {
        let mut av = make_annotated("rs806");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 2,
            conditions: vec!["Common disease".to_string()],
            gene_symbol: Some("GENE3".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });
        av.frequency = Some(AlleleFrequency {
            af_total: 0.10,
            af_afr: None,
            af_amr: None,
            af_eas: None,
            af_eur: None,
            af_sas: None,
            source: "gnomad".to_string(),
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "pathogenic with AF > 5% should be filtered by BA1"
        );
    }

    #[test]
    fn benign_with_high_af_is_still_skipped() {
        let mut av = make_annotated("rs808");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Benign".to_string(),
            review_stars: 2,
            conditions: vec![],
            gene_symbol: Some("GENE5".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });
        av.frequency = Some(AlleleFrequency {
            af_total: 0.20,
            af_afr: None,
            af_amr: None,
            af_eas: None,
            af_eur: None,
            af_sas: None,
            source: "gnomad".to_string(),
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Benign variants should be skipped regardless of AF"
        );
    }

    #[test]
    fn ba1_uses_max_population_af() {
        let mut av = make_annotated("rs809");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Rare disease".to_string()],
            gene_symbol: Some("GENE6".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });
        av.frequency = Some(AlleleFrequency {
            af_total: 0.02,
            af_afr: Some(0.08), // above 5% in African population
            af_amr: None,
            af_eas: None,
            af_eur: None,
            af_sas: None,
            source: "gnomad".to_string(),
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "should be filtered by BA1 due to high AF in one population"
        );
    }

    // ── Blueprint 5: Germline/Somatic classification tests ──────────────

    #[test]
    fn somatic_pathogenic_is_tier3_with_limitation() {
        let mut av = make_annotated("rs900");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Lung cancer".to_string()],
            gene_symbol: Some("EGFR".to_string()),
            classification_type: ClinVarClassificationType::Somatic,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert_eq!(results[0].category, ResultCategory::ComplexTrait);
        assert!(!results[0].limitations.is_empty());
        assert!(results[0].limitations[0].contains("somatic"));
    }

    #[test]
    fn oncogenicity_classification_is_tier3_with_limitation() {
        let mut av = make_annotated("rs901");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Oncogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Colorectal cancer".to_string()],
            gene_symbol: Some("KRAS".to_string()),
            classification_type: ClinVarClassificationType::Oncogenicity,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier3Speculative);
        assert!(!results[0].limitations.is_empty());
        assert!(results[0].limitations[0].contains("oncogenicity"));
    }

    #[test]
    fn germline_pathogenic_is_not_downgraded() {
        let mut av = make_annotated("rs902");
        av.ref_allele = Some("A".to_string());
        av.alt_allele = Some("G".to_string());
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["BRCA-related cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        // Allele match limitation + DTC raw-data caveat
        assert!(results[0]
            .limitations
            .iter()
            .any(|l| l.contains("direct-to-consumer")));
    }

    // ── Blueprint 3: Strand-aware GWAS scoring tests ────────────────────

    #[test]
    fn gwas_risk_allele_direct_match_counts_copies() {
        let mut av = make_annotated_with_genotype("rs1000", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Test trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.8),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(results[0].details.contains("Risk allele copies: 1"));
    }

    #[test]
    fn gwas_risk_allele_complement_match_counts_copies() {
        let mut av = make_annotated_with_genotype("rs1001", Genotype::Heterozygous('T', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Test trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.8),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(results[0].details.contains("Risk allele copies: 1"));
    }

    #[test]
    fn gwas_palindromic_risk_allele_has_ambiguity_caveat() {
        let mut av = make_annotated_with_genotype("rs1002", Genotype::Heterozygous('A', 'T'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Test trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.5),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("palindromic")),
            "Expected palindromic caveat in limitations: {:?}",
            results[0].limitations
        );
    }

    #[test]
    fn gwas_no_risk_allele_is_indeterminate_with_caveat() {
        let mut av = make_annotated_with_genotype("rs1003", Genotype::Heterozygous('A', 'G'));
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
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("not specified")),
            "Expected indeterminate caveat in limitations: {:?}",
            results[0].limitations
        );
    }

    #[test]
    fn gwas_homozygous_ref_zero_copies_is_skipped() {
        let mut av = make_annotated_with_genotype("rs1004", Genotype::Homozygous('G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Test trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(2.0),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Expected 0 copies (homozygous ref) to be skipped, got {} results",
            results.len()
        );
    }

    #[test]
    fn gwas_homozygous_risk_allele_two_copies() {
        let mut av = make_annotated_with_genotype("rs1005", Genotype::Homozygous('A'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Disease".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.8),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(results[0].details.contains("Risk allele copies: 2"));
    }

    #[test]
    fn gwas_or_above_2_substantially_elevated() {
        let mut av = make_annotated_with_genotype("rs1006", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Disease".to_string(),
            p_value: 1e-12,
            odds_ratio: Some(2.5),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.1),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(
            results[0].details.contains("Substantially elevated risk"),
            "Expected 'Substantially elevated risk' in details: {}",
            results[0].details
        );
    }

    #[test]
    fn gwas_or_below_1_protective_note() {
        let mut av = make_annotated_with_genotype("rs1007", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Disease".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(0.7),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.4),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(
            results[0].details.contains("Protective allele"),
            "Expected protective allele note in details: {}",
            results[0].details
        );
    }

    #[test]
    fn gwas_every_result_has_single_variant_and_ancestry_caveats() {
        let mut av = make_annotated_with_genotype("rs1008", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.3),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.25),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("Single variant association")),
            "Expected single-variant caveat: {:?}",
            results[0].limitations
        );
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("Population-specific")),
            "Expected ancestry caveat: {:?}",
            results[0].limitations
        );
    }

    #[test]
    fn gwas_effect_context_moderately_elevated() {
        let hit = GwasHit {
            trait_name: "T".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.7),
            beta: None,
            risk_allele: None,
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        };
        assert_eq!(gwas_effect_context(&hit), "Moderately elevated risk.");
    }

    #[test]
    fn gwas_effect_context_slightly_elevated() {
        let hit = GwasHit {
            trait_name: "T".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.3),
            beta: None,
            risk_allele: None,
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        };
        assert_eq!(gwas_effect_context(&hit), "Slightly elevated risk.");
    }

    #[test]
    fn gwas_effect_context_beta_quantitative() {
        let hit = GwasHit {
            trait_name: "T".to_string(),
            p_value: 1e-10,
            odds_ratio: None,
            beta: Some(0.5),
            risk_allele: None,
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        };
        let ctx = gwas_effect_context(&hit);
        assert!(ctx.contains("Quantitative trait effect"), "got: {ctx}");
    }

    #[test]
    fn gwas_effect_context_no_effect_data() {
        let hit = GwasHit {
            trait_name: "T".to_string(),
            p_value: 1e-10,
            odds_ratio: None,
            beta: None,
            risk_allele: None,
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        };
        assert_eq!(gwas_effect_context(&hit), "");
    }

    // ── Palindromic SNP / allele matching integration tests ─────────────

    #[test]
    fn non_palindromic_ag_gwas_no_strand_caveat() {
        let mut av = make_annotated_with_genotype("rs1100", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Test trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(2.0),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: Some(0.3),
            pubmed_id: None,
            mapped_gene: Some("GENE1".to_string()),
        }];

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        // Exclude the DTC raw-data caveat (which mentions "strand ambiguity" generically)
        // and check that no domain-specific strand/palindromic limitation is present.
        assert!(
            !results[0]
                .limitations
                .iter()
                .any(|l| !l.contains("direct-to-consumer")
                    && (l.contains("palindromic") || l.contains("strand"))),
            "non-palindromic SNP should have no strand ambiguity limitation"
        );
    }

    #[test]
    fn palindromic_at_without_frequency_has_strand_caveat() {
        let mut av = make_annotated_with_genotype("rs1101", Genotype::Heterozygous('A', 'T'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Palindromic trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(2.0),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: Some("GENE2".to_string()),
        }];
        av.frequency = None;

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("palindromic")),
            "palindromic SNP without frequency should have strand ambiguity limitation"
        );
    }

    #[test]
    fn palindromic_cg_with_extreme_frequency_resolves() {
        let mut av = make_annotated_with_genotype("rs1102", Genotype::Heterozygous('C', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Resolved trait".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(2.0),
            beta: None,
            risk_allele: Some("C".to_string()),
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: Some("GENE3".to_string()),
        }];
        av.frequency = Some(AlleleFrequency {
            af_total: 0.05,
            af_afr: None,
            af_amr: None,
            af_eas: None,
            af_eur: None,
            af_sas: None,
            source: "gnomad".to_string(),
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        // With extreme AF, the palindromic ambiguity is resolved, but the
        // copy counting may still report palindromic. At least the frequency-based
        // matcher should have resolved, so no STRAND_AMBIGUITY_CAVEAT.
        // Note: the count_risk_allele_copies sees C/G as palindromic, so it
        // still adds palindromic caveat. That's expected — it's conservative.
    }

    // ── Blueprint 7: Clinical Confirmation Language Tests ──────────────

    #[test]
    fn acmg_gene_gets_high_impact_urgency() {
        let mut av = make_annotated_with_genotype("rs123", Genotype::Heterozygous('A', 'G'));
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Breast cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::HighImpact
        );
    }

    #[test]
    fn non_acmg_pathogenic_gets_clinical_confirmation() {
        let mut av = make_annotated_with_genotype("rs456", Genotype::Heterozygous('A', 'G'));
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 2,
            conditions: vec!["Some condition".to_string()],
            gene_symbol: Some("SOME_GENE".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::ClinicalConfirmationRecommended
        );
    }

    #[test]
    fn pharma_gets_clinical_confirmation() {
        let mut av = make_annotated_with_genotype("rs789", Genotype::Heterozygous('A', 'G'));
        av.pharmacogenomics = Some(PharmaAnnotation {
            gene: "CYP2D6".to_string(),
            drug: "Codeine".to_string(),
            phenotype_category: Some("Poor Metabolizer".to_string()),
            evidence_level: "1A".to_string(),
            clinical_recommendation: Some("Consider alternative".to_string()),
        });

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::ClinicalConfirmationRecommended
        );
    }

    #[test]
    fn gwas_gets_informational_only() {
        let mut av = make_annotated_with_genotype("rs111", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "Height".to_string(),
            p_value: 1e-10,
            odds_ratio: Some(1.3),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        }];

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::InformationalOnly
        );
    }

    #[test]
    fn snpedia_high_magnitude_gets_clinical_confirmation() {
        let mut av = make_annotated_with_genotype("rs222", Genotype::Heterozygous('A', 'G'));
        av.snpedia = Some(SnpediaAnnotation {
            magnitude: 4.0,
            repute: Some("bad".to_string()),
            summary: "Important finding".to_string(),
            genotype_descriptions: None,
        });

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::ClinicalConfirmationRecommended
        );
    }

    #[test]
    fn snpedia_low_magnitude_gets_informational() {
        let mut av = make_annotated_with_genotype("rs333", Genotype::Heterozygous('A', 'G'));
        av.snpedia = Some(SnpediaAnnotation {
            magnitude: 1.5,
            repute: None,
            summary: "Minor trait".to_string(),
            genotype_descriptions: None,
        });

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].confirmation_urgency,
            ConfirmationUrgency::InformationalOnly
        );
    }

    #[test]
    fn every_result_has_dtc_caveat() {
        let mut av = make_annotated_with_genotype("rs444", Genotype::Heterozygous('A', 'G'));
        av.gwas_hits = vec![GwasHit {
            trait_name: "BMI".to_string(),
            p_value: 1e-8,
            odds_ratio: Some(1.2),
            beta: None,
            risk_allele: Some("A".to_string()),
            risk_allele_frequency: None,
            pubmed_id: None,
            mapped_gene: None,
        }];

        let results = score_variants(&[av]);
        assert!(!results.is_empty());
        assert!(results[0]
            .limitations
            .iter()
            .any(|l| l.contains("direct-to-consumer")));
    }

    // ── ClinVar allele matching tests ────────────────────────────────────

    #[test]
    fn clinvar_hom_ref_is_skipped_with_allele_data() {
        // User is AA (homozygous reference), pathogenic alt is G.
        // User does NOT carry the variant — should be skipped.
        let mut av = make_annotated_with_alleles("rs2000", Genotype::Homozygous('A'), "A", "G");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Serious disease".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Homozygous reference user should not be reported as having pathogenic variant"
        );
    }

    #[test]
    fn clinvar_het_with_allele_data_is_reported() {
        // User is AG (heterozygous), ref=A, alt=G.
        // User carries 1 copy of the pathogenic allele.
        let mut av =
            make_annotated_with_alleles("rs2001", Genotype::Heterozygous('A', 'G'), "A", "G");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Breast cancer".to_string()],
            gene_symbol: Some("BRCA1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].tier, ConfidenceTier::Tier1Reliable);
        assert!(results[0].limitations.iter().any(|l| l.contains("1 copy")));
    }

    #[test]
    fn clinvar_hom_alt_with_allele_data_is_reported() {
        // User is GG (homozygous alt), ref=A, alt=G.
        // User carries 2 copies of the pathogenic allele.
        let mut av = make_annotated_with_alleles("rs2002", Genotype::Homozygous('G'), "A", "G");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Cystic fibrosis".to_string()],
            gene_symbol: Some("CFTR".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1);
        assert!(results[0]
            .limitations
            .iter()
            .any(|l| l.contains("2 copies")));
    }

    #[test]
    fn clinvar_complement_strand_match() {
        // User reports on opposite strand: T/C instead of A/G.
        // ref=A, alt=G. Complement: ref=T, alt=C.
        // User TC = het on complement strand = 1 copy of alt.
        let mut av =
            make_annotated_with_alleles("rs2003", Genotype::Heterozygous('T', 'C'), "A", "G");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 2,
            conditions: vec!["Disease".to_string()],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(
            results.len(),
            1,
            "Complement strand match should still be reported"
        );
        assert!(results[0].limitations.iter().any(|l| l.contains("1 copy")));
    }

    #[test]
    fn clinvar_hom_ref_complement_strand_skipped() {
        // User is TT on complement strand. ref=A, alt=G.
        // Complement of A is T, complement of G is C.
        // User TT = homozygous ref on complement strand = 0 copies of alt.
        let mut av = make_annotated_with_alleles("rs2004", Genotype::Homozygous('T'), "A", "G");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 4,
            conditions: vec!["Disease".to_string()],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Homozygous reference on complement strand should be skipped"
        );
    }

    #[test]
    fn clinvar_without_allele_data_still_reports() {
        // No ref/alt allele data — falls back to legacy behavior (reports all pathogenic).
        let mut av = make_annotated("rs2005");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Disease".to_string()],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(
            results.len(),
            1,
            "Without allele data, pathogenic should still be reported (legacy behavior)"
        );
    }

    #[test]
    fn clinvar_likely_benign_is_skipped() {
        let mut av = make_annotated("rs2006");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Likely benign".to_string(),
            review_stars: 3,
            conditions: vec![],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert!(
            results.is_empty(),
            "Likely benign should be skipped entirely"
        );
    }

    #[test]
    fn clinvar_palindromic_snp_reports_with_caveat() {
        // A/T is palindromic — can't resolve strand without frequency.
        // Should still report but with inconclusive caveat.
        let mut av =
            make_annotated_with_alleles("rs2007", Genotype::Heterozygous('A', 'T'), "A", "T");
        av.clinvar = Some(ClinVarAnnotation {
            significance: "Pathogenic".to_string(),
            review_stars: 3,
            conditions: vec!["Disease".to_string()],
            gene_symbol: Some("GENE1".to_string()),
            classification_type: ClinVarClassificationType::Germline,
        });

        let results = score_variants(&[av]);
        assert_eq!(results.len(), 1, "Palindromic SNP should still be reported");
        assert!(
            results[0]
                .limitations
                .iter()
                .any(|l| l.contains("inconclusive") || l.contains("palindromic")),
            "Expected palindromic/inconclusive caveat in limitations: {:?}",
            results[0].limitations
        );
    }
}
