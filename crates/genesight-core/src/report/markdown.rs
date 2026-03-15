//! Markdown report renderer.
//!
//! Generates a structured Markdown report with results grouped by
//! confidence tier and category.

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::models::confidence::ConfidenceTier;
use crate::models::report::{ConfirmationUrgency, Report, ResultCategory, ScoredResult};

use super::ReportError;

/// Render a report as Markdown.
///
/// The output includes:
/// - Header with title and generation date
/// - Medical disclaimer
/// - Summary statistics
/// - Results grouped by tier, then by category
/// - Data source attributions
pub fn render(report: &Report) -> Result<String, ReportError> {
    let mut out = String::with_capacity(4096);

    write_header(&mut out);
    write_disclaimer(&mut out, &report.disclaimer);
    write_dtc_context(&mut out, &report.dtc_context);
    write_summary(&mut out, report);
    write_fda_pgx_disclaimer(&mut out, &report.results);
    write_results(&mut out, &report.results);
    write_attributions(&mut out, &report.attributions);

    Ok(out)
}

fn write_header(out: &mut String) {
    out.push_str("# GeneSight Analysis Report\n\n");
    out.push_str("---\n\n");
}

fn write_disclaimer(out: &mut String, disclaimer: &str) {
    out.push_str("> **Medical Disclaimer**\n>\n");
    for line in disclaimer.lines() {
        let _ = writeln!(out, "> {line}");
    }
    out.push_str("\n---\n\n");
}

/// Render the DTC context statement if non-empty.
fn write_dtc_context(out: &mut String, dtc_context: &str) {
    if dtc_context.is_empty() {
        return;
    }
    out.push_str("> **Direct-to-Consumer Data Context**\n>\n");
    for line in dtc_context.lines() {
        let _ = writeln!(out, "> {line}");
    }
    out.push_str("\n---\n\n");
}

/// Render the FDA PGx disclaimer if any pharmacogenomic results exist.
fn write_fda_pgx_disclaimer(out: &mut String, results: &[ScoredResult]) {
    let has_pgx = results
        .iter()
        .any(|r| r.category == ResultCategory::Pharmacogenomics);
    if !has_pgx {
        return;
    }
    out.push_str("> **FDA Notice: Pharmacogenomic Results**\n>\n");
    out.push_str(
        "> Pharmacogenomic results from consumer genotyping arrays have NOT been reviewed or \
         approved by the U.S. Food and Drug Administration (FDA) for clinical use. Do not alter \
         any medication regimen based solely on these results. Consult a healthcare provider or \
         clinical pharmacogenomics service for validated testing.\n\n",
    );
}

/// Format urgency level for markdown output.
fn urgency_text(urgency: ConfirmationUrgency) -> &'static str {
    match urgency {
        ConfirmationUrgency::HighImpact => {
            "**High Impact:** Clinical-grade confirmation strongly recommended (ACMG actionable gene)."
        }
        ConfirmationUrgency::ClinicalConfirmationRecommended => {
            "**Clinical Confirmation Recommended:** This finding should be confirmed through clinical-grade testing before any medical decisions."
        }
        ConfirmationUrgency::InformationalOnly => {
            "**Informational Only:** No clinical action warranted from DTC data alone."
        }
    }
}

fn write_summary(out: &mut String, report: &Report) {
    out.push_str("## Summary\n\n");
    let _ = writeln!(out, "| Metric | Value |");
    let _ = writeln!(out, "|--------|------:|");
    let _ = writeln!(
        out,
        "| Total variants in file | {} |",
        report.total_variants
    );
    let _ = writeln!(
        out,
        "| Annotated variants | {} |",
        report.annotated_variants
    );
    let _ = writeln!(out, "| Scored results | {} |", report.results.len());

    // Count by tier
    let tier1 = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier1Reliable)
        .count();
    let tier2 = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier2Probable)
        .count();
    let tier3 = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier3Speculative)
        .count();

    let _ = writeln!(out, "| Tier 1 (Reliable) | {tier1} |");
    let _ = writeln!(out, "| Tier 2 (Probable) | {tier2} |");
    let _ = writeln!(out, "| Tier 3 (Speculative) | {tier3} |");
    let _ = writeln!(out, "| Input assembly | {} |", report.input_assembly);
    let _ = writeln!(out, "| Database assembly | {} |", report.db_assembly);
    out.push('\n');

    write_assembly_warnings(out, &report.assembly_warnings);
}

fn write_assembly_warnings(out: &mut String, warnings: &[String]) {
    if warnings.is_empty() {
        return;
    }

    out.push_str("### Assembly Warnings\n\n");
    for warning in warnings {
        let _ = writeln!(out, "> **Warning:** {warning}");
        out.push_str(">\n");
    }
    out.push('\n');
}

fn write_results(out: &mut String, results: &[ScoredResult]) {
    if results.is_empty() {
        out.push_str("## Results\n\nNo significant findings.\n\n");
        return;
    }

    out.push_str("## Results\n\n");

    // Group by tier, then by category
    let grouped = group_results(results);

    for (tier, categories) in &grouped {
        let _ = writeln!(out, "### {tier}\n");

        for (category, items) in categories {
            let _ = writeln!(out, "#### {category}\n");

            for result in items {
                let rsid = result.variant.variant.rsid.as_deref().unwrap_or("unknown");
                let genotype = result.variant.variant.genotype.to_string();
                let gene = extract_gene(result);
                let tier_badge = tier_badge(*tier);

                let _ = writeln!(out, "**{rsid}** {tier_badge}");
                if !gene.is_empty() {
                    let _ = writeln!(out, "- **Gene:** {gene}");
                }
                let _ = writeln!(out, "- **Genotype:** {genotype}");
                let _ = writeln!(out, "- **Summary:** {}", result.summary);
                let _ = writeln!(out, "- **Details:** {}", result.details);
                // Urgency indicator
                match result.confirmation_urgency {
                    ConfirmationUrgency::HighImpact
                    | ConfirmationUrgency::ClinicalConfirmationRecommended => {
                        let _ =
                            writeln!(out, "\n> {}\n", urgency_text(result.confirmation_urgency));
                    }
                    ConfirmationUrgency::InformationalOnly => {
                        let _ = writeln!(out, "- {}", urgency_text(result.confirmation_urgency));
                    }
                }
                if !result.limitations.is_empty() {
                    let _ = writeln!(out, "- **Limitations:**");
                    for limitation in &result.limitations {
                        let _ = writeln!(out, "  - {limitation}");
                    }
                }
                out.push('\n');
            }
        }
    }
}

fn write_attributions(out: &mut String, attributions: &[String]) {
    out.push_str("---\n\n");
    out.push_str("## Data Sources & Attributions\n\n");
    for attr in attributions {
        let _ = writeln!(out, "- {attr}");
    }
    out.push('\n');
    out.push_str("---\n\n");
    out.push_str("*Generated by GeneSight*\n");
}

/// Group results by tier (ordered) and then by category (ordered).
fn group_results(
    results: &[ScoredResult],
) -> BTreeMap<ConfidenceTier, BTreeMap<ResultCategory, Vec<&ScoredResult>>> {
    let mut grouped: BTreeMap<ConfidenceTier, BTreeMap<ResultCategory, Vec<&ScoredResult>>> =
        BTreeMap::new();

    for result in results {
        grouped
            .entry(result.tier)
            .or_default()
            .entry(result.category)
            .or_default()
            .push(result);
    }

    grouped
}

/// Create a tier badge string.
fn tier_badge(tier: ConfidenceTier) -> &'static str {
    match tier {
        ConfidenceTier::Tier1Reliable => "`[Tier 1: Reliable]`",
        ConfidenceTier::Tier2Probable => "`[Tier 2: Probable]`",
        ConfidenceTier::Tier3Speculative => "`[Tier 3: Speculative]`",
    }
}

/// Extract the most relevant gene name from a scored result.
fn extract_gene(result: &ScoredResult) -> String {
    if let Some(clinvar) = &result.variant.clinvar {
        if let Some(gene) = &clinvar.gene_symbol {
            return gene.clone();
        }
    }
    if let Some(pharma) = &result.variant.pharmacogenomics {
        return pharma.gene.clone();
    }
    for hit in &result.variant.gwas_hits {
        if let Some(gene) = &hit.mapped_gene {
            return gene.clone();
        }
    }
    String::new()
}

// Implement Ord for ResultCategory so BTreeMap ordering works
impl PartialOrd for ResultCategory {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ResultCategory {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.sort_key().cmp(&other.sort_key())
    }
}

impl ResultCategory {
    /// Assign a sort order for report display.
    fn sort_key(&self) -> u8 {
        match self {
            ResultCategory::MonogenicDisease => 0,
            ResultCategory::CarrierStatus => 1,
            ResultCategory::Pharmacogenomics => 2,
            ResultCategory::GwasAssociation => 3,
            ResultCategory::PhysicalTrait => 4,
            ResultCategory::ComplexTrait => 5,
            ResultCategory::Ancestry => 6,
            ResultCategory::ClinVarConflicting => 7,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::annotation::*;
    use crate::models::assembly::GenomeAssembly;
    use crate::models::report::ConfirmationUrgency;
    use crate::models::variant::{Genotype, SourceFormat, Variant};

    fn make_test_report() -> Report {
        let variant = Variant {
            rsid: Some("rs123".to_string()),
            chromosome: "17".to_string(),
            position: 43093802,
            genotype: Genotype::Heterozygous('A', 'G'),
            source_format: SourceFormat::TwentyThreeAndMe,
        };

        let annotated = AnnotatedVariant {
            variant,
            clinvar: Some(ClinVarAnnotation {
                significance: "Pathogenic".to_string(),
                review_stars: 3,
                conditions: vec!["Breast cancer".to_string()],
                gene_symbol: Some("BRCA1".to_string()),
                classification_type: crate::models::annotation::ClinVarClassificationType::Germline,
            }),
            snpedia: None,
            gwas_hits: Vec::new(),
            frequency: None,
            pharmacogenomics: None,
            ref_allele: None,
            alt_allele: None,
        };

        Report {
            total_variants: 600000,
            annotated_variants: 1500,
            results: vec![ScoredResult {
                variant: annotated,
                tier: ConfidenceTier::Tier1Reliable,
                category: ResultCategory::MonogenicDisease,
                confirmation_urgency: ConfirmationUrgency::HighImpact,
                summary: "BRCA1 (rs123) — Pathogenic (3-star review)".to_string(),
                details: "Genotype: AG. Classification: Pathogenic.".to_string(),
                limitations: Vec::new(),
            }],
            attributions: vec![
                "ClinVar: NCBI/NLM (public domain)".to_string(),
                "GWAS Catalog: EMBL-EBI (open access)".to_string(),
            ],
            disclaimer: "This is not medical advice. Consult a healthcare professional."
                .to_string(),
            dtc_context: String::new(),
            input_assembly: GenomeAssembly::GRCh37,
            db_assembly: GenomeAssembly::GRCh37,
            assembly_warnings: Vec::new(),
        }
    }

    #[test]
    fn render_contains_header() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("# GeneSight Analysis Report"));
    }

    #[test]
    fn render_contains_disclaimer() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("Medical Disclaimer"));
        assert!(md.contains("not medical advice"));
    }

    #[test]
    fn render_contains_summary_stats() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("600000"));
        assert!(md.contains("1500"));
    }

    #[test]
    fn render_contains_results() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("rs123"));
        assert!(md.contains("BRCA1"));
        assert!(md.contains("Tier 1: Reliable"));
    }

    #[test]
    fn render_contains_attributions() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("ClinVar"));
        assert!(md.contains("GWAS Catalog"));
    }

    #[test]
    fn render_empty_results() {
        let report = Report {
            total_variants: 0,
            annotated_variants: 0,
            results: vec![],
            attributions: vec![],
            disclaimer: "Disclaimer.".to_string(),
            dtc_context: String::new(),
            input_assembly: GenomeAssembly::Unknown,
            db_assembly: GenomeAssembly::Unknown,
            assembly_warnings: Vec::new(),
        };
        let md = render(&report).expect("render");
        assert!(md.contains("No significant findings"));
    }

    #[test]
    fn render_contains_assembly_info() {
        let report = make_test_report();
        let md = render(&report).expect("render");
        assert!(md.contains("GRCh37 (hg19)"));
        assert!(md.contains("Input assembly"));
        assert!(md.contains("Database assembly"));
    }

    #[test]
    fn render_contains_assembly_warnings() {
        let mut report = make_test_report();
        report.assembly_warnings = vec!["Assembly mismatch detected".to_string()];
        let md = render(&report).expect("render");
        assert!(md.contains("Assembly Warnings"));
        assert!(md.contains("Assembly mismatch detected"));
    }
}
