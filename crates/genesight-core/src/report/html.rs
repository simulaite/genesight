//! HTML report renderer.
//!
//! Generates a self-contained HTML page with embedded CSS for displaying
//! GeneSight analysis results. No external dependencies or assets required.

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::models::confidence::ConfidenceTier;
use crate::models::report::{Report, ResultCategory, ScoredResult};

use super::ReportError;

/// Render a report as a self-contained HTML page.
///
/// The output is a complete HTML document with embedded CSS styling.
/// No external stylesheets, scripts, or assets are referenced.
pub fn render(report: &Report) -> Result<String, ReportError> {
    let mut out = String::with_capacity(8192);

    write_html_head(&mut out);
    write_body_open(&mut out);
    write_header(&mut out);
    write_disclaimer(&mut out, &report.disclaimer);
    write_summary(&mut out, report);
    write_assembly_warnings(&mut out, &report.assembly_warnings);
    write_results(&mut out, &report.results);
    write_attributions(&mut out, &report.attributions);
    write_body_close(&mut out);

    Ok(out)
}

fn write_html_head(out: &mut String) {
    out.push_str(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GeneSight Analysis Report</title>
<style>
:root {
    --color-bg: #fafafa;
    --color-surface: #ffffff;
    --color-text: #1a1a2e;
    --color-muted: #6b7280;
    --color-border: #e5e7eb;
    --color-tier1: #065f46;
    --color-tier1-bg: #d1fae5;
    --color-tier2: #92400e;
    --color-tier2-bg: #fef3c7;
    --color-tier3: #6b7280;
    --color-tier3-bg: #f3f4f6;
    --color-accent: #1e40af;
    --color-disclaimer-bg: #fef2f2;
    --color-disclaimer-border: #fca5a5;
}

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    background: var(--color-bg);
    color: var(--color-text);
    line-height: 1.6;
    max-width: 900px;
    margin: 0 auto;
    padding: 2rem 1.5rem;
}

h1 { font-size: 1.75rem; font-weight: 700; margin-bottom: 0.5rem; }
h2 { font-size: 1.35rem; font-weight: 600; margin: 2rem 0 1rem; border-bottom: 2px solid var(--color-border); padding-bottom: 0.4rem; }
h3 { font-size: 1.15rem; font-weight: 600; margin: 1.5rem 0 0.75rem; color: var(--color-accent); }
h4 { font-size: 1rem; font-weight: 600; margin: 1rem 0 0.5rem; }

.disclaimer {
    background: var(--color-disclaimer-bg);
    border-left: 4px solid var(--color-disclaimer-border);
    padding: 1rem 1.25rem;
    margin: 1.5rem 0;
    border-radius: 0 6px 6px 0;
    font-size: 0.9rem;
}
.disclaimer strong { display: block; margin-bottom: 0.3rem; }

table { width: 100%; border-collapse: collapse; margin: 1rem 0; }
th, td { padding: 0.5rem 0.75rem; text-align: left; border-bottom: 1px solid var(--color-border); }
th { background: var(--color-surface); font-weight: 600; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.03em; color: var(--color-muted); }
td:last-child { text-align: right; font-variant-numeric: tabular-nums; }

.result-card {
    background: var(--color-surface);
    border: 1px solid var(--color-border);
    border-radius: 8px;
    padding: 1rem 1.25rem;
    margin: 0.75rem 0;
    box-shadow: 0 1px 2px rgba(0,0,0,0.04);
}
.result-card .rsid { font-weight: 700; font-size: 1rem; }
.result-card .meta { font-size: 0.85rem; color: var(--color-muted); margin: 0.25rem 0; }
.result-card .summary { margin: 0.5rem 0; }
.result-card .details { font-size: 0.9rem; color: var(--color-muted); margin-top: 0.5rem; }

.badge {
    display: inline-block;
    font-size: 0.75rem;
    font-weight: 600;
    padding: 0.15rem 0.5rem;
    border-radius: 9999px;
    vertical-align: middle;
    margin-left: 0.4rem;
}
.badge-tier1 { background: var(--color-tier1-bg); color: var(--color-tier1); }
.badge-tier2 { background: var(--color-tier2-bg); color: var(--color-tier2); }
.badge-tier3 { background: var(--color-tier3-bg); color: var(--color-tier3); }

.attributions { font-size: 0.85rem; color: var(--color-muted); margin-top: 2rem; }
.attributions ul { list-style: disc; padding-left: 1.5rem; }
.attributions li { margin: 0.2rem 0; }

.footer { margin-top: 2rem; padding-top: 1rem; border-top: 1px solid var(--color-border); font-size: 0.8rem; color: var(--color-muted); text-align: center; }

/* Summary split view — same pattern as findings-split */
.summary-split {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 0;
    margin: 1rem 0;
    background: var(--color-surface);
    border: 1px solid var(--color-border);
    border-radius: 8px;
    overflow: hidden;
    box-shadow: 0 1px 2px rgba(0,0,0,0.04);
}
@media (max-width: 700px) { .summary-split { grid-template-columns: 1fr; } }

.summary-left {
    border-right: 1px solid var(--color-border);
    padding: 0;
}
@media (max-width: 700px) { .summary-left { border-right: none; border-bottom: 1px solid var(--color-border); } }
.summary-left table { margin: 0; }
.summary-left td:first-child { font-weight: 600; }
.summary-left td:last-child { text-align: right; font-weight: 700; font-size: 1.05rem; font-variant-numeric: tabular-nums; }
.summary-left .tier-row-1 td { border-left: 3px solid var(--color-tier1); }
.summary-left .tier-row-2 td { border-left: 3px solid var(--color-tier2); }
.summary-left .tier-row-3 td { border-left: 3px solid var(--color-tier3); }
.summary-left .tier-row td:first-child .badge { margin-left: 0; margin-right: 0.4rem; }

.summary-right {
    padding: 1rem 1.25rem;
    overflow-y: auto;
    max-height: 500px;
}
.summary-right h4 { font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--color-muted); margin: 0 0 0.5rem; }
.summary-section { margin-bottom: 1rem; }
.summary-section:last-child { margin-bottom: 0; }
.summary-section p { font-size: 0.88rem; color: var(--color-muted); margin: 0.2rem 0; line-height: 1.5; }

.coverage-bar-container { background: var(--color-border); border-radius: 4px; height: 8px; margin: 0.4rem 0 0.1rem; overflow: hidden; }
.coverage-bar { background: var(--color-accent); height: 100%; border-radius: 4px; }

.tier-explain { padding: 0.4rem 0.6rem; border-radius: 6px; font-size: 0.85rem; margin: 0.3rem 0; }
.tier-explain-1 { background: var(--color-tier1-bg); color: var(--color-tier1); }
.tier-explain-2 { background: var(--color-tier2-bg); color: var(--color-tier2); }
.tier-explain-3 { background: var(--color-tier3-bg); color: var(--color-tier3); }
.tier-explain strong { display: block; font-size: 0.8rem; margin-bottom: 0.1rem; }

.cat-bar-row { display: flex; align-items: center; gap: 0.5rem; font-size: 0.85rem; margin: 0.2rem 0; }
.cat-bar-label { width: 120px; flex-shrink: 0; }
.cat-bar-bg { flex: 1; background: var(--color-border); border-radius: 3px; height: 8px; }
.cat-bar { height: 100%; border-radius: 3px; }
.cat-bar-count { width: 30px; text-align: right; font-weight: 600; font-variant-numeric: tabular-nums; }

/* Tier 1 findings split view */
.findings-split {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 0;
    margin: 1rem 0;
    background: var(--color-surface);
    border: 1px solid var(--color-border);
    border-radius: 8px;
    overflow: hidden;
    box-shadow: 0 1px 2px rgba(0,0,0,0.04);
}
@media (max-width: 700px) { .findings-split { grid-template-columns: 1fr; } }

.findings-list {
    border-right: 1px solid var(--color-border);
    max-height: 600px;
    overflow-y: auto;
}
@media (max-width: 700px) { .findings-list { border-right: none; border-bottom: 1px solid var(--color-border); max-height: 400px; } }

.findings-list table { margin: 0; }
.findings-list th { position: sticky; top: 0; background: var(--color-surface); z-index: 1; }
.findings-list td { font-size: 0.9rem; cursor: default; }
.findings-list tr:hover td { background: var(--color-tier1-bg); }
.findings-list .cat-label { font-size: 0.75rem; color: var(--color-muted); }
.findings-list .gene-name { font-weight: 700; }
.findings-list .rsid-col { font-family: monospace; font-size: 0.85rem; color: var(--color-muted); }

.findings-detail {
    padding: 1rem 1.25rem;
    max-height: 600px;
    overflow-y: auto;
}
.findings-detail h4 { margin: 0 0 0.75rem; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--color-muted); }

.finding-entry { margin-bottom: 1.25rem; padding-bottom: 1.25rem; border-bottom: 1px solid var(--color-border); }
.finding-entry:last-child { border-bottom: none; margin-bottom: 0; padding-bottom: 0; }
.finding-entry .fe-header { display: flex; align-items: baseline; gap: 0.5rem; margin-bottom: 0.4rem; }
.finding-entry .fe-gene { font-weight: 700; font-size: 1rem; }
.finding-entry .fe-rsid { font-family: monospace; font-size: 0.85rem; color: var(--color-muted); }
.finding-entry .fe-summary { margin: 0.3rem 0; font-size: 0.92rem; }

/* Annotation detail blocks within findings */
.anno-block { margin: 0.5rem 0 0; padding: 0.5rem 0.75rem; background: var(--color-bg); border-radius: 6px; font-size: 0.85rem; }
.anno-block summary { cursor: pointer; font-weight: 600; font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.03em; color: var(--color-accent); }
.anno-block[open] summary { margin-bottom: 0.4rem; }
.anno-row { display: flex; justify-content: space-between; padding: 0.15rem 0; }
.anno-label { color: var(--color-muted); }
.anno-value { font-weight: 500; text-align: right; }
.anno-list { list-style: disc; padding-left: 1.2rem; margin: 0.2rem 0; font-size: 0.85rem; }
.anno-list li { margin: 0.1rem 0; }

/* Stars for ClinVar review status */
.stars { color: #f59e0b; letter-spacing: 0.1em; font-size: 0.95rem; }

/* Frequency bar */
.freq-bar-bg { background: var(--color-border); border-radius: 3px; height: 6px; width: 100%; margin: 0.15rem 0; }
.freq-bar { background: var(--color-accent); height: 100%; border-radius: 3px; min-width: 1px; }
</style>
</head>
"#,
    );
}

fn write_body_open(out: &mut String) {
    out.push_str("<body>\n");
}

fn write_header(out: &mut String) {
    out.push_str("<h1>GeneSight Analysis Report</h1>\n");
    out.push_str("<hr>\n");
}

fn write_disclaimer(out: &mut String, disclaimer: &str) {
    out.push_str("<div class=\"disclaimer\">\n<strong>Medical Disclaimer</strong>\n");
    let escaped = html_escape(disclaimer);
    for line in escaped.lines() {
        let _ = writeln!(out, "<p>{line}</p>");
    }
    out.push_str("</div>\n");
}

fn write_summary(out: &mut String, report: &Report) {
    use std::collections::BTreeMap;

    let tier1_count = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier1Reliable)
        .count();
    let tier2_count = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier2Probable)
        .count();
    let tier3_count = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier3Speculative)
        .count();

    // Count results by category
    let mut category_counts: BTreeMap<ResultCategory, usize> = BTreeMap::new();
    for result in &report.results {
        *category_counts.entry(result.category).or_default() += 1;
    }

    // Compute coverage percentage
    let coverage_pct = if report.total_variants > 0 {
        (report.annotated_variants as f64 / report.total_variants as f64) * 100.0
    } else {
        0.0
    };

    let max_cat_count = category_counts.values().copied().max().unwrap_or(1);

    out.push_str("<h2>Summary</h2>\n");

    // --- Split view: left = metrics table, right = context details ---
    out.push_str("<div class=\"summary-split\">\n");

    // LEFT PANEL: metrics table
    out.push_str("<div class=\"summary-left\">\n<table>\n");
    out.push_str("<thead><tr><th>Metric</th><th>Value</th></tr></thead>\n<tbody>\n");
    let _ = writeln!(
        out,
        "<tr><td>Total variants</td><td>{}</td></tr>",
        report.total_variants
    );
    let _ = writeln!(
        out,
        "<tr><td>Annotated</td><td>{}</td></tr>",
        report.annotated_variants
    );
    let _ = writeln!(
        out,
        "<tr><td>Scored results</td><td>{}</td></tr>",
        report.results.len()
    );
    let _ = writeln!(
        out,
        "<tr class=\"tier-row tier-row-1\"><td><span class=\"badge badge-tier1\">T1</span> Reliable</td><td>{tier1_count}</td></tr>"
    );
    let _ = writeln!(
        out,
        "<tr class=\"tier-row tier-row-2\"><td><span class=\"badge badge-tier2\">T2</span> Probable</td><td>{tier2_count}</td></tr>"
    );
    let _ = writeln!(
        out,
        "<tr class=\"tier-row tier-row-3\"><td><span class=\"badge badge-tier3\">T3</span> Speculative</td><td>{tier3_count}</td></tr>"
    );
    let _ = writeln!(
        out,
        "<tr><td>Input assembly</td><td>{}</td></tr>",
        html_escape(&report.input_assembly.to_string())
    );
    let _ = writeln!(
        out,
        "<tr><td>Database assembly</td><td>{}</td></tr>",
        html_escape(&report.db_assembly.to_string())
    );
    out.push_str("</tbody>\n</table>\n</div>\n");

    // RIGHT PANEL: contextual details
    out.push_str("<div class=\"summary-right\">\n");

    // Coverage section
    out.push_str("<div class=\"summary-section\">\n");
    out.push_str("<h4>Annotation Coverage</h4>\n");
    let _ = writeln!(
        out,
        "<p>{:.1}% of variants matched at least one database entry \
         (ClinVar, GWAS Catalog, dbSNP/gnomAD, PharmGKB, or SNPedia).</p>",
        coverage_pct
    );
    let _ = writeln!(
        out,
        "<div class=\"coverage-bar-container\"><div class=\"coverage-bar\" style=\"width:{coverage_pct:.1}%\"></div></div>"
    );
    out.push_str("</div>\n");

    // Confidence tier explanations
    out.push_str("<div class=\"summary-section\">\n");
    out.push_str("<h4>Confidence Tiers</h4>\n");
    let _ = writeln!(
        out,
        "<div class=\"tier-explain tier-explain-1\"><strong>Tier 1: Reliable ({tier1_count})</strong>\
         &gt;95% predictive value. Monogenic diseases, carrier status, pharmacogenetics. \
         Sources: ClinVar (\u{2265}2-star review), PharmGKB (Level 1A/1B).</div>"
    );
    let _ = writeln!(
        out,
        "<div class=\"tier-explain tier-explain-2\"><strong>Tier 2: Probable ({tier2_count})</strong>\
         60\u{2013}85% predictive value. Polygenic risk scores, physical traits. \
         Sources: GWAS Catalog (genome-wide significant), SNPedia (magnitude 2\u{2013}3.9).</div>"
    );
    let _ = writeln!(
        out,
        "<div class=\"tier-explain tier-explain-3\"><strong>Tier 3: Speculative ({tier3_count})</strong>\
         50\u{2013}65% predictive value. Complex diseases, personality traits. \
         Sources: GWAS Catalog (low effect size), SNPedia (magnitude &lt;2).</div>"
    );
    out.push_str("</div>\n");

    // Category breakdown with bars
    if !category_counts.is_empty() {
        out.push_str("<div class=\"summary-section\">\n");
        out.push_str("<h4>Category Breakdown</h4>\n");
        for (cat, count) in &category_counts {
            let bar_pct = (*count as f64 / max_cat_count as f64) * 100.0;
            let bar_color = category_bar_color(*cat);
            let _ = writeln!(
                out,
                "<div class=\"cat-bar-row\">\
                 <span class=\"cat-bar-label\">{cat}</span>\
                 <div class=\"cat-bar-bg\"><div class=\"cat-bar\" style=\"width:{bar_pct:.0}%;background:{bar_color}\"></div></div>\
                 <span class=\"cat-bar-count\">{count}</span>\
                 </div>"
            );
        }
        out.push_str("</div>\n");
    }

    out.push_str("</div>\n"); // .summary-right
    out.push_str("</div>\n"); // .summary-split

    // --- Tier 1: Clinically Actionable Findings (split-screen table) ---
    let tier1_results: Vec<&ScoredResult> = report
        .results
        .iter()
        .filter(|r| r.tier == ConfidenceTier::Tier1Reliable)
        .collect();

    if !tier1_results.is_empty() {
        out.push_str("<h3>Clinically Actionable Findings</h3>\n");
        out.push_str("<div class=\"findings-split\">\n");

        // Left panel: compact results table
        out.push_str("<div class=\"findings-list\">\n<table>\n");
        out.push_str("<thead><tr><th>Gene</th><th>rsID</th><th>Category</th><th>Genotype</th></tr></thead>\n<tbody>\n");

        for result in &tier1_results {
            let rsid = result.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");
            let gene = extract_gene(result);
            let gene_display = if gene.is_empty() { "\u{2014}" } else { &gene };
            let genotype = result.variant.variant.genotype.to_string();
            let cat = short_category_label(result.category);

            let _ = writeln!(
                out,
                "<tr><td class=\"gene-name\">{}</td><td class=\"rsid-col\">{}</td><td class=\"cat-label\">{cat}</td><td>{}</td></tr>",
                html_escape(gene_display),
                html_escape(rsid),
                html_escape(&genotype)
            );
        }

        out.push_str("</tbody>\n</table>\n</div>\n");

        // Right panel: full details for each finding
        out.push_str("<div class=\"findings-detail\">\n");
        out.push_str("<h4>Details</h4>\n");

        for result in &tier1_results {
            let rsid = result.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");
            let gene = extract_gene(result);
            let gene_display = if gene.is_empty() {
                "\u{2014}".to_string()
            } else {
                gene
            };

            out.push_str("<div class=\"finding-entry\">\n");

            // Header: Gene + rsID
            let _ = writeln!(
                out,
                "<div class=\"fe-header\"><span class=\"fe-gene\">{}</span><span class=\"fe-rsid\">{}</span>\
                 <span class=\"badge badge-tier1\">Tier 1</span></div>",
                html_escape(&gene_display),
                html_escape(rsid)
            );

            // Summary line
            let _ = writeln!(
                out,
                "<div class=\"fe-summary\">{}</div>",
                html_escape(&result.summary)
            );

            // ClinVar annotation block
            if let Some(ref cv) = result.variant.clinvar {
                out.push_str("<details class=\"anno-block\" open>\n");
                out.push_str("<summary>ClinVar</summary>\n");
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Classification</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    html_escape(&cv.significance)
                );
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Review Status</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    render_stars(cv.review_stars)
                );
                if !cv.conditions.is_empty() {
                    out.push_str("<div class=\"anno-row\"><span class=\"anno-label\">Conditions</span></div>\n");
                    out.push_str("<ul class=\"anno-list\">\n");
                    for cond in &cv.conditions {
                        let _ = writeln!(out, "<li>{}</li>", html_escape(cond));
                    }
                    out.push_str("</ul>\n");
                }
                out.push_str("</details>\n");
            }

            // Pharmacogenomics annotation block
            if let Some(ref pharma) = result.variant.pharmacogenomics {
                out.push_str("<details class=\"anno-block\" open>\n");
                out.push_str("<summary>Pharmacogenomics</summary>\n");
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Gene</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    html_escape(&pharma.gene)
                );
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Drug</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    html_escape(&pharma.drug)
                );
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Evidence Level</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    html_escape(&pharma.evidence_level)
                );
                if let Some(ref pheno) = pharma.phenotype_category {
                    let _ = writeln!(
                        out,
                        "<div class=\"anno-row\"><span class=\"anno-label\">Phenotype</span>\
                         <span class=\"anno-value\">{}</span></div>",
                        html_escape(pheno)
                    );
                }
                if let Some(ref rec) = pharma.clinical_recommendation {
                    let _ = writeln!(
                        out,
                        "<div class=\"anno-row\"><span class=\"anno-label\">Recommendation</span></div>\
                         <div style=\"font-size:0.85rem;margin:0.2rem 0;\">{}</div>",
                        html_escape(rec)
                    );
                }
                out.push_str("</details>\n");
            }

            // Allele frequency block
            if let Some(ref freq) = result.variant.frequency {
                out.push_str("<details class=\"anno-block\">\n");
                out.push_str("<summary>Allele Frequency</summary>\n");
                write_freq_row(out, "Overall", freq.af_total);
                if let Some(af) = freq.af_afr {
                    write_freq_row(out, "African", af);
                }
                if let Some(af) = freq.af_amr {
                    write_freq_row(out, "American", af);
                }
                if let Some(af) = freq.af_eas {
                    write_freq_row(out, "East Asian", af);
                }
                if let Some(af) = freq.af_eur {
                    write_freq_row(out, "European", af);
                }
                if let Some(af) = freq.af_sas {
                    write_freq_row(out, "South Asian", af);
                }
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Source</span>\
                     <span class=\"anno-value\">{}</span></div>",
                    html_escape(&freq.source)
                );
                out.push_str("</details>\n");
            }

            // GWAS hits block
            if !result.variant.gwas_hits.is_empty() {
                out.push_str("<details class=\"anno-block\">\n");
                let _ = writeln!(
                    out,
                    "<summary>GWAS Associations ({})</summary>",
                    result.variant.gwas_hits.len()
                );
                for hit in &result.variant.gwas_hits {
                    let _ = writeln!(
                        out,
                        "<div class=\"anno-row\"><span class=\"anno-label\">Trait</span>\
                         <span class=\"anno-value\">{}</span></div>",
                        html_escape(&hit.trait_name)
                    );
                    let _ = writeln!(
                        out,
                        "<div class=\"anno-row\"><span class=\"anno-label\">p-value</span>\
                         <span class=\"anno-value\">{:.2e}</span></div>",
                        hit.p_value
                    );
                    if let Some(or) = hit.odds_ratio {
                        let _ = writeln!(
                            out,
                            "<div class=\"anno-row\"><span class=\"anno-label\">Odds Ratio</span>\
                             <span class=\"anno-value\">{or:.2}</span></div>"
                        );
                    }
                    if let Some(ref pubmed) = hit.pubmed_id {
                        let _ = writeln!(
                            out,
                            "<div class=\"anno-row\"><span class=\"anno-label\">PubMed</span>\
                             <span class=\"anno-value\">{}</span></div>",
                            html_escape(pubmed)
                        );
                    }
                }
                out.push_str("</details>\n");
            }

            // SNPedia block
            if let Some(ref snp) = result.variant.snpedia {
                out.push_str("<details class=\"anno-block\">\n");
                out.push_str("<summary>SNPedia</summary>\n");
                let _ = writeln!(
                    out,
                    "<div class=\"anno-row\"><span class=\"anno-label\">Magnitude</span>\
                     <span class=\"anno-value\">{:.1}</span></div>",
                    snp.magnitude
                );
                if let Some(ref repute) = snp.repute {
                    let _ = writeln!(
                        out,
                        "<div class=\"anno-row\"><span class=\"anno-label\">Repute</span>\
                         <span class=\"anno-value\">{}</span></div>",
                        html_escape(repute)
                    );
                }
                let _ = writeln!(
                    out,
                    "<div style=\"font-size:0.85rem;margin:0.3rem 0;\">{}</div>",
                    html_escape(&snp.summary)
                );
                out.push_str("</details>\n");
            }

            out.push_str("</div>\n"); // .finding-entry
        }

        out.push_str("</div>\n"); // .findings-detail
        out.push_str("</div>\n"); // .findings-split
    }
}

fn write_assembly_warnings(out: &mut String, warnings: &[String]) {
    if warnings.is_empty() {
        return;
    }

    out.push_str(
        "<div class=\"disclaimer\" style=\"border-left-color:#f59e0b;background:#fffbeb;\">\n",
    );
    out.push_str("<strong>Assembly Warnings</strong>\n");
    for warning in warnings {
        let _ = writeln!(out, "<p>{}</p>", html_escape(warning));
    }
    out.push_str("</div>\n");
}

/// Render ClinVar review stars as HTML.
fn render_stars(count: u8) -> String {
    let filled = "\u{2605}".repeat(count as usize);
    let empty = "\u{2606}".repeat(4_usize.saturating_sub(count as usize));
    format!("<span class=\"stars\">{filled}{empty}</span> ({count}/4)")
}

/// Write a frequency row with a tiny inline bar.
fn write_freq_row(out: &mut String, label: &str, af: f64) {
    let pct = af * 100.0;
    // Scale bar: cap at 50% for visibility (most AFs are small)
    let bar_pct = (pct * 2.0).min(100.0);
    let _ = writeln!(
        out,
        "<div class=\"anno-row\"><span class=\"anno-label\">{label}</span>\
         <span class=\"anno-value\">{pct:.2}%</span></div>\
         <div class=\"freq-bar-bg\"><div class=\"freq-bar\" style=\"width:{bar_pct:.1}%\"></div></div>"
    );
}

/// Short category label for compact table display.
fn short_category_label(cat: ResultCategory) -> &'static str {
    match cat {
        ResultCategory::MonogenicDisease => "Disease",
        ResultCategory::CarrierStatus => "Carrier",
        ResultCategory::Pharmacogenomics => "Pharma",
        ResultCategory::PolygenicRiskScore => "PRS",
        ResultCategory::PhysicalTrait => "Trait",
        ResultCategory::ComplexTrait => "Complex",
        ResultCategory::Ancestry => "Ancestry",
    }
}

/// CSS color for category bar charts (matches TUI category colors).
fn category_bar_color(cat: ResultCategory) -> &'static str {
    match cat {
        ResultCategory::MonogenicDisease => "#dc2626",
        ResultCategory::CarrierStatus => "#2563eb",
        ResultCategory::Pharmacogenomics => "#9333ea",
        ResultCategory::PolygenicRiskScore => "#d97706",
        ResultCategory::PhysicalTrait => "#0891b2",
        ResultCategory::ComplexTrait => "#6b7280",
        ResultCategory::Ancestry => "#374151",
    }
}

fn write_results(out: &mut String, results: &[ScoredResult]) {
    out.push_str("<h2>Results</h2>\n");

    if results.is_empty() {
        out.push_str("<p>No significant findings.</p>\n");
        return;
    }

    let grouped = group_results(results);

    for (tier, categories) in &grouped {
        let _ = writeln!(out, "<h3>{tier}</h3>");

        for (category, items) in categories {
            let _ = writeln!(out, "<h4>{category}</h4>");

            for result in items {
                let rsid = result.variant.variant.rsid.as_deref().unwrap_or("unknown");
                let genotype = result.variant.variant.genotype.to_string();
                let gene = extract_gene(result);
                let badge_class = tier_badge_class(*tier);
                let badge_text = tier_badge_text(*tier);

                out.push_str("<div class=\"result-card\">\n");
                let _ = writeln!(
                    out,
                    "<span class=\"rsid\">{}</span> <span class=\"badge {badge_class}\">{badge_text}</span>",
                    html_escape(rsid)
                );

                let mut meta_parts = Vec::new();
                if !gene.is_empty() {
                    meta_parts.push(format!("Gene: {}", html_escape(&gene)));
                }
                meta_parts.push(format!("Genotype: {}", html_escape(&genotype)));

                let _ = writeln!(out, "<div class=\"meta\">{}</div>", meta_parts.join(" | "));
                let _ = writeln!(
                    out,
                    "<div class=\"summary\">{}</div>",
                    html_escape(&result.summary)
                );
                let _ = writeln!(
                    out,
                    "<div class=\"details\">{}</div>",
                    html_escape(&result.details)
                );
                if !result.limitations.is_empty() {
                    out.push_str(
                        "<div class=\"details\" style=\"color:#b45309;margin-top:0.5rem;\">\n",
                    );
                    out.push_str(
                        "<strong>Limitations:</strong><ul style=\"margin:0.25rem 0 0 1rem;\">\n",
                    );
                    for limitation in &result.limitations {
                        let _ = writeln!(out, "<li>{}</li>", html_escape(limitation));
                    }
                    out.push_str("</ul></div>\n");
                }
                out.push_str("</div>\n");
            }
        }
    }
}

fn write_attributions(out: &mut String, attributions: &[String]) {
    out.push_str("<div class=\"attributions\">\n");
    out.push_str("<h2>Data Sources &amp; Attributions</h2>\n<ul>\n");
    for attr in attributions {
        let _ = writeln!(out, "<li>{}</li>", html_escape(attr));
    }
    out.push_str("</ul>\n</div>\n");
}

fn write_body_close(out: &mut String) {
    out.push_str("<div class=\"footer\">Generated by GeneSight</div>\n");
    out.push_str("</body>\n</html>\n");
}

/// Group results by tier then category for display.
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

fn tier_badge_class(tier: ConfidenceTier) -> &'static str {
    match tier {
        ConfidenceTier::Tier1Reliable => "badge-tier1",
        ConfidenceTier::Tier2Probable => "badge-tier2",
        ConfidenceTier::Tier3Speculative => "badge-tier3",
    }
}

fn tier_badge_text(tier: ConfidenceTier) -> &'static str {
    match tier {
        ConfidenceTier::Tier1Reliable => "Tier 1: Reliable",
        ConfidenceTier::Tier2Probable => "Tier 2: Probable",
        ConfidenceTier::Tier3Speculative => "Tier 3: Speculative",
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

/// Escape HTML special characters to prevent XSS and rendering issues.
fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::annotation::*;
    use crate::models::assembly::GenomeAssembly;
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
            }),
            snpedia: None,
            gwas_hits: Vec::new(),
            frequency: None,
            pharmacogenomics: None,
        };

        Report {
            total_variants: 600000,
            annotated_variants: 1500,
            results: vec![ScoredResult {
                variant: annotated,
                tier: ConfidenceTier::Tier1Reliable,
                category: ResultCategory::MonogenicDisease,
                summary: "BRCA1 (rs123) — Pathogenic (3-star review)".to_string(),
                details: "Genotype: AG. Classification: Pathogenic.".to_string(),
                limitations: Vec::new(),
            }],
            attributions: vec!["ClinVar: NCBI/NLM (public domain)".to_string()],
            disclaimer: "This is not medical advice.".to_string(),
            input_assembly: GenomeAssembly::GRCh37,
            db_assembly: GenomeAssembly::GRCh37,
            assembly_warnings: Vec::new(),
        }
    }

    #[test]
    fn render_produces_valid_html_structure() {
        let report = make_test_report();
        let html = render(&report).expect("render");
        assert!(html.starts_with("<!DOCTYPE html>"));
        assert!(html.contains("</html>"));
        assert!(html.contains("<title>GeneSight Analysis Report</title>"));
    }

    #[test]
    fn render_contains_disclaimer() {
        let report = make_test_report();
        let html = render(&report).expect("render");
        assert!(html.contains("Medical Disclaimer"));
        assert!(html.contains("not medical advice"));
    }

    #[test]
    fn render_contains_results() {
        let report = make_test_report();
        let html = render(&report).expect("render");
        assert!(html.contains("rs123"));
        assert!(html.contains("BRCA1"));
        assert!(html.contains("Tier 1: Reliable"));
    }

    #[test]
    fn render_escapes_html() {
        let escaped = html_escape("<script>alert('xss')</script>");
        assert!(!escaped.contains('<'));
        assert!(!escaped.contains('>'));
        assert!(escaped.contains("&lt;"));
        assert!(escaped.contains("&gt;"));
    }

    #[test]
    fn render_empty_results() {
        let report = Report {
            total_variants: 0,
            annotated_variants: 0,
            results: vec![],
            attributions: vec![],
            disclaimer: "Disclaimer.".to_string(),
            input_assembly: GenomeAssembly::Unknown,
            db_assembly: GenomeAssembly::Unknown,
            assembly_warnings: Vec::new(),
        };
        let html = render(&report).expect("render");
        assert!(html.contains("No significant findings"));
    }

    #[test]
    fn render_contains_assembly_info() {
        let report = make_test_report();
        let html = render(&report).expect("render");
        assert!(html.contains("Input assembly"));
        assert!(html.contains("Database assembly"));
        assert!(html.contains("GRCh37 (hg19)"));
    }

    #[test]
    fn render_contains_assembly_warnings() {
        let mut report = make_test_report();
        report.assembly_warnings = vec!["Assembly mismatch detected".to_string()];
        let html = render(&report).expect("render");
        assert!(html.contains("Assembly Warnings"));
        assert!(html.contains("Assembly mismatch detected"));
    }
}
