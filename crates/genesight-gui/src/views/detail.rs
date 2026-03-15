use egui::RichText;
use genesight_core::models::ScoredResult;

use crate::state::{gene_name, short_category};
use crate::theme;

/// Draw the detail panel for a selected result.
pub fn draw(ui: &mut egui::Ui, result: &ScoredResult) {
    egui::ScrollArea::vertical()
        .auto_shrink(false)
        .show(ui, |ui| {
            let gene = gene_name(result);
            let rsid = result.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");

            // ── Header ───────────────────────────────────────────────
            ui.vertical(|ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new(&gene)
                            .size(20.0)
                            .strong()
                            .color(theme::TEXT_PRIMARY),
                    );
                    tier_badge(ui, result.tier);
                });
                ui.add_space(2.0);
                ui.label(
                    RichText::new(rsid)
                        .size(13.0)
                        .color(theme::TEXT_SECONDARY)
                        .family(egui::FontFamily::Monospace),
                );
            });

            ui.add_space(12.0);

            // ── Category + Location bar ──────────────────────────────
            ui.horizontal_wrapped(|ui| {
                ui.spacing_mut().item_spacing.x = 6.0;
                // Category badge
                let cat = short_category(result.category);
                let cat_color = theme::category_color(result.category);
                let cat_bg = theme::category_bg(result.category);
                egui::Frame::NONE
                    .fill(cat_bg)
                    .corner_radius(egui::CornerRadius::same(4))
                    .inner_margin(egui::Margin::symmetric(8, 3))
                    .show(ui, |ui| {
                        ui.label(RichText::new(cat).size(11.0).strong().color(cat_color));
                    });

                // Location badge
                egui::Frame::NONE
                    .fill(theme::BG_SIDEBAR)
                    .corner_radius(egui::CornerRadius::same(4))
                    .inner_margin(egui::Margin::symmetric(8, 3))
                    .show(ui, |ui| {
                        ui.label(
                            RichText::new(format!(
                                "chr{}:{}",
                                result.variant.variant.chromosome, result.variant.variant.position
                            ))
                            .size(11.0)
                            .family(egui::FontFamily::Monospace)
                            .color(theme::TEXT_SECONDARY),
                        );
                    });

                // Genotype badge
                egui::Frame::NONE
                    .fill(theme::BG_SIDEBAR)
                    .corner_radius(egui::CornerRadius::same(4))
                    .inner_margin(egui::Margin::symmetric(8, 3))
                    .show(ui, |ui| {
                        ui.label(
                            RichText::new(result.variant.variant.genotype.to_string())
                                .size(11.0)
                                .family(egui::FontFamily::Monospace)
                                .strong()
                                .color(theme::TEXT_PRIMARY),
                        );
                    });
            });

            ui.add_space(12.0);

            // ── Summary ──────────────────────────────────────────────
            theme::card_frame()
                .inner_margin(egui::Margin::same(16))
                .shadow(egui::Shadow::NONE)
                .show(ui, |ui| {
                    ui.label(
                        RichText::new(&result.summary)
                            .size(13.0)
                            .color(theme::TEXT_PRIMARY),
                    );
                    if !result.details.is_empty() {
                        ui.add_space(6.0);
                        ui.label(
                            RichText::new(&result.details)
                                .size(12.0)
                                .color(theme::TEXT_SECONDARY),
                        );
                    }
                });

            ui.add_space(8.0);

            // ── ClinVar ──────────────────────────────────────────────
            if let Some(ref cv) = result.variant.clinvar {
                collapsible_section(ui, "\u{1F3E5}  ClinVar", "clinvar", true, |ui| {
                    info_row(ui, "Classification", &cv.significance);
                    info_row(ui, "Review Status", &format!("{}/4 stars", cv.review_stars));
                    if let Some(ref gene) = cv.gene_symbol {
                        info_row(ui, "Gene", gene);
                    }
                    if !cv.conditions.is_empty() {
                        ui.add_space(6.0);
                        ui.label(
                            RichText::new("Conditions")
                                .size(11.0)
                                .strong()
                                .color(theme::TEXT_SECONDARY),
                        );
                        for cond in &cv.conditions {
                            ui.horizontal(|ui| {
                                ui.add_space(8.0);
                                ui.label(
                                    RichText::new("\u{2022}")
                                        .size(12.0)
                                        .color(theme::TEXT_MUTED),
                                );
                                ui.label(RichText::new(cond).size(12.0).color(theme::TEXT_PRIMARY));
                            });
                        }
                    }
                });
            }

            // ── Pharmacogenomics ─────────────────────────────────────
            if let Some(ref pharma) = result.variant.pharmacogenomics {
                collapsible_section(ui, "\u{1F48A}  Pharmacogenomics", "pharma", true, |ui| {
                    info_row(ui, "Gene", &pharma.gene);
                    info_row(ui, "Drug", &pharma.drug);
                    info_row(ui, "Evidence", &pharma.evidence_level);
                    if let Some(ref pheno) = pharma.phenotype_category {
                        info_row(ui, "Phenotype", pheno);
                    }
                    if let Some(ref rec) = pharma.clinical_recommendation {
                        ui.add_space(6.0);
                        egui::Frame::NONE
                            .fill(theme::ACCENT_LIGHT)
                            .corner_radius(egui::CornerRadius::same(6))
                            .inner_margin(egui::Margin::same(10))
                            .show(ui, |ui| {
                                ui.label(
                                    RichText::new("Recommendation")
                                        .size(10.0)
                                        .strong()
                                        .color(theme::ACCENT),
                                );
                                ui.label(RichText::new(rec).size(12.0).color(theme::TEXT_PRIMARY));
                            });
                    }
                });
            }

            // ── Allele Frequency ─────────────────────────────────────
            if let Some(ref freq) = result.variant.frequency {
                collapsible_section(ui, "\u{1F4CA}  Allele Frequency", "freq", false, |ui| {
                    // Overall frequency with visual bar
                    let pct = freq.af_total * 100.0;
                    ui.horizontal(|ui| {
                        ui.label(
                            RichText::new("Overall")
                                .size(11.0)
                                .color(theme::TEXT_SECONDARY),
                        );
                        ui.label(
                            RichText::new(format!("{pct:.4}%"))
                                .size(12.0)
                                .strong()
                                .color(theme::TEXT_PRIMARY),
                        );
                    });
                    // Mini frequency bar
                    let bar_width = ui.available_width().min(200.0);
                    let bar_frac = (freq.af_total as f32).clamp(0.0, 1.0);
                    ui.add(
                        egui::ProgressBar::new(bar_frac)
                            .desired_width(bar_width)
                            .corner_radius(egui::CornerRadius::same(3)),
                    );
                    ui.add_space(6.0);

                    // Per-population
                    let pops: &[(&str, Option<f64>)] = &[
                        ("African", freq.af_afr),
                        ("American", freq.af_amr),
                        ("East Asian", freq.af_eas),
                        ("European", freq.af_eur),
                        ("South Asian", freq.af_sas),
                    ];
                    for (name, af) in pops {
                        if let Some(af) = af {
                            freq_row(ui, name, *af);
                        }
                    }
                    ui.add_space(4.0);
                    info_row(ui, "Source", &freq.source);
                });
            }

            // ── GWAS ─────────────────────────────────────────────────
            if !result.variant.gwas_hits.is_empty() {
                let label = format!(
                    "\u{1F52C}  GWAS Associations ({})",
                    result.variant.gwas_hits.len()
                );
                collapsible_section(ui, &label, "gwas", false, |ui| {
                    for (i, hit) in result.variant.gwas_hits.iter().enumerate() {
                        if i > 0 {
                            ui.add_space(4.0);
                            ui.separator();
                            ui.add_space(4.0);
                        }
                        ui.label(
                            RichText::new(&hit.trait_name)
                                .size(12.0)
                                .strong()
                                .color(theme::TEXT_PRIMARY),
                        );
                        ui.add_space(2.0);
                        egui::Grid::new(format!("gwas_{i}"))
                            .num_columns(2)
                            .spacing([12.0, 3.0])
                            .show(ui, |ui| {
                                grid_row(ui, "p-value", &format!("{:.2e}", hit.p_value));
                                if let Some(or) = hit.odds_ratio {
                                    grid_row(ui, "Odds Ratio", &format!("{or:.2}"));
                                }
                                if let Some(ref allele) = hit.risk_allele {
                                    grid_row(ui, "Risk Allele", allele);
                                }
                                if let Some(ref pmid) = hit.pubmed_id {
                                    grid_row(ui, "PubMed", pmid);
                                }
                            });
                    }
                });
            }

            // ── SNPedia ──────────────────────────────────────────────
            if let Some(ref snp) = result.variant.snpedia {
                collapsible_section(ui, "\u{1F4D6}  SNPedia", "snpedia", false, |ui| {
                    ui.horizontal(|ui| {
                        // Magnitude badge
                        let mag_color = if snp.magnitude >= 3.0 {
                            theme::DANGER
                        } else if snp.magnitude >= 2.0 {
                            theme::WARNING
                        } else {
                            theme::TEXT_SECONDARY
                        };
                        let mag_bg = if snp.magnitude >= 3.0 {
                            theme::DANGER_BG
                        } else if snp.magnitude >= 2.0 {
                            theme::WARNING_BG
                        } else {
                            theme::BG_SIDEBAR
                        };
                        egui::Frame::NONE
                            .fill(mag_bg)
                            .corner_radius(egui::CornerRadius::same(4))
                            .inner_margin(egui::Margin::symmetric(8, 3))
                            .show(ui, |ui| {
                                ui.label(
                                    RichText::new(format!("Mag {:.1}", snp.magnitude))
                                        .size(11.0)
                                        .strong()
                                        .color(mag_color),
                                );
                            });

                        if let Some(ref rep) = snp.repute {
                            egui::Frame::NONE
                                .fill(theme::BG_SIDEBAR)
                                .corner_radius(egui::CornerRadius::same(4))
                                .inner_margin(egui::Margin::symmetric(8, 3))
                                .show(ui, |ui| {
                                    ui.label(
                                        RichText::new(rep).size(11.0).color(theme::TEXT_SECONDARY),
                                    );
                                });
                        }
                    });
                    ui.add_space(6.0);
                    ui.label(
                        RichText::new(&snp.summary)
                            .size(12.0)
                            .color(theme::TEXT_PRIMARY),
                    );
                });
            }

            ui.add_space(20.0);
        });
}

// ── Helper widgets ───────────────────────────────────────────────────────

fn info_row(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(11.0).color(theme::TEXT_MUTED));
        ui.label(
            RichText::new(value)
                .size(12.0)
                .strong()
                .color(theme::TEXT_PRIMARY),
        );
    });
}

fn grid_row(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.label(RichText::new(label).size(11.0).color(theme::TEXT_SECONDARY));
    ui.label(
        RichText::new(value)
            .size(12.0)
            .strong()
            .color(theme::TEXT_PRIMARY),
    );
    ui.end_row();
}

fn freq_row(ui: &mut egui::Ui, label: &str, af: f64) {
    ui.horizontal(|ui| {
        ui.add_space(8.0);
        ui.label(RichText::new(label).size(11.0).color(theme::TEXT_SECONDARY));
        ui.label(
            RichText::new(format!("{:.4}%", af * 100.0))
                .size(11.0)
                .color(theme::TEXT_PRIMARY),
        );
    });
}

fn tier_badge(ui: &mut egui::Ui, tier: genesight_core::models::ConfidenceTier) {
    let color = theme::tier_color(tier);
    let bg = theme::tier_bg(tier);
    let label = theme::tier_label(tier);
    egui::Frame::NONE
        .fill(bg)
        .corner_radius(egui::CornerRadius::same(12))
        .inner_margin(egui::Margin::symmetric(10, 3))
        .show(ui, |ui| {
            ui.label(RichText::new(label).size(11.0).strong().color(color));
        });
}

fn collapsible_section(
    ui: &mut egui::Ui,
    title: &str,
    id: &str,
    default_open: bool,
    content: impl FnOnce(&mut egui::Ui),
) {
    ui.add_space(4.0);
    theme::section_frame().show(ui, |ui| {
        egui::CollapsingHeader::new(
            RichText::new(title)
                .size(12.0)
                .strong()
                .color(theme::ACCENT),
        )
        .id_salt(id)
        .default_open(default_open)
        .show(ui, |ui| {
            ui.add_space(4.0);
            content(ui);
        });
    });
}
