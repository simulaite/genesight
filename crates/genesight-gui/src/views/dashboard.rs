use egui::{Color32, RichText, Vec2};
use genesight_core::models::ConfidenceTier;

use crate::state::{gene_name, short_category, AppData};
use crate::theme;
use crate::views::detail;

/// Draw the main dashboard view.
pub fn draw(ui: &mut egui::Ui, data: &mut AppData) {
    if data.report.is_none() {
        return;
    }

    // Stats bar at top
    draw_stats_bar(ui, data);

    // Main split: sidebar | results table | detail panel
    let available = ui.available_size();

    // Dynamic sidebar width: 15% of window, clamped between 180–260px
    let sidebar_width = (available.x * 0.15).clamp(180.0, 260.0);

    egui::SidePanel::left("sidebar")
        .resizable(true)
        .default_width(sidebar_width)
        .width_range(160.0..=320.0)
        .show_inside(ui, |ui| {
            draw_sidebar(ui, data);
        });

    // Right detail panel: 30-40% of window, responsive
    let detail_result = data.selected_result().cloned();
    let detail_width = (available.x * 0.33).clamp(280.0, 500.0);

    egui::SidePanel::right("detail_panel")
        .resizable(true)
        .default_width(detail_width)
        .min_width(240.0)
        .show_inside(ui, |ui| {
            egui::Frame::NONE
                .fill(theme::BG_SURFACE)
                .inner_margin(egui::Margin::same(16))
                .show(ui, |ui| {
                    if let Some(ref result) = detail_result {
                        detail::draw(ui, result);
                    } else {
                        ui.with_layout(
                            egui::Layout::centered_and_justified(egui::Direction::TopDown),
                            |ui| {
                                ui.vertical_centered(|ui| {
                                    ui.add_space(ui.available_height() * 0.3);
                                    ui.label(RichText::new("\u{1F50D}").size(32.0));
                                    ui.add_space(8.0);
                                    ui.label(
                                        RichText::new("Select a result")
                                            .size(14.0)
                                            .color(theme::TEXT_MUTED),
                                    );
                                    ui.label(
                                        RichText::new("Click a row to view details")
                                            .size(12.0)
                                            .color(theme::TEXT_MUTED),
                                    );
                                });
                            },
                        );
                    }
                });
        });

    // Center: results table
    egui::CentralPanel::default().show_inside(ui, |ui| {
        egui::Frame::NONE
            .inner_margin(egui::Margin::symmetric(8, 4))
            .show(ui, |ui| {
                draw_results_table(ui, data);
            });
    });

    // Keyboard navigation
    if ui.input(|i| i.key_pressed(egui::Key::ArrowDown)) && !data.filtered_indices.is_empty() {
        data.selected_index = (data.selected_index + 1).min(data.filtered_indices.len() - 1);
    }
    if ui.input(|i| i.key_pressed(egui::Key::ArrowUp)) && !data.filtered_indices.is_empty() {
        data.selected_index = data.selected_index.saturating_sub(1);
    }
}

fn draw_stats_bar(ui: &mut egui::Ui, data: &AppData) {
    let report = data.report.as_ref().unwrap();
    let (t1, t2, t3) = data.tier_counts();

    egui::Frame::NONE
        .fill(theme::BG_SURFACE)
        .stroke(egui::Stroke::new(1.0, theme::BORDER_LIGHT))
        .inner_margin(egui::Margin::symmetric(16, 10))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 6.0;

                stat_chip(
                    ui,
                    "Variants",
                    &report.total_variants.to_string(),
                    theme::TEXT_PRIMARY,
                    theme::BG_SIDEBAR,
                );
                stat_chip(
                    ui,
                    "Annotated",
                    &report.annotated_variants.to_string(),
                    theme::TEXT_PRIMARY,
                    theme::BG_SIDEBAR,
                );
                stat_chip(
                    ui,
                    "Results",
                    &report.results.len().to_string(),
                    theme::ACCENT,
                    theme::ACCENT_LIGHT,
                );

                ui.add_space(4.0);
                ui.separator();
                ui.add_space(4.0);

                stat_chip(ui, "T1", &t1.to_string(), theme::TIER1, theme::TIER1_BG);
                stat_chip(ui, "T2", &t2.to_string(), theme::TIER2, theme::TIER2_BG);
                stat_chip(ui, "T3", &t3.to_string(), theme::TIER3, theme::TIER3_BG);

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    egui::Frame::NONE
                        .fill(theme::BG_SIDEBAR)
                        .corner_radius(egui::CornerRadius::same(12))
                        .inner_margin(egui::Margin::symmetric(10, 4))
                        .show(ui, |ui| {
                            ui.label(
                                RichText::new(format!("{} shown", data.filtered_indices.len()))
                                    .size(11.0)
                                    .color(theme::TEXT_SECONDARY),
                            );
                        });
                });
            });
        });
}

fn stat_chip(ui: &mut egui::Ui, label: &str, value: &str, color: Color32, bg: Color32) {
    egui::Frame::NONE
        .fill(bg)
        .corner_radius(egui::CornerRadius::same(8))
        .inner_margin(egui::Margin::symmetric(10, 5))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 4.0;
                ui.label(RichText::new(label).size(11.0).color(theme::TEXT_MUTED));
                ui.label(RichText::new(value).size(13.0).strong().color(color));
            });
        });
}

fn draw_sidebar(ui: &mut egui::Ui, data: &mut AppData) {
    egui::Frame::NONE.fill(theme::BG_SIDEBAR).show(ui, |ui| {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add_space(12.0);

            // File info
            if let Some(ref path) = data.dna_file {
                let name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default();

                egui::Frame::NONE
                    .fill(theme::BG_SURFACE)
                    .corner_radius(egui::CornerRadius::same(8))
                    .inner_margin(egui::Margin::same(10))
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(RichText::new("\u{1F4C4}").size(14.0));
                            ui.vertical(|ui| {
                                ui.spacing_mut().item_spacing.y = 1.0;
                                ui.label(
                                    RichText::new("Loaded File")
                                        .size(10.0)
                                        .color(theme::TEXT_MUTED),
                                );
                                ui.label(
                                    RichText::new(&name)
                                        .size(12.0)
                                        .strong()
                                        .color(theme::TEXT_PRIMARY),
                                );
                            });
                        });
                    });
                ui.add_space(12.0);
            }

            // Search
            sidebar_heading(ui, "Search");
            ui.add_space(4.0);
            let search_changed = ui
                .add(
                    egui::TextEdit::singleline(&mut data.search_query)
                        .hint_text("\u{1F50D} rsID, gene, keyword...")
                        .desired_width(ui.available_width() - 16.0)
                        .margin(egui::Margin::symmetric(8, 6)),
                )
                .changed();

            ui.add_space(16.0);

            // Tier filter
            sidebar_heading(ui, "Confidence");
            ui.add_space(6.0);

            let (t1, t2, t3) = data.tier_counts();
            let mut filter_changed = false;

            filter_changed |= tier_checkbox(
                ui,
                &mut data.tier_filter[0],
                "Reliable",
                t1,
                theme::TIER1,
                theme::TIER1_BG,
            );
            filter_changed |= tier_checkbox(
                ui,
                &mut data.tier_filter[1],
                "Probable",
                t2,
                theme::TIER2,
                theme::TIER2_BG,
            );
            filter_changed |= tier_checkbox(
                ui,
                &mut data.tier_filter[2],
                "Speculative",
                t3,
                theme::TIER3,
                theme::TIER3_BG,
            );

            ui.add_space(16.0);

            // Dataset toggles
            sidebar_heading(ui, "Datasets");
            ui.add_space(6.0);

            egui::Frame::NONE
                .fill(theme::BG_SURFACE)
                .corner_radius(egui::CornerRadius::same(8))
                .inner_margin(egui::Margin::same(10))
                .show(ui, |ui| {
                    ui.spacing_mut().item_spacing.y = 4.0;
                    ui.checkbox(&mut data.annotation_config.clinvar, "ClinVar");
                    ui.checkbox(&mut data.annotation_config.gwas, "GWAS Catalog");
                    ui.checkbox(&mut data.annotation_config.frequencies, "Allele Freq.");
                    ui.checkbox(&mut data.annotation_config.pharmacogenomics, "PharmGKB");
                    ui.checkbox(&mut data.snpedia_enabled, "SNPedia");
                });

            ui.add_space(6.0);
            ui.horizontal_wrapped(|ui| {
                ui.add_space(4.0);
                ui.label(
                    RichText::new("Re-analyze to apply dataset changes")
                        .size(10.0)
                        .color(theme::TEXT_MUTED),
                );
            });

            if search_changed || filter_changed {
                data.apply_filter();
            }

            ui.add_space(16.0);
        });
    });
}

fn sidebar_heading(ui: &mut egui::Ui, text: &str) {
    ui.horizontal(|ui| {
        ui.add_space(4.0);
        ui.label(
            RichText::new(text.to_uppercase())
                .size(10.0)
                .strong()
                .color(theme::TEXT_MUTED),
        );
    });
}

fn tier_checkbox(
    ui: &mut egui::Ui,
    value: &mut bool,
    label: &str,
    count: usize,
    color: Color32,
    bg: Color32,
) -> bool {
    let mut changed = false;

    egui::Frame::NONE
        .fill(theme::BG_SURFACE)
        .corner_radius(egui::CornerRadius::same(6))
        .inner_margin(egui::Margin::symmetric(8, 4))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                changed = ui.checkbox(value, "").changed();

                // Colored dot
                egui::Frame::NONE
                    .fill(bg)
                    .corner_radius(egui::CornerRadius::same(4))
                    .inner_margin(egui::Margin::symmetric(4, 2))
                    .show(ui, |ui| {
                        ui.label(RichText::new("\u{25CF}").size(8.0).color(color));
                    });

                ui.label(RichText::new(label).size(12.0).color(theme::TEXT_PRIMARY));

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    egui::Frame::NONE
                        .fill(bg)
                        .corner_radius(egui::CornerRadius::same(8))
                        .inner_margin(egui::Margin::symmetric(6, 2))
                        .show(ui, |ui| {
                            ui.label(
                                RichText::new(format!("{count}"))
                                    .size(11.0)
                                    .strong()
                                    .color(color),
                            );
                        });
                });
            });
        });

    changed
}

fn draw_results_table(ui: &mut egui::Ui, data: &mut AppData) {
    let report = match &data.report {
        Some(r) => r,
        None => return,
    };

    if data.filtered_indices.is_empty() {
        ui.with_layout(
            egui::Layout::centered_and_justified(egui::Direction::TopDown),
            |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(ui.available_height() * 0.3);
                    ui.label(
                        RichText::new("No results match the current filters")
                            .color(theme::TEXT_MUTED)
                            .size(14.0),
                    );
                    ui.add_space(4.0);
                    ui.label(
                        RichText::new("Try adjusting your search or confidence filters")
                            .color(theme::TEXT_MUTED)
                            .size(12.0),
                    );
                });
            },
        );
        return;
    }

    let row_height = 32.0;
    let num_rows = data.filtered_indices.len();

    // Header row
    egui::Frame::NONE
        .fill(theme::BG_SIDEBAR)
        .inner_margin(egui::Margin::symmetric(4, 6))
        .corner_radius(egui::CornerRadius {
            nw: 6,
            ne: 6,
            sw: 0,
            se: 0,
        })
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                let widths = compute_column_widths(ui.available_width());
                ui.add_space(4.0);
                header_cell(ui, "rsID", widths[0]);
                header_cell(ui, "Gene", widths[1]);
                header_cell(ui, "Category", widths[2]);
                header_cell(ui, "Summary", widths[3]);
                header_cell(ui, "Tier", widths[4]);
            });
        });

    // Virtual scrolling table body
    egui::ScrollArea::vertical().auto_shrink(false).show_rows(
        ui,
        row_height,
        num_rows,
        |ui, row_range| {
            for row_idx in row_range {
                let result_idx = data.filtered_indices[row_idx];
                let result = &report.results[result_idx];
                let is_selected = row_idx == data.selected_index;

                let bg = if is_selected {
                    theme::SELECTED_ROW
                } else if row_idx % 2 == 0 {
                    theme::BG_SURFACE
                } else {
                    theme::BG_PRIMARY
                };

                let response = egui::Frame::NONE
                    .fill(bg)
                    .inner_margin(egui::Margin::symmetric(4, 4))
                    .show(ui, |ui| {
                        ui.set_min_height(row_height);
                        ui.horizontal(|ui| {
                            let widths = compute_column_widths(ui.available_width());

                            // rsID
                            let rsid = result.variant.variant.rsid.as_deref().unwrap_or("\u{2014}");
                            table_cell(
                                ui,
                                rsid,
                                widths[0],
                                egui::FontFamily::Monospace,
                                theme::TEXT_PRIMARY,
                            );

                            // Gene
                            let gene = gene_name(result);
                            table_cell(
                                ui,
                                &gene,
                                widths[1],
                                egui::FontFamily::Proportional,
                                theme::TEXT_PRIMARY,
                            );

                            // Category badge
                            let cat = short_category(result.category);
                            let cat_color = theme::category_color(result.category);
                            let cat_bg = theme::category_bg(result.category);
                            ui.allocate_ui(Vec2::new(widths[2], 24.0), |ui| {
                                egui::Frame::NONE
                                    .fill(cat_bg)
                                    .corner_radius(egui::CornerRadius::same(4))
                                    .inner_margin(egui::Margin::symmetric(6, 2))
                                    .show(ui, |ui| {
                                        ui.label(RichText::new(cat).size(11.0).color(cat_color));
                                    });
                            });

                            // Summary (truncated)
                            let summary = truncate(&result.summary, 50);
                            table_cell(
                                ui,
                                &summary,
                                widths[3],
                                egui::FontFamily::Proportional,
                                theme::TEXT_SECONDARY,
                            );

                            // Tier badge
                            let tier_color = theme::tier_color(result.tier);
                            let tier_bg = theme::tier_bg(result.tier);
                            let tier_text = match result.tier {
                                ConfidenceTier::Tier1Reliable => "T1",
                                ConfidenceTier::Tier2Probable => "T2",
                                ConfidenceTier::Tier3Speculative => "T3",
                            };
                            ui.allocate_ui(Vec2::new(widths[4], 24.0), |ui| {
                                egui::Frame::NONE
                                    .fill(tier_bg)
                                    .corner_radius(egui::CornerRadius::same(10))
                                    .inner_margin(egui::Margin::symmetric(8, 2))
                                    .show(ui, |ui| {
                                        ui.label(
                                            RichText::new(tier_text)
                                                .size(11.0)
                                                .strong()
                                                .color(tier_color),
                                        );
                                    });
                            });
                        });
                    })
                    .response;

                let response = response.interact(egui::Sense::click());
                if response.clicked() {
                    data.selected_index = row_idx;
                }
                if response.hovered() {
                    ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
                }
            }
        },
    );
}

fn compute_column_widths(total: f32) -> [f32; 5] {
    let total = (total - 24.0).max(200.0);
    [
        total * 0.13, // rsID
        total * 0.13, // Gene
        total * 0.12, // Category
        total * 0.50, // Summary
        total * 0.08, // Tier
    ]
}

fn header_cell(ui: &mut egui::Ui, text: &str, width: f32) {
    ui.allocate_ui(Vec2::new(width, 20.0), |ui| {
        ui.label(
            RichText::new(text)
                .size(11.0)
                .strong()
                .color(theme::TEXT_MUTED),
        );
    });
}

fn table_cell(ui: &mut egui::Ui, text: &str, width: f32, family: egui::FontFamily, color: Color32) {
    ui.allocate_ui(Vec2::new(width, 24.0), |ui| {
        ui.label(RichText::new(text).size(12.0).family(family).color(color));
    });
}

fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len.saturating_sub(3)])
    }
}
