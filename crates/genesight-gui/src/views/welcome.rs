use egui::{Align, Color32, Layout, RichText, Vec2};

use crate::state::AppData;
use crate::theme;

/// Draw the welcome/file picker screen.
pub fn draw(ui: &mut egui::Ui, data: &AppData) -> Option<std::path::PathBuf> {
    let mut selected_file = None;
    let available_width = ui.available_width();

    egui::ScrollArea::vertical()
        .auto_shrink(false)
        .show(ui, |ui| {
            ui.vertical_centered(|ui| {
                // Dynamic top spacing based on window height
                let top_space = (ui.available_height() * 0.08).clamp(20.0, 80.0);
                ui.add_space(top_space);

                // App icon / DNA helix visual
                ui.label(
                    RichText::new("\u{1F9EC}")
                        .size(48.0),
                );
                ui.add_space(12.0);

                // Title
                ui.label(
                    RichText::new("GeneSight")
                        .size(38.0)
                        .strong()
                        .color(theme::TEXT_PRIMARY),
                );
                ui.add_space(4.0);
                ui.label(
                    RichText::new("Privacy-First DNA Analysis")
                        .size(16.0)
                        .color(theme::TEXT_SECONDARY),
                );

                ui.add_space(32.0);

                // Card width scales with window
                let card_width = available_width.clamp(320.0, 560.0) - 40.0;

                // ── File picker card ──────────────────────────────────────
                ui.allocate_ui_with_layout(
                    Vec2::new(card_width, 0.0),
                    Layout::top_down(Align::Center),
                    |ui| {
                        theme::card_frame()
                            .inner_margin(egui::Margin::same(32))
                            .show(ui, |ui| {
                                ui.vertical_centered(|ui| {
                                    ui.label(
                                        RichText::new("Get Started")
                                            .size(20.0)
                                            .strong()
                                            .color(theme::TEXT_PRIMARY),
                                    );
                                    ui.add_space(4.0);
                                    ui.label(
                                        RichText::new(
                                            "Select a raw DNA data file to begin analysis",
                                        )
                                        .size(14.0)
                                        .color(theme::TEXT_SECONDARY),
                                    );

                                    ui.add_space(24.0);

                                    // Open file button
                                    let btn = egui::Button::new(
                                        RichText::new("Choose File...")
                                            .size(15.0)
                                            .strong()
                                            .color(Color32::WHITE),
                                    )
                                    .fill(theme::ACCENT)
                                    .stroke(egui::Stroke::NONE)
                                    .corner_radius(egui::CornerRadius::same(8))
                                    .min_size(Vec2::new(220.0, 44.0));

                                    let response = ui.add(btn);
                                    if response.clicked() {
                                        if let Some(path) = rfd::FileDialog::new()
                                            .set_title("Open DNA Raw Data File")
                                            .add_filter(
                                                "DNA Files",
                                                &["txt", "csv", "tsv", "vcf", "gz"],
                                            )
                                            .add_filter("All Files", &["*"])
                                            .pick_file()
                                        {
                                            selected_file = Some(path);
                                        }
                                    }

                                    // Hover feedback on button
                                    if response.hovered() {
                                        ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
                                    }

                                    ui.add_space(20.0);

                                    // Supported formats as pills
                                    ui.label(
                                        RichText::new("Supported formats")
                                            .size(11.0)
                                            .color(theme::TEXT_MUTED),
                                    );
                                    ui.add_space(6.0);
                                    ui.horizontal_wrapped(|ui| {
                                        ui.spacing_mut().item_spacing.x = 6.0;
                                        for fmt in &["23andMe", "AncestryDNA", "VCF"] {
                                            format_pill(ui, fmt);
                                        }
                                    });
                                });
                            });
                    },
                );

                ui.add_space(16.0);

                // ── Feature highlights ────────────────────────────────────
                // Show on wider screens only
                if card_width > 400.0 {
                    ui.allocate_ui_with_layout(
                        Vec2::new(card_width, 0.0),
                        Layout::top_down(Align::Center),
                        |ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.spacing_mut().item_spacing = egui::vec2(8.0, 8.0);
                                let pill_w = ((card_width - 24.0) / 3.0).max(100.0);
                                feature_chip(ui, "\u{1F512}", "100% Local", "No data uploaded", pill_w);
                                feature_chip(ui, "\u{1F4CA}", "5 Databases", "ClinVar, GWAS & more", pill_w);
                                feature_chip(ui, "\u{1F4CB}", "Export", "HTML & PDF reports", pill_w);
                            });
                        },
                    );
                    ui.add_space(16.0);
                }

                // ── Database status card ──────────────────────────────────
                ui.allocate_ui_with_layout(
                    Vec2::new(card_width, 0.0),
                    Layout::top_down(Align::LEFT),
                    |ui| {
                        theme::card_frame()
                            .inner_margin(egui::Margin::same(20))
                            .show(ui, |ui| {
                                ui.label(
                                    RichText::new("Database Status")
                                        .size(13.0)
                                        .strong()
                                        .color(theme::TEXT_PRIMARY),
                                );
                                ui.add_space(10.0);

                                let db_ok = data.main_db_ok();
                                let snp_ok = data
                                    .snpedia_db_path
                                    .as_ref()
                                    .is_some_and(|p| p.exists());

                                db_status_row(
                                    ui,
                                    db_ok,
                                    if db_ok { "Main database ready" } else { "Main database not found" },
                                    if db_ok {
                                        Some(data.main_db_path.display().to_string())
                                    } else {
                                        None
                                    },
                                );
                                ui.add_space(4.0);
                                db_status_row(
                                    ui,
                                    snp_ok,
                                    if snp_ok { "SNPedia database ready" } else { "SNPedia not available" },
                                    None,
                                );

                                if !db_ok {
                                    ui.add_space(10.0);
                                    egui::Frame::NONE
                                        .fill(theme::WARNING_BG)
                                        .corner_radius(egui::CornerRadius::same(6))
                                        .inner_margin(egui::Margin::same(10))
                                        .show(ui, |ui| {
                                            ui.label(
                                                RichText::new(
                                                    "Place genesight.db in ~/.genesight/ or next to the executable.",
                                                )
                                                .size(12.0)
                                                .color(theme::WARNING),
                                            );
                                        });
                                }
                            });
                    },
                );

                ui.add_space(24.0);

                // Privacy notice
                ui.horizontal(|ui| {
                    ui.add_space((available_width - card_width) / 2.0);
                    ui.label(
                        RichText::new("\u{1F512}")
                            .size(12.0),
                    );
                    ui.label(
                        RichText::new(
                            "All processing happens locally. No data ever leaves your machine.",
                        )
                        .size(12.0)
                        .color(theme::TEXT_MUTED),
                    );
                });

                ui.add_space(24.0);
            });
        });

    selected_file
}

fn format_pill(ui: &mut egui::Ui, label: &str) {
    egui::Frame::NONE
        .fill(theme::ACCENT_LIGHT)
        .corner_radius(egui::CornerRadius::same(12))
        .inner_margin(egui::Margin::symmetric(10, 4))
        .show(ui, |ui| {
            ui.label(
                RichText::new(label)
                    .size(11.0)
                    .color(theme::ACCENT),
            );
        });
}

fn feature_chip(ui: &mut egui::Ui, icon: &str, title: &str, subtitle: &str, width: f32) {
    ui.allocate_ui(Vec2::new(width, 52.0), |ui| {
        theme::section_frame()
            .inner_margin(egui::Margin::same(10))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(icon).size(16.0));
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 1.0;
                        ui.label(RichText::new(title).size(12.0).strong().color(theme::TEXT_PRIMARY));
                        ui.label(RichText::new(subtitle).size(10.0).color(theme::TEXT_MUTED));
                    });
                });
            });
    });
}

fn db_status_row(ui: &mut egui::Ui, ok: bool, label: &str, detail: Option<String>) {
    ui.horizontal(|ui| {
        let (icon, bg, fg) = if ok {
            ("\u{2713}", theme::SUCCESS_BG, theme::SUCCESS)
        } else {
            ("\u{2717}", theme::DANGER_BG, theme::DANGER)
        };
        // Status dot
        egui::Frame::NONE
            .fill(bg)
            .corner_radius(egui::CornerRadius::same(10))
            .inner_margin(egui::Margin::symmetric(4, 2))
            .show(ui, |ui| {
                ui.label(RichText::new(icon).size(11.0).strong().color(fg));
            });
        ui.label(RichText::new(label).size(12.0).color(theme::TEXT_PRIMARY));
        if let Some(path) = detail {
            ui.label(RichText::new(path).size(10.0).color(theme::TEXT_MUTED));
        }
    });
}
