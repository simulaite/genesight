use egui::{Align, Layout, RichText, Vec2};

use crate::theme;

/// Draw the analysis progress screen.
pub fn draw(ui: &mut egui::Ui, file_name: &str, stages: &[(String, bool)], current_stage: &str) {
    let available_width = ui.available_width();

    ui.vertical_centered(|ui| {
        let top_space = (ui.available_height() * 0.15).clamp(40.0, 100.0);
        ui.add_space(top_space);

        // Spinner icon
        ui.label(RichText::new("\u{1F9EC}").size(40.0));
        ui.add_space(12.0);

        ui.label(
            RichText::new("Analyzing DNA Data")
                .size(24.0)
                .strong()
                .color(theme::TEXT_PRIMARY),
        );
        ui.add_space(4.0);
        ui.label(
            RichText::new(file_name)
                .size(13.0)
                .color(theme::TEXT_SECONDARY),
        );

        ui.add_space(32.0);

        // Progress card
        let card_width = available_width.clamp(320.0, 500.0) - 40.0;
        ui.allocate_ui_with_layout(
            Vec2::new(card_width, 0.0),
            Layout::top_down(Align::Center),
            |ui| {
                theme::card_frame()
                    .inner_margin(egui::Margin::same(28))
                    .show(ui, |ui| {
                        ui.set_min_width(card_width - 56.0);

                        let total_steps =
                            stages.len() + if current_stage.is_empty() { 0 } else { 1 };
                        let completed = stages.len();

                        // Completed stages
                        for (i, (desc, _done)) in stages.iter().enumerate() {
                            draw_step(ui, i + 1, desc, StepState::Done);
                            if i < total_steps - 1 {
                                draw_connector(ui, true);
                            }
                        }

                        // Current stage with spinner
                        if !current_stage.is_empty() {
                            draw_step(ui, completed + 1, current_stage, StepState::Active);
                        }

                        ui.add_space(20.0);

                        // Progress bar
                        let progress_frac = if total_steps > 0 {
                            completed as f32 / (total_steps as f32 + 1.0)
                        } else {
                            0.0
                        };

                        egui::Frame::NONE
                            .fill(theme::BG_SIDEBAR)
                            .corner_radius(egui::CornerRadius::same(4))
                            .show(ui, |ui| {
                                ui.set_min_width(card_width - 56.0);
                                ui.add(
                                    egui::ProgressBar::new(progress_frac)
                                        .animate(true)
                                        .desired_width(card_width - 56.0)
                                        .corner_radius(egui::CornerRadius::same(4)),
                                );
                            });
                    });
            },
        );
    });
}

enum StepState {
    Done,
    Active,
}

fn draw_step(ui: &mut egui::Ui, _step_num: usize, label: &str, state: StepState) {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 10.0;
        match state {
            StepState::Done => {
                // Green circle with check
                egui::Frame::NONE
                    .fill(theme::SUCCESS_BG)
                    .corner_radius(egui::CornerRadius::same(12))
                    .inner_margin(egui::Margin::symmetric(6, 3))
                    .show(ui, |ui| {
                        ui.label(
                            RichText::new("\u{2713}")
                                .size(12.0)
                                .strong()
                                .color(theme::SUCCESS),
                        );
                    });
                ui.label(RichText::new(label).size(13.0).color(theme::TEXT_PRIMARY));
            }
            StepState::Active => {
                // Accent circle with number
                egui::Frame::NONE
                    .fill(theme::ACCENT_LIGHT)
                    .corner_radius(egui::CornerRadius::same(12))
                    .inner_margin(egui::Margin::symmetric(6, 3))
                    .show(ui, |ui| {
                        ui.spinner();
                    });
                ui.label(
                    RichText::new(label)
                        .size(13.0)
                        .strong()
                        .color(theme::ACCENT),
                );
            }
        }
    });
}

fn draw_connector(ui: &mut egui::Ui, _completed: bool) {
    ui.horizontal(|ui| {
        ui.add_space(14.0);
        ui.label(RichText::new("\u{2502}").size(10.0).color(theme::BORDER));
    });
}
