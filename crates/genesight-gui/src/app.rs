use std::path::PathBuf;
use std::sync::mpsc;

use eframe::egui;
use genesight_core::models::{AnnotationConfig, ConfidenceTier};

use crate::export;
use crate::state::{AnalysisProgress, AppData, AppState};
use crate::theme;
use crate::views;

/// Main application struct.
pub struct GeneSightApp {
    pub state: AppState,
    pub data: AppData,
}

impl GeneSightApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        theme::apply_theme(&cc.egui_ctx);
        Self {
            state: AppState::Welcome,
            data: AppData::new(),
        }
    }

    /// Start the analysis pipeline in a background thread.
    fn start_analysis(&mut self, file_path: PathBuf) {
        let file_name = file_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".into());

        self.data.dna_file = Some(file_path.clone());

        let (tx, rx) = mpsc::channel();
        let main_db_path = self.data.main_db_path.clone();
        let snpedia_db_path = if self.data.snpedia_enabled {
            self.data.snpedia_db_path.clone()
        } else {
            None
        };
        let config = self.data.annotation_config.clone();

        self.state = AppState::Analyzing {
            file_name: file_name.clone(),
            stages: Vec::new(),
            current_stage: "Reading file...".into(),
            rx,
        };

        std::thread::spawn(move || {
            run_analysis(file_path, main_db_path, snpedia_db_path, config, tx);
        });
    }
}

impl eframe::App for GeneSightApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Menu bar
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open DNA File...").clicked() {
                        ui.close_menu();
                        if let Some(path) = pick_dna_file() {
                            self.start_analysis(path);
                        }
                    }
                    ui.separator();
                    if let Some(ref report) = self.data.report {
                        if ui.button("Export as HTML...").clicked() {
                            ui.close_menu();
                            if let Some(path) = rfd::FileDialog::new()
                                .set_title("Save HTML Report")
                                .add_filter("HTML", &["html"])
                                .set_file_name("genesight_report.html")
                                .save_file()
                            {
                                match export::export_html(report, &path) {
                                    Ok(()) => {
                                        self.data.error_message = None;
                                        tracing::info!("HTML exported to {}", path.display());
                                    }
                                    Err(e) => {
                                        self.data.error_message =
                                            Some(format!("HTML export failed: {e}"));
                                    }
                                }
                            }
                        }
                        if ui.button("Export as PDF...").clicked() {
                            ui.close_menu();
                            if let Some(path) = rfd::FileDialog::new()
                                .set_title("Save PDF Report")
                                .add_filter("PDF", &["pdf"])
                                .set_file_name("genesight_report.pdf")
                                .save_file()
                            {
                                match export::export_pdf(report, &path) {
                                    Ok(()) => {
                                        self.data.error_message = None;
                                        tracing::info!("PDF exported to {}", path.display());
                                    }
                                    Err(e) => {
                                        self.data.error_message =
                                            Some(format!("PDF export failed: {e}"));
                                    }
                                }
                            }
                        }
                        ui.separator();
                    }
                    #[cfg(not(target_arch = "wasm32"))]
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });

                if self.data.report.is_some() {
                    ui.menu_button("View", |ui| {
                        if ui.button("Back to Welcome").clicked() {
                            ui.close_menu();
                            self.state = AppState::Welcome;
                            self.data.report = None;
                            self.data.dna_file = None;
                            self.data.filtered_indices.clear();
                            self.data.search_query.clear();
                            self.data.selected_index = 0;
                        }
                    });

                    ui.menu_button("Analyze", |ui| {
                        if ui.button("Re-analyze with current settings").clicked() {
                            ui.close_menu();
                            if let Some(path) = self.data.dna_file.clone() {
                                self.start_analysis(path);
                            }
                        }
                    });
                }
            });
        });

        // Error toast
        if let Some(ref msg) = self.data.error_message.clone() {
            egui::TopBottomPanel::bottom("error_bar").show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new(msg)
                            .color(theme::DANGER)
                            .size(12.0),
                    );
                    if ui.small_button("Dismiss").clicked() {
                        self.data.error_message = None;
                    }
                });
            });
        }

        // Main content based on state
        egui::CentralPanel::default().show(ctx, |ui| {
            match &mut self.state {
                AppState::Welcome => {
                    if let Some(path) = views::welcome::draw(ui, &self.data) {
                        // File was selected — defer start_analysis to avoid borrow conflict
                        let p = path;
                        // We need to break out and call start_analysis
                        // Use a temporary to avoid borrow issues
                        self.data.error_message = None;
                        self.start_analysis(p);
                        return;
                    }
                }
                AppState::Analyzing {
                    file_name,
                    stages,
                    current_stage,
                    rx,
                } => {
                    // Poll for progress
                    while let Ok(msg) = rx.try_recv() {
                        match msg {
                            AnalysisProgress::Stage(desc) => {
                                if !current_stage.is_empty() {
                                    stages.push((current_stage.clone(), true));
                                }
                                *current_stage = desc;
                            }
                            AnalysisProgress::ParsedVariants(count) => {
                                stages.push((
                                    format!("Parsed {count} variants"),
                                    true,
                                ));
                            }
                            AnalysisProgress::AnnotationComplete { annotated } => {
                                stages.push((
                                    format!("Annotated {annotated} variants"),
                                    true,
                                ));
                            }
                            AnalysisProgress::Complete(report) => {
                                self.data.report = Some(*report);
                                self.data.apply_filter();
                                self.state = AppState::Dashboard;
                                return;
                            }
                            AnalysisProgress::Error(err) => {
                                self.data.error_message = Some(err);
                                self.state = AppState::Welcome;
                                return;
                            }
                        }
                    }

                    views::progress::draw(ui, file_name, stages, current_stage);
                    ctx.request_repaint();
                }
                AppState::Dashboard => {
                    views::dashboard::draw(ui, &mut self.data);
                }
            }
        });
    }
}

fn pick_dna_file() -> Option<PathBuf> {
    rfd::FileDialog::new()
        .set_title("Open DNA Raw Data File")
        .add_filter("DNA Files", &["txt", "csv", "tsv", "vcf", "gz"])
        .add_filter("All Files", &["*"])
        .pick_file()
}

/// Run the analysis pipeline on a background thread, sending progress updates.
fn run_analysis(
    file_path: PathBuf,
    main_db_path: PathBuf,
    snpedia_db_path: Option<PathBuf>,
    config: AnnotationConfig,
    tx: mpsc::Sender<AnalysisProgress>,
) {
    let send = |msg: AnalysisProgress| {
        let _ = tx.send(msg);
    };

    // Read and parse the file
    send(AnalysisProgress::Stage("Reading DNA file...".into()));

    let file_data = match std::fs::read_to_string(&file_path) {
        Ok(d) => d,
        Err(e) => {
            send(AnalysisProgress::Error(format!(
                "Failed to read file: {e}"
            )));
            return;
        }
    };

    send(AnalysisProgress::Stage("Parsing variants...".into()));

    let variants = match genesight_core::parser::parse_auto(&file_data) {
        Ok(v) => v,
        Err(e) => {
            send(AnalysisProgress::Error(format!("Parse error: {e}")));
            return;
        }
    };

    send(AnalysisProgress::ParsedVariants(variants.len()));

    // Open databases
    send(AnalysisProgress::Stage("Opening databases...".into()));

    let main_db = match rusqlite::Connection::open(&main_db_path) {
        Ok(db) => db,
        Err(e) => {
            send(AnalysisProgress::Error(format!(
                "Failed to open main database at {}: {e}",
                main_db_path.display()
            )));
            return;
        }
    };

    let snpedia_conn = snpedia_db_path.and_then(|p| {
        if p.exists() {
            rusqlite::Connection::open(&p).ok()
        } else {
            None
        }
    });

    // Run analysis
    send(AnalysisProgress::Stage("Annotating variants...".into()));

    let tiers = vec![
        ConfidenceTier::Tier1Reliable,
        ConfidenceTier::Tier2Probable,
        ConfidenceTier::Tier3Speculative,
    ];

    match genesight_core::analyze_with_config(
        &variants,
        &main_db,
        snpedia_conn.as_ref(),
        &tiers,
        &config,
    ) {
        Ok(report) => {
            send(AnalysisProgress::AnnotationComplete {
                annotated: report.annotated_variants,
            });
            send(AnalysisProgress::Stage("Building report...".into()));
            send(AnalysisProgress::Complete(Box::new(report)));
        }
        Err(e) => {
            send(AnalysisProgress::Error(format!("Analysis failed: {e}")));
        }
    }
}
