use std::time::Instant;

use genesight_core::models::{ConfidenceTier, Report, ScoredResult};

/// Progress messages sent from the analysis thread to the TUI.
pub enum AnalysisProgress {
    /// Description of the current processing stage.
    Stage(String),
    /// Number of variants parsed from the DNA file.
    FileRead(usize),
    /// Which database is currently being queried.
    Annotating(String),
    /// Analysis finished successfully.
    Complete(Report),
    /// Analysis failed with an error message.
    Error(String),
}

/// Which panel currently has keyboard focus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Panel {
    Results,
    Details,
}

/// Which view is active in the dashboard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewMode {
    /// High-level summary with key findings and category breakdown.
    Summary,
    /// Full scrollable table of all results.
    Table,
}

/// Two-phase state: loading while analysis runs, dashboard when done.
pub enum AppPhase {
    /// Shown while the background analysis thread is running.
    Loading {
        /// Ordered list of stage descriptions and whether each is completed.
        stages: Vec<(String, bool)>,
        /// Description of the stage currently in progress.
        current_stage: String,
        /// When the loading phase started (for elapsed time display).
        elapsed: Instant,
        /// Tick counter driving the spinner animation.
        spinner_tick: usize,
    },
    /// Interactive dashboard shown after analysis completes.
    Dashboard {
        /// Which view is currently active.
        view_mode: ViewMode,
        /// Scroll offset for the summary view.
        summary_scroll: u16,
        /// Indices into `report.results` matching the current filter.
        filtered_indices: Vec<usize>,
        /// Currently selected index within `filtered_indices`.
        selected_index: usize,
        /// Scroll offset for the details panel.
        detail_scroll: u16,
        /// Which panel has keyboard focus.
        active_panel: Panel,
        /// Optional tier filter. `None` means show all tiers.
        tier_filter: Option<ConfidenceTier>,
        /// Whether the search input bar is active.
        search_mode: bool,
        /// Current search query text.
        search_query: String,
    },
}

/// Application state for the TUI dashboard.
pub struct App {
    /// The full analysis report. `None` during loading, `Some` after completion.
    pub report: Option<Report>,
    /// Current phase of the application.
    pub phase: AppPhase,
    /// Whether the application should quit.
    pub should_quit: bool,
    /// Error message from the analysis thread, if any.
    pub error_message: Option<String>,
}

impl App {
    /// Create a new App in loading state (used with `run_with_analysis`).
    pub fn new_loading() -> Self {
        Self {
            report: None,
            phase: AppPhase::Loading {
                stages: Vec::new(),
                current_stage: "Initializing...".to_string(),
                elapsed: Instant::now(),
                spinner_tick: 0,
            },
            should_quit: false,
            error_message: None,
        }
    }

    /// Create a new App directly from an analysis report (legacy path).
    #[allow(dead_code)]
    pub fn new(report: Report) -> Self {
        let filtered_indices: Vec<usize> = (0..report.results.len()).collect();
        Self {
            report: Some(report),
            phase: AppPhase::Dashboard {
                view_mode: ViewMode::Summary,
                summary_scroll: 0,
                filtered_indices,
                selected_index: 0,
                detail_scroll: 0,
                active_panel: Panel::Results,
                tier_filter: None,
                search_mode: false,
                search_query: String::new(),
            },
            should_quit: false,
            error_message: None,
        }
    }

    /// Handle an incoming progress message from the analysis thread.
    pub fn handle_progress(&mut self, progress: AnalysisProgress) {
        match progress {
            AnalysisProgress::Stage(desc) => {
                if let AppPhase::Loading {
                    stages,
                    current_stage,
                    ..
                } = &mut self.phase
                {
                    // Mark the previous current stage as completed
                    if !current_stage.is_empty() && current_stage != "Initializing..." {
                        // Check if this stage is already tracked
                        if let Some(entry) = stages.iter_mut().find(|(s, _)| s == current_stage) {
                            entry.1 = true;
                        } else {
                            stages.push((current_stage.clone(), true));
                        }
                    }
                    *current_stage = desc;
                }
            }
            AnalysisProgress::FileRead(count) => {
                if let AppPhase::Loading {
                    stages,
                    current_stage,
                    ..
                } = &mut self.phase
                {
                    // Mark previous stage as done
                    if !current_stage.is_empty() && current_stage != "Initializing..." {
                        if let Some(entry) = stages.iter_mut().find(|(s, _)| s == current_stage) {
                            entry.1 = true;
                        } else {
                            stages.push((current_stage.clone(), true));
                        }
                    }
                    let msg = format!("Parsed {count} variants");
                    stages.push((msg.clone(), true));
                    *current_stage = "Opening database...".to_string();
                }
            }
            AnalysisProgress::Annotating(db_name) => {
                if let AppPhase::Loading {
                    stages,
                    current_stage,
                    ..
                } = &mut self.phase
                {
                    // Mark previous stage as done
                    if !current_stage.is_empty() {
                        if let Some(entry) = stages.iter_mut().find(|(s, _)| s == current_stage) {
                            entry.1 = true;
                        } else {
                            stages.push((current_stage.clone(), true));
                        }
                    }
                    let msg = format!("Querying {db_name}...");
                    *current_stage = msg;
                }
            }
            AnalysisProgress::Complete(report) => {
                // Transition from Loading to Dashboard
                if let AppPhase::Loading {
                    stages,
                    current_stage,
                    ..
                } = &mut self.phase
                {
                    // Mark the last stage as complete
                    if !current_stage.is_empty() {
                        if let Some(entry) = stages.iter_mut().find(|(s, _)| s == current_stage) {
                            entry.1 = true;
                        } else {
                            stages.push((current_stage.clone(), true));
                        }
                    }
                    let summary = format!(
                        "Analysis complete: {} annotated, {} scored results",
                        report.annotated_variants,
                        report.results.len()
                    );
                    stages.push((summary, true));
                }

                let filtered_indices: Vec<usize> = (0..report.results.len()).collect();
                self.phase = AppPhase::Dashboard {
                    view_mode: ViewMode::Summary,
                    summary_scroll: 0,
                    filtered_indices,
                    selected_index: 0,
                    detail_scroll: 0,
                    active_panel: Panel::Results,
                    tier_filter: None,
                    search_mode: false,
                    search_query: String::new(),
                };
                self.report = Some(report);
            }
            AnalysisProgress::Error(msg) => {
                self.error_message = Some(msg);
            }
        }
    }

    /// Increment the spinner tick counter (only in loading phase).
    pub fn tick_spinner(&mut self) {
        if let AppPhase::Loading { spinner_tick, .. } = &mut self.phase {
            *spinner_tick = spinner_tick.wrapping_add(1);
        }
    }

    /// Toggle between Summary and Table views.
    pub fn toggle_view(&mut self) {
        if let AppPhase::Dashboard { view_mode, .. } = &mut self.phase {
            *view_mode = match *view_mode {
                ViewMode::Summary => ViewMode::Table,
                ViewMode::Table => ViewMode::Summary,
            };
        }
    }

    /// Get the current view mode.
    pub fn view_mode(&self) -> ViewMode {
        if let AppPhase::Dashboard { view_mode, .. } = &self.phase {
            *view_mode
        } else {
            ViewMode::Summary
        }
    }

    /// Scroll the summary view down.
    pub fn scroll_summary_down(&mut self) {
        if let AppPhase::Dashboard { summary_scroll, .. } = &mut self.phase {
            *summary_scroll = summary_scroll.saturating_add(1);
        }
    }

    /// Scroll the summary view up.
    pub fn scroll_summary_up(&mut self) {
        if let AppPhase::Dashboard { summary_scroll, .. } = &mut self.phase {
            *summary_scroll = summary_scroll.saturating_sub(1);
        }
    }

    /// Get summary scroll position.
    pub fn summary_scroll(&self) -> u16 {
        if let AppPhase::Dashboard { summary_scroll, .. } = &self.phase {
            *summary_scroll
        } else {
            0
        }
    }

    /// Collect category counts from the report.
    pub fn category_counts(&self) -> Vec<(genesight_core::models::ResultCategory, usize)> {
        use genesight_core::models::ResultCategory;
        let report = match &self.report {
            Some(r) => r,
            None => return Vec::new(),
        };
        let categories = [
            ResultCategory::MonogenicDisease,
            ResultCategory::CarrierStatus,
            ResultCategory::Pharmacogenomics,
            ResultCategory::GwasAssociation,
            ResultCategory::PhysicalTrait,
            ResultCategory::ComplexTrait,
            ResultCategory::ClinVarConflicting,
        ];
        categories
            .iter()
            .filter_map(|&cat| {
                let count = report.results.iter().filter(|r| r.category == cat).count();
                if count > 0 {
                    Some((cat, count))
                } else {
                    None
                }
            })
            .collect()
    }

    // --- Dashboard navigation methods ----------------------------------------
    // These are no-ops if not in Dashboard phase or if no report is loaded.

    /// Move selection down by one.
    pub fn next(&mut self) {
        if let AppPhase::Dashboard {
            filtered_indices,
            selected_index,
            detail_scroll,
            ..
        } = &mut self.phase
        {
            if !filtered_indices.is_empty() {
                *selected_index = (*selected_index + 1).min(filtered_indices.len() - 1);
                *detail_scroll = 0;
            }
        }
    }

    /// Move selection up by one.
    pub fn previous(&mut self) {
        if let AppPhase::Dashboard {
            selected_index,
            detail_scroll,
            ..
        } = &mut self.phase
        {
            *selected_index = selected_index.saturating_sub(1);
            *detail_scroll = 0;
        }
    }

    /// Jump to the first item.
    pub fn go_to_top(&mut self) {
        if let AppPhase::Dashboard {
            selected_index,
            detail_scroll,
            ..
        } = &mut self.phase
        {
            *selected_index = 0;
            *detail_scroll = 0;
        }
    }

    /// Jump to the last item.
    pub fn go_to_bottom(&mut self) {
        if let AppPhase::Dashboard {
            filtered_indices,
            selected_index,
            detail_scroll,
            ..
        } = &mut self.phase
        {
            if !filtered_indices.is_empty() {
                *selected_index = filtered_indices.len() - 1;
            }
            *detail_scroll = 0;
        }
    }

    /// Scroll the detail panel down.
    pub fn scroll_detail_down(&mut self) {
        if let AppPhase::Dashboard { detail_scroll, .. } = &mut self.phase {
            *detail_scroll = detail_scroll.saturating_add(1);
        }
    }

    /// Scroll the detail panel up.
    pub fn scroll_detail_up(&mut self) {
        if let AppPhase::Dashboard { detail_scroll, .. } = &mut self.phase {
            *detail_scroll = detail_scroll.saturating_sub(1);
        }
    }

    /// Toggle the tier filter. If the same tier is already active, clear it.
    pub fn toggle_tier(&mut self, tier: ConfidenceTier) {
        if let AppPhase::Dashboard { tier_filter, .. } = &mut self.phase {
            if *tier_filter == Some(tier) {
                *tier_filter = None;
            } else {
                *tier_filter = Some(tier);
            }
        }
        self.apply_filter();
    }

    /// Show all tiers (clear filter).
    pub fn show_all(&mut self) {
        if let AppPhase::Dashboard { tier_filter, .. } = &mut self.phase {
            *tier_filter = None;
        }
        self.apply_filter();
    }

    /// Recompute `filtered_indices` based on current tier filter and search query.
    pub fn apply_filter(&mut self) {
        let report = match &self.report {
            Some(r) => r,
            None => return,
        };

        if let AppPhase::Dashboard {
            filtered_indices,
            selected_index,
            detail_scroll,
            tier_filter,
            search_query,
            ..
        } = &mut self.phase
        {
            let query_lower = search_query.to_lowercase();

            *filtered_indices = report
                .results
                .iter()
                .enumerate()
                .filter(|(_, result)| {
                    // Tier filter
                    if let Some(tier) = tier_filter {
                        if result.tier != *tier {
                            return false;
                        }
                    }
                    // Search filter
                    if !query_lower.is_empty() {
                        let matches_rsid = result
                            .variant
                            .variant
                            .rsid
                            .as_deref()
                            .is_some_and(|rsid| rsid.to_lowercase().contains(&query_lower));

                        let matches_gene = gene_name_for(report, result)
                            .to_lowercase()
                            .contains(&query_lower);

                        let matches_summary = result.summary.to_lowercase().contains(&query_lower);

                        if !matches_rsid && !matches_gene && !matches_summary {
                            return false;
                        }
                    }
                    true
                })
                .map(|(i, _)| i)
                .collect();

            // Clamp selected index
            if filtered_indices.is_empty() {
                *selected_index = 0;
            } else {
                *selected_index = (*selected_index).min(filtered_indices.len() - 1);
            }
            *detail_scroll = 0;
        }
    }

    /// Get the currently selected scored result, if any.
    pub fn selected_result(&self) -> Option<&ScoredResult> {
        let report = self.report.as_ref()?;
        if let AppPhase::Dashboard {
            filtered_indices,
            selected_index,
            ..
        } = &self.phase
        {
            filtered_indices
                .get(*selected_index)
                .and_then(|&idx| report.results.get(idx))
        } else {
            None
        }
    }

    /// Extract the best gene name from a scored result.
    pub fn gene_name(&self, result: &ScoredResult) -> String {
        self.report
            .as_ref()
            .map(|r| gene_name_for(r, result))
            .unwrap_or_else(|| "\u{2014}".to_string())
    }

    /// Count results by tier (across all results, not filtered).
    pub fn tier_counts(&self) -> (usize, usize, usize) {
        let report = match &self.report {
            Some(r) => r,
            None => return (0, 0, 0),
        };
        let mut t1 = 0;
        let mut t2 = 0;
        let mut t3 = 0;
        for r in &report.results {
            match r.tier {
                ConfidenceTier::Tier1Reliable => t1 += 1,
                ConfidenceTier::Tier2Probable => t2 += 1,
                ConfidenceTier::Tier3Speculative => t3 += 1,
            }
        }
        (t1, t2, t3)
    }

    /// Access dashboard-specific fields. Returns None if in loading phase.
    pub fn dashboard_state(&self) -> Option<DashboardView<'_>> {
        if let AppPhase::Dashboard {
            view_mode,
            summary_scroll: _,
            filtered_indices,
            selected_index,
            detail_scroll,
            active_panel,
            tier_filter,
            search_mode,
            search_query,
        } = &self.phase
        {
            Some(DashboardView {
                view_mode: *view_mode,
                filtered_indices,
                selected_index: *selected_index,
                detail_scroll: *detail_scroll,
                active_panel: *active_panel,
                tier_filter: *tier_filter,
                search_mode: *search_mode,
                search_query,
            })
        } else {
            None
        }
    }
}

/// Borrowed view into dashboard-phase state, used by the UI layer.
pub struct DashboardView<'a> {
    pub view_mode: ViewMode,
    pub filtered_indices: &'a [usize],
    pub selected_index: usize,
    pub detail_scroll: u16,
    pub active_panel: Panel,
    pub tier_filter: Option<ConfidenceTier>,
    pub search_mode: bool,
    pub search_query: &'a str,
}

/// Extract the best gene name from a scored result (free function for use in filtering).
fn gene_name_for(_report: &Report, result: &ScoredResult) -> String {
    // Try ClinVar gene_symbol first
    if let Some(ref cv) = result.variant.clinvar {
        if let Some(ref gene) = cv.gene_symbol {
            return gene.clone();
        }
    }
    // Try pharmacogenomics gene
    if let Some(ref pharma) = result.variant.pharmacogenomics {
        return pharma.gene.clone();
    }
    // Try GWAS mapped gene
    for hit in &result.variant.gwas_hits {
        if let Some(ref gene) = hit.mapped_gene {
            return gene.clone();
        }
    }
    "\u{2014}".to_string()
}
