use std::path::PathBuf;
use std::sync::mpsc;

use genesight_core::models::{AnnotationConfig, ConfidenceTier, Report, ScoredResult};

/// Progress messages sent from the analysis background thread.
#[derive(Debug)]
pub enum AnalysisProgress {
    Stage(String),
    ParsedVariants(usize),
    AnnotationComplete { annotated: usize },
    Complete(Box<Report>),
    Error(String),
}

/// Top-level application state.
pub enum AppState {
    /// No file loaded, show welcome/file picker.
    Welcome,
    /// Analysis is running in a background thread.
    Analyzing {
        file_name: String,
        stages: Vec<(String, bool)>,
        current_stage: String,
        rx: mpsc::Receiver<AnalysisProgress>,
    },
    /// Results are ready, show the dashboard.
    Dashboard,
}

/// Persistent app data that survives state transitions.
pub struct AppData {
    /// Path to the loaded DNA file.
    pub dna_file: Option<PathBuf>,
    /// The analysis report (set after analysis completes).
    pub report: Option<Report>,
    /// Main database path.
    pub main_db_path: PathBuf,
    /// SNPedia database path (optional).
    pub snpedia_db_path: Option<PathBuf>,
    /// Which datasets are enabled.
    pub annotation_config: AnnotationConfig,
    /// Whether SNPedia is enabled (separate DB).
    pub snpedia_enabled: bool,
    /// Tier filter checkboxes.
    pub tier_filter: [bool; 3],
    /// Search query in results table.
    pub search_query: String,
    /// Index of selected result in filtered list.
    pub selected_index: usize,
    /// Filtered result indices.
    pub filtered_indices: Vec<usize>,
    /// Error message to display.
    pub error_message: Option<String>,
    /// Show settings panel.
    pub show_settings: bool,
}

impl AppData {
    pub fn new() -> Self {
        // Look for databases in multiple locations
        let main_db_path = find_database("genesight.db");
        let snpedia_db_path = {
            let p = find_database("snpedia.db");
            if p.exists() { Some(p) } else { None }
        };

        Self {
            dna_file: None,
            report: None,
            main_db_path,
            snpedia_db_path: snpedia_db_path.clone(),
            annotation_config: AnnotationConfig::default(),
            snpedia_enabled: snpedia_db_path.is_some(),
            tier_filter: [true, true, true],
            search_query: String::new(),
            selected_index: 0,
            filtered_indices: Vec::new(),
            error_message: None,
            show_settings: false,
        }
    }

    /// Recompute filtered indices based on current tier filter and search.
    pub fn apply_filter(&mut self) {
        let report = match &self.report {
            Some(r) => r,
            None => {
                self.filtered_indices.clear();
                return;
            }
        };

        let query = self.search_query.to_lowercase();

        self.filtered_indices = report
            .results
            .iter()
            .enumerate()
            .filter(|(_, r)| {
                // Tier filter
                let tier_ok = match r.tier {
                    ConfidenceTier::Tier1Reliable => self.tier_filter[0],
                    ConfidenceTier::Tier2Probable => self.tier_filter[1],
                    ConfidenceTier::Tier3Speculative => self.tier_filter[2],
                };
                if !tier_ok {
                    return false;
                }

                // Search filter
                if !query.is_empty() {
                    let rsid_match = r
                        .variant
                        .variant
                        .rsid
                        .as_deref()
                        .is_some_and(|s| s.to_lowercase().contains(&query));
                    let gene_match = gene_name(r).to_lowercase().contains(&query);
                    let summary_match = r.summary.to_lowercase().contains(&query);
                    if !rsid_match && !gene_match && !summary_match {
                        return false;
                    }
                }

                true
            })
            .map(|(i, _)| i)
            .collect();

        // Clamp selection
        if self.filtered_indices.is_empty() {
            self.selected_index = 0;
        } else {
            self.selected_index = self.selected_index.min(self.filtered_indices.len() - 1);
        }
    }

    /// Get the currently selected result.
    pub fn selected_result(&self) -> Option<&ScoredResult> {
        let report = self.report.as_ref()?;
        self.filtered_indices
            .get(self.selected_index)
            .and_then(|&i| report.results.get(i))
    }

    /// Tier counts from the report.
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

    /// Check if the main database exists and is accessible.
    pub fn main_db_ok(&self) -> bool {
        self.main_db_path.exists()
    }
}

/// Extract the best gene name from a scored result.
pub fn gene_name(result: &ScoredResult) -> String {
    if let Some(ref cv) = result.variant.clinvar {
        if let Some(ref gene) = cv.gene_symbol {
            return gene.clone();
        }
    }
    if let Some(ref pharma) = result.variant.pharmacogenomics {
        return pharma.gene.clone();
    }
    for hit in &result.variant.gwas_hits {
        if let Some(ref gene) = hit.mapped_gene {
            return gene.clone();
        }
    }
    "\u{2014}".to_string()
}

/// Short category label.
pub fn short_category(cat: genesight_core::models::ResultCategory) -> &'static str {
    use genesight_core::models::ResultCategory;
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

/// Find a database file in standard search paths.
fn find_database(name: &str) -> PathBuf {
    // 1. Next to the executable
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let candidate = dir.join(name);
            if candidate.exists() {
                return candidate;
            }
            // Also check data/ subdirectory next to exe
            let candidate = dir.join("data").join(name);
            if candidate.exists() {
                return candidate;
            }
        }
    }

    // 2. data/seed/ in working directory (development)
    let dev = PathBuf::from("data/seed").join(name);
    if dev.exists() {
        return dev;
    }

    // 3. ~/.genesight/
    if let Some(home) = dirs::home_dir() {
        let candidate = home.join(".genesight").join(name);
        if candidate.exists() {
            return candidate;
        }
    }

    // Default path (may not exist yet)
    dirs::home_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(".genesight")
        .join(name)
}
