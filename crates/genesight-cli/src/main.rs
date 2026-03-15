mod tui;

use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::thread;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use genesight_core::models::ConfidenceTier;
use genesight_core::report::OutputFormat;
use rusqlite::Connection;

use tui::AnalysisProgress;

#[derive(Parser)]
#[command(name = "genesight")]
#[command(about = "Open-source, privacy-first DNA analysis tool")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(clap::Subcommand)]
enum Commands {
    /// Analyze a DNA file and generate a report
    Analyze {
        /// Path to the DNA raw data file (23andMe, AncestryDNA, or VCF)
        file: PathBuf,

        /// Output format
        #[arg(long, short, default_value = "markdown")]
        format: Format,

        /// Which tiers to include (comma-separated, e.g., "1,2" or "1,2,3")
        #[arg(long, short, default_value = "1,2")]
        tiers: String,

        /// Path to the main GeneSight database
        #[arg(long, default_value = "~/.genesight/genesight.db")]
        db: PathBuf,

        /// Path to the optional SNPedia database (CC-BY-NC-SA 3.0)
        #[arg(long)]
        snpedia_db: Option<PathBuf>,

        /// Output file (defaults to stdout)
        #[arg(long, short)]
        output: Option<PathBuf>,

        /// Show all annotated variants, not just notable ones
        #[arg(long, short)]
        verbose: bool,

        /// Launch interactive TUI dashboard instead of printing to stdout
        #[arg(long, short)]
        interactive: bool,
    },

    /// Download and update reference databases
    Fetch {
        /// Download all databases
        #[arg(long)]
        all: bool,

        /// Download ClinVar
        #[arg(long)]
        clinvar: bool,

        /// Download GWAS Catalog
        #[arg(long)]
        gwas: bool,

        /// Download SNPedia (CC-BY-NC-SA 3.0)
        #[arg(long)]
        snpedia: bool,

        /// Download gnomAD allele frequencies
        #[arg(long)]
        gnomad: bool,

        /// Download PharmGKB
        #[arg(long)]
        pharmgkb: bool,

        /// Database directory
        #[arg(long, default_value = "~/.genesight/")]
        db_dir: PathBuf,
    },

    /// Show database status and statistics
    Info {
        /// Path to the main database
        #[arg(long, default_value = "~/.genesight/genesight.db")]
        db: PathBuf,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Format {
    Markdown,
    Json,
    Html,
}

impl From<Format> for OutputFormat {
    fn from(f: Format) -> Self {
        match f {
            Format::Markdown => OutputFormat::Markdown,
            Format::Json => OutputFormat::Json,
            Format::Html => OutputFormat::Html,
        }
    }
}

/// Parse the tiers string (e.g., "1,2,3") into a Vec of ConfidenceTier.
fn parse_tiers(tiers: &str) -> Result<Vec<ConfidenceTier>> {
    tiers
        .split(',')
        .map(|s| match s.trim() {
            "1" => Ok(ConfidenceTier::Tier1Reliable),
            "2" => Ok(ConfidenceTier::Tier2Probable),
            "3" => Ok(ConfidenceTier::Tier3Speculative),
            other => anyhow::bail!("invalid tier '{other}', expected 1, 2, or 3"),
        })
        .collect()
}

/// Expand `~` to the user's home directory.
fn expand_tilde(path: &Path) -> PathBuf {
    if let Some(s) = path.to_str() {
        if let Some(rest) = s.strip_prefix("~/") {
            if let Some(home) = dirs::home_dir() {
                return home.join(rest);
            }
        }
    }
    path.to_path_buf()
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Analyze {
            file,
            format,
            tiers,
            db,
            snpedia_db,
            output,
            verbose,
            interactive,
        } => {
            let db_path = expand_tilde(&db);
            let _ = verbose;

            let tier_filter = parse_tiers(&tiers)?;

            if interactive {
                // Launch TUI immediately with a loading screen.
                // Analysis runs in a background thread, sending progress via mpsc.
                let (tx, rx) = mpsc::channel::<AnalysisProgress>();

                let snpedia_path = snpedia_db.map(|p| expand_tilde(&p));

                thread::spawn(move || {
                    run_analysis_background(tx, &file, &db_path, snpedia_path, &tier_filter);
                });

                tui::run_with_analysis(rx)?;
            } else {
                // Non-interactive path: run analysis synchronously and print output.
                eprintln!("Reading {}...", file.display());
                let content = std::fs::read_to_string(&file)
                    .with_context(|| format!("cannot read {}", file.display()))?;

                let parsed = genesight_core::parser::parse_auto_with_metadata(&content)?;
                eprintln!(
                    "Parsed {} variants (assembly: {})",
                    parsed.variants.len(),
                    parsed.assembly
                );

                let main_conn = Connection::open(&db_path)
                    .with_context(|| format!("cannot open database {}", db_path.display()))?;

                let db_assembly = genesight_core::db::query_db_assembly(&main_conn);

                let snpedia_conn = snpedia_db
                    .map(|p| {
                        let p = expand_tilde(&p);
                        Connection::open(&p).with_context(|| {
                            format!("cannot open SNPedia database {}", p.display())
                        })
                    })
                    .transpose()?;

                eprintln!("Analyzing against database...");
                let report = genesight_core::analyze_with_assembly(
                    &parsed.variants,
                    &main_conn,
                    snpedia_conn.as_ref(),
                    &tier_filter,
                    parsed.assembly,
                    db_assembly,
                )?;

                eprintln!(
                    "Found {} annotated variants, {} scored results",
                    report.annotated_variants,
                    report.results.len()
                );

                let rendered = genesight_core::report::render(&report, format.into())?;

                match output {
                    Some(path) => {
                        std::fs::write(&path, &rendered)
                            .with_context(|| format!("cannot write {}", path.display()))?;
                        eprintln!("Report written to {}", path.display());
                    }
                    None => {
                        println!("{rendered}");
                    }
                }
            }
        }
        Commands::Fetch {
            all,
            clinvar,
            gwas,
            snpedia,
            gnomad,
            pharmgkb,
            db_dir,
        } => {
            let _ = (all, clinvar, gwas, snpedia, gnomad, pharmgkb);
            println!(
                "Database fetch not yet implemented. Target dir: {}",
                db_dir.display()
            );
        }
        Commands::Info { db } => {
            let db_path = expand_tilde(&db);
            let conn = Connection::open(&db_path)
                .with_context(|| format!("cannot open database {}", db_path.display()))?;

            println!("GeneSight Database: {}", db_path.display());
            println!("---");

            let tables = [
                ("variants", "Reference variants"),
                ("clinvar", "ClinVar clinical classifications"),
                ("gwas", "GWAS Catalog associations"),
                ("frequencies", "Allele frequencies"),
                ("pharmacogenomics", "Pharmacogenomic annotations"),
            ];

            for (table, desc) in &tables {
                match conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |r| {
                    r.get::<_, i64>(0)
                }) {
                    Ok(count) => println!("  {desc}: {count} entries"),
                    Err(_) => println!("  {desc}: (table not found)"),
                }
            }
        }
    }

    Ok(())
}

/// Run the full analysis pipeline in a background thread, sending progress
/// updates through the provided channel.
///
/// This function is designed to be called from `thread::spawn`. It catches all
/// errors and sends them as `AnalysisProgress::Error` rather than panicking.
fn run_analysis_background(
    tx: mpsc::Sender<AnalysisProgress>,
    file: &Path,
    db_path: &Path,
    snpedia_path: Option<PathBuf>,
    tier_filter: &[ConfidenceTier],
) {
    // Helper: send a message, returning false if the receiver is gone.
    let send = |msg: AnalysisProgress| -> bool { tx.send(msg).is_ok() };

    // Step 1: Read the DNA file
    if !send(AnalysisProgress::Stage("Reading DNA file...".to_string())) {
        return;
    }
    let content = match std::fs::read_to_string(file) {
        Ok(c) => c,
        Err(e) => {
            let _ = send(AnalysisProgress::Error(format!(
                "Cannot read {}: {e}",
                file.display()
            )));
            return;
        }
    };

    // Step 2: Parse variants
    if !send(AnalysisProgress::Stage("Parsing variants...".to_string())) {
        return;
    }
    let parsed = match genesight_core::parser::parse_auto_with_metadata(&content) {
        Ok(p) => p,
        Err(e) => {
            let _ = send(AnalysisProgress::Error(format!("Parse error: {e}")));
            return;
        }
    };
    let variants = parsed.variants;
    let input_assembly = parsed.assembly;
    if !send(AnalysisProgress::FileRead(variants.len())) {
        return;
    }

    // Step 3: Open databases
    if !send(AnalysisProgress::Stage("Opening database...".to_string())) {
        return;
    }
    let main_conn = match Connection::open(db_path) {
        Ok(c) => c,
        Err(e) => {
            let _ = send(AnalysisProgress::Error(format!(
                "Cannot open database {}: {e}",
                db_path.display()
            )));
            return;
        }
    };

    let snpedia_conn = match snpedia_path {
        Some(ref p) => match Connection::open(p) {
            Ok(c) => Some(c),
            Err(e) => {
                let _ = send(AnalysisProgress::Error(format!(
                    "Cannot open SNPedia database {}: {e}",
                    p.display()
                )));
                return;
            }
        },
        None => None,
    };

    // Step 4: Annotate
    if !send(AnalysisProgress::Annotating("ClinVar".to_string())) {
        return;
    }
    // The core library runs annotation in one call; we send progress markers
    // for the TUI but the actual work is a single batch operation.
    if !send(AnalysisProgress::Annotating("GWAS Catalog".to_string())) {
        return;
    }
    if snpedia_conn.is_some() && !send(AnalysisProgress::Annotating("SNPedia".to_string())) {
        return;
    }
    if !send(AnalysisProgress::Stage(
        "Annotating variants...".to_string(),
    )) {
        return;
    }

    let db_assembly = genesight_core::db::query_db_assembly(&main_conn);

    let report = match genesight_core::analyze_with_assembly(
        &variants,
        &main_conn,
        snpedia_conn.as_ref(),
        tier_filter,
        input_assembly,
        db_assembly,
    ) {
        Ok(r) => r,
        Err(e) => {
            let _ = send(AnalysisProgress::Error(format!("Analysis failed: {e}")));
            return;
        }
    };

    // Step 5: Scoring (already done inside analyze, but signal the TUI)
    if !send(AnalysisProgress::Stage("Scoring variants...".to_string())) {
        return;
    }

    // Step 6: Complete
    let _ = send(AnalysisProgress::Complete(report));
}
