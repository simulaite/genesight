use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, ValueEnum};

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

        /// Which tiers to include (comma-separated, e.g., "1,2,3")
        #[arg(long, short, default_value = "1,2,3")]
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

        /// Omit disclaimer (for piping/scripting)
        #[arg(long)]
        no_disclaimer: bool,

        /// Show all annotated variants, not just notable ones
        #[arg(long, short)]
        verbose: bool,
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
            no_disclaimer,
            verbose,
        } => {
            tracing::info!(?file, ?format, ?tiers, ?db, "Starting analysis");

            let content =
                std::fs::read_to_string(&file).map_err(|e| anyhow::anyhow!("{}: {e}", file.display()))?;

            let variants = genesight_core::parser::parse_auto(&content)?;
            tracing::info!(count = variants.len(), "Parsed variants");

            // TODO: Open database, annotate, score, generate report
            let _ = (db, snpedia_db, output, no_disclaimer, verbose, format, tiers);

            println!(
                "Parsed {} variants from {}",
                variants.len(),
                file.display()
            );
            println!("Database annotation not yet implemented — coming soon!");
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
            println!(
                "Database info not yet implemented. Target: {}",
                db.display()
            );
        }
    }

    Ok(())
}
