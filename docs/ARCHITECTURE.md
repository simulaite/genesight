# GeneSight — Architecture

## Design Principles

1. **Privacy first** — DNA data never leaves the user's machine
2. **Offline-capable** — No network connection required after initial database download
3. **Modular** — Core library is independent of CLI/Desktop/Web
4. **Transparent** — Every result shows its source and confidence level

---

## Crate Structure

```
genesight (workspace)
│
├── genesight-core          # Library — no IO, no network, pure logic
│   ├── parser::*           # DNA files → Vec<Variant>
│   ├── db::*               # SQLite queries (takes Connection as parameter)
│   ├── annotator::*        # Variant + DB → AnnotatedVariant
│   ├── scorer::*           # AnnotatedVariant → ScoredResult + ConfidenceTier
│   ├── report::*           # Vec<ScoredResult> → Report (MD/JSON/HTML)
│   └── models::*           # Shared types
│
├── genesight-cli           # Binary — IO, arg parsing, DB opening
│   └── main.rs             # clap, rusqlite::Connection, calls into core
│
├── genesight-server        # Binary — Axum API (Phase 3)
│   └── main.rs             # Upload endpoint, calls into core
│
└── genesight-gui           # Binary — egui desktop app (Phase 2)
    └── ...                 # Thin wrapper around core
```

### Why This Separation?

`genesight-core` has **no** dependencies on:
- Filesystem IO (receives `&[u8]` or `&str`, not paths)
- Network (receives `rusqlite::Connection`, does not open it itself)
- CLI framework (no clap)
- Async runtime (everything is synchronous)

This makes the library:
- **Testable** — Unit tests do not need a filesystem
- **Reusable** — CLI, egui, Axum, WASM can all use the same library
- **Compilable for WASM** — Future option: run in the browser

---

## Data Flow

```
┌─────────────┐     ┌──────────┐     ┌───────────┐     ┌──────────┐     ┌──────────┐
│  DNA File   │ ──► │  Parser  │ ──► │ Annotator │ ──► │  Scorer  │ ──► │  Report  │
│ (23andMe/   │     │          │     │           │     │          │     │ (MD/JSON │
│  VCF/etc.)  │     │ → Vec<   │     │ → Vec<    │     │ → Vec<   │     │  /HTML)  │
│             │     │ Variant> │     │ Annotated │     │ Scored   │     │          │
│             │     │          │     │ Variant>  │     │ Result>  │     │          │
└─────────────┘     └──────────┘     └─────┬─────┘     └──────────┘     └──────────┘
                                           │
                                    ┌──────┴──────┐
                                    │  SQLite DB  │
                                    │ (genesight  │
                                    │    .db)     │
                                    └─────────────┘
```

### Step by Step

1. **Parser** reads the DNA file and produces `Vec<Variant>`
   - Automatically detects the format (23andMe vs AncestryDNA vs VCF)
   - Normalizes to a unified `Variant` struct
   - Filters invalid lines (--calls, no-calls)

2. **Annotator** takes each `Variant` and looks it up in the local DB
   - Batch queries for performance (not 600K individual queries)
   - Produces `AnnotatedVariant` with all discovered annotations

3. **Scorer** evaluates each annotated variant
   - Assigns a `ConfidenceTier` (Tier 1/2/3)
   - Calculates polygenic risk scores (summed across many variants)
   - Determines pharmacogenomic phenotypes (metabolizer status)

4. **Report** formats the results
   - Groups by tier and category
   - Adds disclaimer and attributions
   - Output: Markdown, JSON, or HTML

---

## Core Data Models

```rust
/// A single DNA variant from the user's file
pub struct Variant {
    pub rsid: Option<String>,       // e.g. "rs1234567"
    pub chromosome: String,          // "1"-"22", "X", "Y", "MT"
    pub position: u64,               // Genomic position
    pub genotype: Genotype,          // e.g. Genotype::Heterozygous('A', 'G')
    pub source_format: SourceFormat, // TwentyThreeAndMe, AncestryDNA, VCF
}

pub enum Genotype {
    Homozygous(char),           // AA, GG, etc.
    Heterozygous(char, char),   // AG, CT, etc.
    NoCall,                      // --
    Indel(String),              // Insertion/Deletion
}

pub enum SourceFormat {
    TwentyThreeAndMe,
    AncestryDNA,
    Vcf,
}

/// A variant with all database annotations
pub struct AnnotatedVariant {
    pub variant: Variant,
    pub clinvar: Option<ClinVarAnnotation>,
    pub snpedia: Option<SnpediaAnnotation>,
    pub gwas_hits: Vec<GwasHit>,
    pub frequency: Option<AlleleFrequency>,
    pub pharmacogenomics: Option<PharmaAnnotation>,
}

/// Scored result with confidence level
pub struct ScoredResult {
    pub variant: AnnotatedVariant,
    pub tier: ConfidenceTier,
    pub category: ResultCategory,
    pub summary: String,            // Human-readable summary
    pub details: String,            // More detailed explanation
}

pub enum ConfidenceTier {
    Tier1Reliable,    // >95% — clinically validated
    Tier2Probable,    // 60-85% — statistical association
    Tier3Speculative, // 50-65% — weak evidence
}

pub enum ResultCategory {
    MonogenicDisease,
    CarrierStatus,
    Pharmacogenomics,
    PolygenicRiskScore,
    PhysicalTrait,
    ComplexTrait,
    Ancestry,
}
```

---

## Performance Considerations

### DNA File Sizes
- 23andMe: ~600K variants, ~15MB text file
- AncestryDNA: ~700K variants, ~20MB
- WGS VCF: ~4-5M variants, ~1-5GB

### SQLite Query Strategy

**Not like this:**
```rust
// 600,000 individual queries
for variant in variants {
    db.query("SELECT * FROM clinvar WHERE rsid = ?", &[&variant.rsid]);
}
```

**But like this:**
```rust
// Batch queries with a temporary table
db.execute("CREATE TEMP TABLE user_variants (rsid TEXT PRIMARY KEY)");
// Bulk-insert the user variants
// Then JOIN against the annotation tables
db.query("
    SELECT uv.rsid, c.*
    FROM user_variants uv
    LEFT JOIN clinvar c ON uv.rsid = c.rsid
");
```

Expected performance:
- Parsing: <2 seconds for 600K variants
- Annotation (batch): <10 seconds against all databases
- Report generation: <1 second
- **Total: <15 seconds** for a complete report

---

## CLI Interface (Phase 1)

```
genesight 0.1.0
Open-source DNA analysis tool

USAGE:
    genesight <COMMAND>

COMMANDS:
    analyze     Analyze a DNA file and generate a report
    fetch       Download/update reference databases
    info        Show database status and statistics
    help        Print help information

ANALYZE:
    genesight analyze <FILE> [OPTIONS]

    OPTIONS:
        -f, --format <FORMAT>    Output format [default: markdown] [possible: markdown, json, html]
        -o, --output <FILE>      Output file [default: stdout]
        -t, --tiers <TIERS>      Which tiers to include [default: 1,2,3]
        --db <PATH>              Path to genesight.db [default: ~/.genesight/genesight.db]
        --snpedia-db <PATH>      Path to snpedia.db (optional)
        --no-disclaimer          Omit disclaimer (for piping/scripting)
        -v, --verbose            Show all annotated variants, not just notable ones

FETCH:
    genesight fetch [OPTIONS]

    OPTIONS:
        --all                    Download all databases
        --clinvar                Download ClinVar
        --gwas                   Download GWAS Catalog
        --snpedia                Download SNPedia (CC-BY-NC-SA 3.0)
        --gnomad                 Download gnomAD frequencies
        --pharmgkb               Download PharmGKB
        --db-dir <PATH>          Database directory [default: ~/.genesight/]

INFO:
    genesight info [OPTIONS]

    OPTIONS:
        --db <PATH>              Path to genesight.db
```

---

## Future: egui Desktop App (Phase 2)

```
┌─────────────────────────────────────────────────┐
│  egui Window (Native)                           │
│  ┌───────────────────────────────────────────┐  │
│  │  egui Frontend (genesight-gui)            │  │
│  │  - File selection (drag & drop)           │  │
│  │  - Report display (tier-based)            │  │
│  │  - Database management                    │  │
│  └──────────────────┬────────────────────────┘  │
│                     │ Direct function calls      │
│  ┌──────────────────┴────────────────────────┐  │
│  │  Rust Backend (genesight-core)            │  │
│  │  - Parser, Annotator, Scorer, Report      │  │
│  │  - SQLite DB Management                   │  │
│  │  - Optional: local LLM (Ollama)           │  │
│  └───────────────────────────────────────────┘  │
└─────────────────────────────────────────────────┘
```

The egui app calls into core directly -- no IPC boundary needed:

```rust
/// Analyzes a DNA file and returns a report.
fn analyze_dna(file_path: &str, db_path: &str) -> Result<Report, AppError> {
    let data = std::fs::read_to_string(file_path)?;
    let conn = Connection::open(db_path)?;

    let variants = genesight_core::parser::parse_auto(&data)?;
    let annotated = genesight_core::annotator::annotate_batch(&conn, &variants)?;
    let scored = genesight_core::scorer::score_all(&annotated)?;
    let report = genesight_core::report::generate(&scored, Format::Json)?;

    Ok(report)
}
```
