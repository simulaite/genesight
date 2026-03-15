# GeneSight – Open-Source DNA Analysis Tool

## Project Identity

**Name:** GeneSight (working title)
**Language:** German & English (code and API in English, documentation bilingual)
**License:** GPL-3.0-or-later (compatible with CC-BY-NC-SA 3.0 from SNPedia)
**Goal:** A privacy-first CLI and desktop tool that annotates personal DNA raw data (23andMe, AncestryDNA, VCF) against public genome databases and generates understandable reports — without data ever leaving the user's machine.

---

## Architecture Overview

```
genesight/
├── crates/
│   ├── genesight-core/       # Library Crate: Parser, Annotator, Scorer, Report Engine
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── parser/       # DNA file parsers (23andMe, AncestryDNA, VCF)
│   │   │   │   ├── mod.rs
│   │   │   │   ├── twentythreeandme.rs
│   │   │   │   ├── ancestry.rs
│   │   │   │   └── vcf.rs
│   │   │   ├── db/           # Database adapters (local SQLite)
│   │   │   │   ├── mod.rs
│   │   │   │   ├── clinvar.rs
│   │   │   │   ├── snpedia.rs
│   │   │   │   ├── gwas.rs
│   │   │   │   ├── dbsnp.rs
│   │   │   │   └── pharmgkb.rs
│   │   │   ├── annotator/    # Variant annotation against databases
│   │   │   │   ├── mod.rs
│   │   │   │   ├── clinical.rs    # ClinVar pathogenicity
│   │   │   │   ├── frequency.rs   # gnomAD/dbSNP allele frequencies
│   │   │   │   ├── pharmacogenomics.rs
│   │   │   │   └── traits.rs      # SNPedia traits & magnitude
│   │   │   ├── scorer/       # Risk scoring & confidence tiers
│   │   │   │   ├── mod.rs
│   │   │   │   ├── monogenic.rs   # Single-gene disorders (Tier 1: >95%)
│   │   │   │   ├── pharmaco.rs    # Pharmacogenetics (Tier 1: >95%)
│   │   │   │   ├── polygenic.rs   # Polygenic risk scores (Tier 2: 60-85%)
│   │   │   │   └── traits.rs      # Traits & lifestyle (Tier 2-3)
│   │   │   ├── report/       # Report generation
│   │   │   │   ├── mod.rs
│   │   │   │   ├── markdown.rs
│   │   │   │   ├── json.rs
│   │   │   │   └── html.rs
│   │   │   └── models/       # Shared types & structs
│   │   │       ├── mod.rs
│   │   │       ├── variant.rs
│   │   │       ├── annotation.rs
│   │   │       ├── confidence.rs  # ConfidenceTier enum
│   │   │       └── report.rs
│   │   └── Cargo.toml
│   ├── genesight-cli/        # CLI tool (clap)
│   │   ├── src/
│   │   │   └── main.rs
│   │   └── Cargo.toml
│   ├── genesight-server/     # Axum API (optional, for web version)
│   │   ├── src/
│   │   │   └── main.rs
│   │   └── Cargo.toml
│   └── genesight-desktop/    # Tauri App (Phase 2)
│       └── ...
├── data/
│   ├── fetch/                # Scripts for downloading the databases
│   │   ├── fetch_clinvar.sh
│   │   ├── fetch_snpedia.py
│   │   ├── fetch_gwas.sh
│   │   ├── fetch_dbsnp.sh
│   │   └── fetch_pharmgkb.sh
│   ├── import/               # Scripts for importing into SQLite
│   │   ├── import_clinvar.rs (or .py)
│   │   ├── import_snpedia.rs
│   │   └── import_gwas.rs
│   └── schema/               # SQLite schema definitions
│       └── schema.sql
├── tests/
│   ├── fixtures/             # Test DNA files (synthetic!)
│   │   ├── sample_23andme.txt
│   │   ├── sample_ancestry.txt
│   │   └── sample.vcf
│   └── integration/
├── docs/
│   ├── ARCHITECTURE.md
│   ├── DATA_SOURCES.md
│   ├── LICENSES.md
│   ├── CONFIDENCE_TIERS.md
│   └── CONTRIBUTING.md
├── Cargo.toml               # Workspace
├── CLAUDE.md                 # This file (Claude Code context)
├── LICENSE                   # GPL-3.0
└── README.md
```

---

## Data Sources & Licenses

### Primary Databases

| Database | Content | License | Access | Size (approx.) |
|----------|---------|---------|--------|-----------------|
| **ClinVar** | Clinically classified variants (pathogenic/benign), >3M variants | Public Domain (US Gov) | FTP: `ftp.ncbi.nlm.nih.gov/pub/clinvar/` + REST API | ~100MB (TSV) |
| **SNPedia** | Wiki with ~112K SNPs, magnitude scores, human-readable summaries | CC-BY-NC-SA 3.0 | MediaWiki API: `snpedia.com/w/api.php` | ~160MB (SQLite dump) |
| **GWAS Catalog** | Genome-wide association studies, polygenic traits | Open Access (EMBL-EBI) | REST API v2: `ebi.ac.uk/gwas/rest/api/v2/` + FTP | ~50MB |
| **dbSNP** | Reference SNP database (rs numbers, allele frequencies) | Public Domain (US Gov) | FTP: `ftp.ncbi.nih.gov/snp/` | ~15GB (complete), subset ~500MB |
| **gnomAD** | Allele frequencies from >250K genomes | Open Access | Download: `gnomad.broadinstitute.org` | Multi-GB, subset ~1GB |
| **PharmGKB** | Pharmacogenetics (drug-gene interactions) | CC-BY-SA 4.0 (academically free) | Download + API: `pharmgkb.org` | ~50MB |

### License Compatibility

- **GPL-3.0** (our project) is compatible with:
  - Public Domain (ClinVar, dbSNP) ✅
  - CC-BY-NC-SA 3.0 (SNPedia) ✅ — as long as we remain non-commercial or treat SNPedia data as a separate, optionally downloadable dataset
  - CC-BY-SA 4.0 (PharmGKB) ✅
  - Open Access (GWAS Catalog, gnomAD) ✅

- **Important:** CC-BY-NC-SA 3.0 from SNPedia means:
  - ✅ Open-source project: no problem
  - ✅ Personal/academic use: no problem
  - ⚠️ If someone wants to commercially fork the project: SNPedia data must be removed or separately licensed
  - → **Solution:** Treat SNPedia data as an optional download, not bundled in the repo

### Attribution Requirements

Every use must correctly attribute:
- ClinVar: "ClinVar data provided by NCBI (National Center for Biotechnology Information)"
- SNPedia: "SNPedia content is licensed under CC-BY-NC-SA 3.0 by SNPedia.com"
- GWAS Catalog: "GWAS Catalog provided by NHGRI-EBI"
- PharmGKB: "PharmGKB data © PharmGKB, CC-BY-SA 4.0"

---

## Confidence Tier System

All results are categorized into three reliability levels:

### Tier 1: Reliable (>95% Accuracy)
- **Monogenic Disorders** — A single variant is directly causal (e.g., BRCA1/2, CFTR, Huntington)
- **Carrier Status** — Carrier status for recessive disorders
- **Pharmacogenetics** — Drug metabolism (CYP2D6, CYP2C19, etc.)
- **Simple Traits** — Lactose tolerance, earwax type, etc.
- Source: Primarily ClinVar (review status ≥ 2 stars), PharmGKB (Level 1-2)

### Tier 2: Probable (60-85% Accuracy)
- **Polygenic Risk Scores** — Diabetes, heart disease, hypertension
- **Physical Traits** — Hair color, freckles, baldness risk
- Source: GWAS Catalog, SNPedia (Magnitude ≥ 3)

### Tier 3: Speculative (50-65% Accuracy)
- **Complex Disorders** — Depression, schizophrenia, autism
- **Personality Traits** — Intelligence, risk-taking propensity
- **Athletic Aptitude** — ACTN3, VO2max predisposition
- Source: GWAS Catalog (low effect size), SNPedia (Magnitude < 3)

**Rule:** Every result MUST be assigned a `ConfidenceTier`. The report displays this prominently.

---

## Supported Input Formats

### 23andMe Raw Data
```
# rsid  chromosome  position  genotype
rs4477212  1  82154  AA
rs3094315  1  752566  AG
```
- Tab-separated, comment lines start with `#`
- Header line: `rsid  chromosome  position  genotype`
- Genotype: 2 characters (e.g., AA, AG, CT), `--` for no-call, `I` or `D` for indels

### AncestryDNA Raw Data
```
rsid  chromosome  position  allele1  allele2
rs4477212  1  82154  A  A
```
- Tab-separated, comment lines start with `#`
- Alleles are separated into two columns

### VCF (Variant Call Format)
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
1  82154  rs4477212  G  A  .  PASS  .  GT  0/1
```
- Standard bioinformatics format
- More complex to parse, but most comprehensive

---

## Development Phases

### Phase 1: Data & CLI (CURRENT)
1. **Data Fetching Scripts** — Download all databases
2. **SQLite Import** — Import data into local, queryable database
3. **DNA Parsers** — Read 23andMe, AncestryDNA, VCF files
4. **Annotation Engine** — Match variants against local DB
5. **CLI Tool** — `genesight analyze my_dna.txt --format markdown`
6. **Report Generator** — Markdown/JSON/HTML output with confidence tiers

### Phase 2: Desktop App
7. **Tauri Integration** — GUI around the core
8. **Auto-Update** — Database updates in the background
9. **LLM Integration** — Optional: Summarize results in plain language via LLM

### Phase 3: Web & Community
10. **Axum API** — For web version (with explicit privacy disclaimer)
11. **Community Reports** — Anonymized, aggregated statistics

---

## Coding Conventions

### Rust
- **Edition:** 2021
- **MSRV:** 1.75+
- **Error Handling:** `thiserror` for library errors, `anyhow` for CLI/App
- **Async:** `tokio` (for data fetching and server), sync for core logic
- **Serialization:** `serde` + `serde_json`
- **CLI:** `clap` v4 (derive API)
- **Database:** `rusqlite` (with bundled SQLite)
- **HTTP Client:** `reqwest` (for data fetching)
- **Testing:** Unit tests in every module, integration tests in `tests/`

### Code Style
- `cargo fmt` and `cargo clippy` must pass cleanly
- All public functions have doc comments
- No `unwrap()` in library code — only in tests and CLI with context
- English variable and function names
- German comments are OK, doc comments in English

### Git
- Conventional Commits: `feat:`, `fix:`, `docs:`, `data:`, `refactor:`
- Branch scheme: `feat/parser-23andme`, `data/clinvar-import`
- No raw DNA data in the repo — only synthetic test data

---

## Important Rules

1. **No real DNA data in the repository.** Test data must be synthetically generated.
2. **No medical diagnoses.** The report is informational, not diagnostic. Every report contains a disclaimer.
3. **Privacy first.** No telemetry, no data uploads, no analytics. Everything local.
4. **Confidence tiers are mandatory.** No result without an assigned reliability level.
5. **Attributions are mandatory.** Every data source must be correctly attributed in the report.
6. **Offline-capable.** After the initial database download, the tool must function completely offline.

---

## Current Focus: Phase 1 – Fetch Data & CLI

### Task 1: Create Data Fetching Scripts
- `data/fetch/fetch_clinvar.sh` — ClinVar VCF + variant_summary.txt from NCBI FTP
- `data/fetch/fetch_snpedia.py` — Scrape SNPedia via MediaWiki API (respect rate limits: 3s delay)
- `data/fetch/fetch_gwas.sh` — GWAS Catalog TSV download
- `data/fetch/fetch_dbsnp.sh` — dbSNP relevant subset data

### Task 2: SQLite Schema & Import
- Unified schema in `data/schema/schema.sql`
- Import scripts that transform downloaded data into SQLite
- Goal: A single `genesight.db` file (~500MB-1GB)

### Task 3: DNA Parsers
- `genesight-core` parsers for 23andMe, AncestryDNA, VCF
- Result: `Vec<Variant>` with rsID, chromosome, position, genotype

### Task 4: Annotation & CLI
- Lookup each variant against the local SQLite
- CLI interface: `genesight analyze <file> [--format json|md|html] [--tier 1|2|3]`
- Report output with confidence tiers
