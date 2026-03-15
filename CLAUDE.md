# GeneSight – Claude Code Context

## Project Overview

GeneSight is a privacy-first CLI and desktop tool for annotating personal DNA raw data
(23andMe, AncestryDNA, VCF) against public genome databases. All processing happens
locally — no data ever leaves the user's machine.

## Build & Test

```bash
cargo build                    # Build all crates
cargo test                     # Run all tests
cargo run -p genesight-cli -- analyze tests/fixtures/sample_23andme.txt
cargo clippy                   # Lint
cargo fmt -- --check           # Format check
```

## Architecture

Rust workspace with four crates:
- `genesight-core` — Library: **no filesystem IO, no network**. Takes `&str`/`&[u8]` and
  `rusqlite::Connection` as parameters. Modules: parser, db, annotator, scorer, pgx, allele,
  normalizer, report, models.
- `genesight-cli` — Binary: `genesight` CLI tool using clap. Handles IO, opens files/DBs.
  Includes TUI for interactive use.
- `genesight-gui` — Binary: Desktop GUI application (eframe/egui).
- `genesight-server` — Binary: Axum web API (Phase 3, not yet active)

See `docs/ARCHITECTURE.md` for full details including data flow diagrams.

## Key Conventions

- Error handling: `thiserror` in library, `anyhow` in binaries
- No `unwrap()` in library code
- All public functions have English doc-comments
- Every analysis result must have a `ConfidenceTier` (Tier1Reliable/Tier2Probable/Tier3Speculative)
- Every result carries a `ConfirmationUrgency` (HighImpact/ClinicalConfirmationRecommended/InformationalOnly)
- No real DNA data in the repository — only synthetic test fixtures
- Data source attributions are mandatory in reports
- Medical disclaimer and FDA PGx disclaimer are mandatory in reports
- ClinVar `classification_type` (germline/somatic/oncogenicity) must be respected — non-germline demoted to Tier3
- Autosomal recessive carrier detection: het pathogenic in AR gene → `CarrierStatus`, not `MonogenicDisease`
- GWAS single-SNP hits use `GwasAssociation` category (not `PolygenicRiskScore` — true PRS deferred to Phase 2)
- Allele verification required before scoring — rsID-only fallback adds limitations and tier downgrade
- Strand-aware allele matching via `allele` module (complement checking for all scoring paths)
- Genome assembly tracking: input/database assembly compatibility checked, warnings added per-result on mismatch

## Database Architecture

Two separate databases:
- `genesight.db` — Main database (ClinVar, GWAS, dbSNP/gnomAD frequencies, PharmGKB, variants). All permissive licenses.
- `snpedia.db` — Optional, separate database (CC-BY-NC-SA 3.0). Downloaded separately via `genesight fetch --snpedia`.

Schema defined in `data/schema/schema.sql`. See `docs/DATA_SOURCES.md` and `docs/LICENSES.md`.

## Data Pipeline

1. Fetch scripts (`data/fetch/`) download public databases
2. Import scripts (`data/import/`) transform data into local SQLite
3. Core library queries SQLite via `rusqlite` (bundled)
4. Batch queries via temp tables for performance (not individual queries per variant)
5. PGx star allele pipeline runs separately via `pgx::StarAlleleCaller`

## Current Phase: Phase 1 (Data & CLI) — Nearing Completion

Phase 1 focus: parsers, database import, annotation engine, scoring, CLI interface.

### Completed
- DNA parsers (23andMe, AncestryDNA, VCF) with assembly detection
- ClinVar import with correct star mapping, classification_type, and allele-aware scoring
- GWAS annotation with OR/beta caveats and OR inversion warnings
- PGx star allele calling (integrated into main pipeline, strand-normalized)
- Confidence tier system with DTC caveats and urgency levels
- Report generation (HTML, Markdown, JSON) with FDA disclaimers, DTC context, urgency banners
- Assembly mismatch detection and per-result warnings
- AR carrier detection (AD/AR mode-of-inheritance)
- Desktop GUI (eframe/egui)

### Open Issues
- **#2** (C2): PGx diplotype/phasing/coverage modules still dead code — `StarAlleleCaller` works but
  `diplotype.rs`, `phasing.rs`, and coverage-aware phenotyping are not integrated
- **#13** (M1): GWAS OR inversion caveat applied unconditionally — `study_date` schema field needed
  for targeted pre-2021 warnings

### Phase 2 (planned)
- True polygenic risk score implementation (`scorer/polygenic.rs` is a stub)
- PGx diplotype calling with phase ambiguity detection and coverage tracking
- Assembly LiftOver (GRCh37 ↔ GRCh38) with `--force-assembly` blocking
- Tauri desktop app / auto-update

See `docs/CLAUDE.md` for full project spec.
