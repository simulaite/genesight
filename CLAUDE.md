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

Rust workspace with three crates:
- `genesight-core` — Library: **no filesystem IO, no network**. Takes `&str`/`&[u8]` and
  `rusqlite::Connection` as parameters. Parsers, DB adapters, annotators, scorers, report engine.
- `genesight-cli` — Binary: `genesight` CLI tool using clap. Handles IO, opens files/DBs.
- `genesight-server` — Binary: Axum web API (Phase 3)

See `docs/ARCHITECTURE.md` for full details including data flow diagrams.

## Key Conventions

- Error handling: `thiserror` in library, `anyhow` in binaries
- No `unwrap()` in library code
- All public functions have English doc-comments
- Every analysis result must have a `ConfidenceTier` (Tier1Reliable/Tier2Probable/Tier3Speculative)
- No real DNA data in the repository — only synthetic test fixtures
- Data source attributions are mandatory in reports
- Medical disclaimer is mandatory in reports

## Database Architecture

Two separate databases:
- `genesight.db` — Main database (ClinVar, GWAS, dbSNP/gnomAD frequencies, PharmGKB). All permissive licenses.
- `snpedia.db` — Optional, separate database (CC-BY-NC-SA 3.0). Downloaded separately via `genesight fetch --snpedia`.

Schema defined in `data/schema/schema.sql`. See `docs/DATA_SOURCES.md` and `docs/LICENSES.md`.

## Data Pipeline

1. Fetch scripts (`data/fetch/`) download public databases
2. Import scripts transform data into local SQLite
3. Core library queries SQLite via `rusqlite` (bundled)
4. Batch queries via temp tables for performance (not individual queries per variant)

## Current Phase: Phase 1 (Data & CLI)

Focus: parsers, database import, annotation engine, CLI interface.
See `docs/CLAUDE.md` for full project spec.
