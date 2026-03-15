# GeneSight

**Privacy-first, open-source DNA analysis tool.**

Analyze your raw DNA data (23andMe, AncestryDNA, VCF) against public genomic databases -- entirely offline, on your own machine. Your DNA never leaves your computer.

> **Medical Disclaimer:** GeneSight is not a medical device and does not provide medical advice. The reports generated are for informational and educational purposes only. Genetic variants are interpreted based on publicly available research databases, which may contain errors, incomplete data, or findings that have not been clinically validated. Do not make medical decisions based on GeneSight output. Always consult a qualified physician or certified genetic counselor before acting on any genetic information.

---

*This project is in active development (v0.1.0). APIs and report formats may change.*

## What It Does

GeneSight reads your raw DNA file, matches your variants against curated scientific databases, and generates a structured report. Every finding is assigned a confidence tier so you know how much weight to give it:

- **Tier 1 -- Reliable (>95%):** Well-established monogenic conditions, pharmacogenomic interactions (e.g., CYP2D6 metabolizer status), and carrier screening results backed by extensive clinical evidence.
- **Tier 2 -- Probable (60--85%):** Polygenic risk scores aggregated from genome-wide association studies, and physical traits with moderate effect sizes.
- **Tier 3 -- Speculative (50--65%):** Complex traits influenced by many genes and environmental factors, where current research provides only weak predictive power.

## Key Features

- **Multi-format parsing** -- Reads 23andMe, AncestryDNA, and VCF files with automatic format detection
- **ClinVar pathogenicity** -- Matches variants against NCBI ClinVar for clinical significance classifications
- **Pharmacogenomics** -- Star-allele calling and metabolizer phenotype prediction for CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, SLCO1B1, and VKORC1
- **GWAS polygenic scores** -- Aggregated risk scores from the NHGRI-EBI GWAS Catalog
- **Allele frequencies** -- Population frequency context from gnomAD (250K+ genomes)
- **Confidence-tiered reports** -- Output in Markdown, JSON, or HTML with every result classified by evidence strength
- **TUI interactive mode** -- Terminal-based interface for browsing results interactively
- **Desktop GUI** -- Native egui-based desktop application for point-and-click analysis
- **Fully offline** -- After initial database download, no network access is required or attempted

## How It Works

```
                         +-------------+
                         |  SQLite DB  |
                         |  (ClinVar,  |
                         |   GWAS,     |
                         |   gnomAD,   |
                         |   PharmGKB) |
                         +------+------+
                                |
  +----------+    +--------+    |    +---------+    +----------+
  | DNA File +--->| Parser +--->+---->Annotator+--->| Scorer   |
  | (raw)    |    |        |        |         |    |          |
  +----------+    +--------+        +---------+    +----+-----+
                                                        |
                                                   +----v-----+
                                                   |  Report  |
                                                   | (MD/JSON |
                                                   |  /HTML)  |
                                                   +----------+
```

1. **Parser** -- Detects the input format and extracts genotype calls (rsID, chromosome, position, alleles).
2. **Annotator** -- Batch-queries the local SQLite database to attach clinical significance, allele frequencies, drug interactions, and GWAS associations to each variant.
3. **Scorer** -- Aggregates annotations into findings, computes polygenic risk scores, resolves pharmacogenomic diplotypes, and assigns a confidence tier to every result.
4. **Report** -- Renders the scored findings into a human-readable report with the mandatory medical disclaimer and data source attributions.

For a detailed architecture overview and data flow diagrams, see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Quick Start

```bash
# Build from source
git clone https://github.com/simulaite/genesight.git
cd genesight
cargo build --release

# Download reference databases (~500 MB)
./target/release/genesight fetch --all

# Analyze your DNA file
./target/release/genesight analyze my_23andme_data.txt --format markdown

# Or use the interactive TUI
./target/release/genesight analyze my_23andme_data.txt --tui

# Launch the desktop GUI
cargo run -p genesight-gui --release
```

## Supported Input Formats

| Format | File Extension | Detection | Notes |
|--------|---------------|-----------|-------|
| 23andMe | `.txt` | Automatic (header detection) | v5 format tested |
| AncestryDNA | `.txt` | Automatic (header detection) | Standard export format |
| VCF | `.vcf`, `.vcf.gz` | Automatic (magic bytes / header) | VCF 4.x supported |

Format detection is fully automatic. GeneSight inspects file headers and structure to determine the correct parser -- no manual format flag is needed.

## Data Sources

| Database | Content | License | Status |
|----------|---------|---------|--------|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinically classified variant pathogenicity | Public Domain | Included |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | Genome-wide association study results | Open Access | Included |
| [gnomAD](https://gnomad.broadinstitute.org/) | Allele frequencies from 250K+ genomes | ODC-ODbL | Included |
| [PharmGKB](https://www.pharmgkb.org/) | Pharmacogenomic annotations and guidelines | CC-BY-SA 4.0 | Included |
| [SNPedia](https://www.snpedia.com/) | Community-curated variant descriptions | CC-BY-NC-SA 3.0 | Optional, separate download |

SNPedia data is kept in a separate database (`snpedia.db`) due to its non-commercial license. Download it with `genesight fetch --snpedia`.

For download instructions, update schedules, and schema details, see [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md).

## Confidence Tier System

Every result in a GeneSight report is assigned one of three confidence tiers:

| Tier | Confidence | Typical Use Cases | Example |
|------|-----------|-------------------|---------|
| **Tier 1 -- Reliable** | >95% | Monogenic conditions, pharmacogenomics, carrier status | CYP2D6 poor metabolizer, BRCA1 pathogenic variant |
| **Tier 2 -- Probable** | 60--85% | Polygenic risk scores, well-studied traits | Type 2 diabetes PRS, lactose intolerance |
| **Tier 3 -- Speculative** | 50--65% | Complex traits, weakly predictive associations | Morning chronotype tendency, muscle fiber composition |

Tiers are determined by the strength of the underlying evidence: number of replicated studies, effect sizes, clinical validation status, and database agreement. See [docs/CONFIDENCE_TIERS.md](docs/CONFIDENCE_TIERS.md) for the full classification methodology.

## Privacy

- **No telemetry.** No analytics, no usage tracking, no crash reporting.
- **No data upload.** Your DNA file is processed entirely on your local machine and is never transmitted anywhere.
- **No account required.** No registration, no email, no login.
- **Fully inspectable.** The entire codebase is open source under GPL-3.0 -- verify the privacy guarantees yourself.
- **Zero network dependencies in core.** The `genesight-core` library crate has no network dependencies whatsoever. Network access exists only in the CLI for the initial database download (`fetch` command).

## Building from Source

Requires Rust 1.75 or later.

```bash
# Build all crates
cargo build

# Run the test suite
cargo test

# Lint
cargo clippy

# Check formatting
cargo fmt -- --check
```

## Project Structure

```
genesight/
├── crates/
│   ├── genesight-core/       # Library: parsers, DB adapters, annotators, scorers, reports
│   ├── genesight-cli/        # Binary: CLI tool with TUI mode (clap)
│   ├── genesight-gui/        # Binary: Desktop GUI application (egui)
│   └── genesight-server/     # Binary: Web API (Axum, planned)
├── data/
│   ├── fetch/                # Database download scripts
│   ├── import/               # Data transformation and import scripts
│   ├── schema/               # SQLite schema definitions
│   └── seed/                 # Seed database builder
├── tests/
│   └── fixtures/             # Synthetic test data (no real DNA)
└── docs/                     # Project documentation
```

The core library (`genesight-core`) performs no filesystem I/O and no network access. It accepts `&str`, `&[u8]`, and `rusqlite::Connection` as parameters, making it straightforward to embed in other applications.

## Documentation

- [Architecture](docs/ARCHITECTURE.md) -- System design, crate responsibilities, and data flow
- [Data Sources](docs/DATA_SOURCES.md) -- Database details, download instructions, and update schedule
- [Confidence Tiers](docs/CONFIDENCE_TIERS.md) -- Tier classification methodology and evidence thresholds
- [Licenses](docs/LICENSES.md) -- License compatibility analysis for all dependencies and data sources
- [PGP Test Data](docs/PGP_TEST_DATA.md) -- Public genome project data used for testing
- [Audit Trail](docs/audit/) -- Development audit logs

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on submitting issues, feature requests, and pull requests.

## License

GeneSight is licensed under the [GNU General Public License v3.0 or later](LICENSE).

```
SPDX-License-Identifier: GPL-3.0-or-later
```

**Note on SNPedia data:** The optional SNPedia database is licensed under CC-BY-NC-SA 3.0, which restricts commercial use. It is distributed separately from the main application and is not required for GeneSight to function.

## Acknowledgments

GeneSight relies on the work of the following public databases and research institutions:

- **ClinVar** -- National Center for Biotechnology Information (NCBI), National Library of Medicine
- **GWAS Catalog** -- National Human Genome Research Institute (NHGRI) and European Bioinformatics Institute (EMBL-EBI)
- **gnomAD** -- Genome Aggregation Database, Broad Institute of MIT and Harvard
- **PharmGKB** -- Pharmacogenomics Knowledgebase, Stanford University
- **SNPedia** -- Community-curated wiki of human genetic variants

This tool is built on top of decades of publicly funded genomics research. We are grateful to the researchers, institutions, and volunteers who make these resources freely available.

## Contact

**STONKS GmbH**
Buber-Neumann-Weg 68
60439 Frankfurt am Main

E-Mail: [info@simulaite.ai](mailto:info@simulaite.ai)
