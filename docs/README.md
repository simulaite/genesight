# 🧬 GeneSight

**Open-source, privacy-first DNA analysis tool.**

Analyze your raw DNA data (23andMe, AncestryDNA, VCF) against public genomic databases — entirely offline, on your own machine. Your DNA never leaves your computer.

> ⚠️ **This is not a medical device.** GeneSight provides informational reports only. Consult a physician or genetic counselor for medical decisions.

## What it does

GeneSight reads your raw DNA file, matches your variants against curated scientific databases, and generates a report with clear confidence tiers:

- 🟢 **Tier 1 — Reliable** (>95%): Monogenic diseases, carrier status, pharmacogenomics
- 🟡 **Tier 2 — Probable** (60-85%): Polygenic risk scores, physical traits
- 🔴 **Tier 3 — Speculative** (50-65%): Complex traits, personality, athletic ability

## Quick Start

```bash
# 1. Install
cargo install genesight-cli

# 2. Download reference databases (~500MB)
genesight fetch --all

# 3. Analyze your DNA
genesight analyze my_23andme_data.txt --format markdown
```

## Data Sources

| Database | Content | License |
|----------|---------|---------|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinically classified variants | Public Domain |
| [SNPedia](https://www.snpedia.com/) | Wiki with human-readable variant descriptions | CC-BY-NC-SA 3.0 |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | Genome-wide association studies | Open Access |
| [gnomAD](https://gnomad.broadinstitute.org/) | Allele frequencies from 250K+ genomes | ODC-ODbL |
| [PharmGKB](https://www.pharmgkb.org/) | Pharmacogenomics | CC-BY-SA 4.0 |

See [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md) for full details.

## Privacy

- **No telemetry.** No analytics. No network calls after database download.
- **No data upload.** Your DNA file is processed locally and never transmitted.
- **No account required.** No registration, no email, no tracking.
- **Inspectable.** The code is GPL-3.0 — verify it yourself.

## Documentation

- [Architecture](docs/ARCHITECTURE.md) — System design and data flow
- [Data Sources](docs/DATA_SOURCES.md) — Databases, APIs, and download details
- [Confidence Tiers](docs/CONFIDENCE_TIERS.md) — How results are classified
- [Licenses](docs/LICENSES.md) — License compatibility analysis

## License

GeneSight is licensed under [GPL-3.0-or-later](LICENSE).

Note: SNPedia data (optional) is licensed under CC-BY-NC-SA 3.0 and distributed separately.

## Contributing

See [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) — contributions welcome!
