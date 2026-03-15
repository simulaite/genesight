# GeneSight Scientific Audit Report

**Date:** 2026-03-15
**Method:** Automated code analysis against `docs/research/deep-research-report.md` + live test with PGP individual huAE4518
**Test file:** `tests/fixtures/pgp/huAE4518_23andme_v4.txt` (601,783 variants)

## Live Test Summary

| Metric | Value |
|--------|-------|
| Variants parsed | 601,783 |
| Annotated | 30 |
| Scored results | 17 (9 Tier1, 1 Tier2, 7 Tier3) |
| PGx false positives | **6 of 8** |
| ClinVar suspicious | **2 of 4** (likely strand-flip FPs) |
| GWAS correct | 7 of 7 |
| Tests passing | 255/255 + 4 doc-tests |

## Issue Index

### Critical (Producing False Results Now)

| ID | Issue | File | Impact |
|----|-------|------|--------|
| [C1](critical/C1-pgx-no-allele-check.md) | PGx scorer has no allele check | `scorer/mod.rs:469-510` | 6 of 8 PGx results are false positives in live test |
| [C2](critical/C2-pgx-pipeline-dead-code.md) | New PGx pipeline (diplotype/phasing/coverage) is dead code | `lib.rs:208-213` | All correct PGx logic unreachable from `analyze()` |
| [C3](critical/C3-clinvar-star-mapping.md) | ClinVar import maps review stars off-by-one | `import_clinvar.py:38-47` | Single-submitter entries promoted to Tier1 |
| [C4](critical/C4-pgx-schema-mismatch.md) | `definitions.rs` queries non-existent column names | `pgx/definitions.rs:63-66` | New PGx pipeline crashes on real DB |

### High (Missing Safety Gates)

| ID | Issue | File | Impact |
|----|-------|------|--------|
| [H1](high/H1-clinvar-no-inheritance-mode.md) | No mode-of-inheritance (AD/AR) in ClinVar scoring | `scorer/mod.rs` | Het carrier of AR disease reported as MonogenicDisease |
| [H2](high/H2-clinvar-germline-import.md) | Import never populates `classification_type` | `import_clinvar.py` | Germline/somatic distinction is dead feature |
| [H3](high/H3-confirmation-urgency-not-rendered.md) | `ConfirmationUrgency` computed but never rendered | `report/html.rs`, `report/markdown.rs` | BRCA1 looks identical to low-mag SNPedia hit |
| [H4](high/H4-fda-pgx-disclaimer.md) | No FDA PGx disclaimer anywhere | entire codebase | Missing required safety language |
| [H5](high/H5-gwas-or-interpretation.md) | OR presented without relative-risk explanation | `scorer/mod.rs:725-750` | Users may interpret OR as absolute probability |
| [H6](high/H6-pgx-strand-normalization.md) | PGx star allele caller uses raw genotype chars | `pgx/mod.rs:150` | Opposite-strand data produces wrong allele counts |
| [H7](high/H7-prs-not-implemented.md) | Polygenic Risk Scores not implemented | `scorer/polygenic.rs` | Single-SNP GWAS hits mislabeled as "PRS" |
| [H8](high/H8-assembly-no-enforcement.md) | Assembly mismatch warns but doesn't block | `lib.rs:260-293` | Pipeline continues with potentially wrong results |

### Medium (Missing Features / Quality)

| ID | Issue | File | Impact |
|----|-------|------|--------|
| [M1](medium/M1-gwas-or-inversion.md) | No historical GWAS OR inversion detection | `scorer/mod.rs` | Pre-2021 inverted ORs silently accepted |
| [M2](medium/M2-dtc-context-not-rendered.md) | `Report.dtc_context` populated but not rendered | `report/html.rs`, `report/markdown.rs` | DTC context statement invisible to user |
| [M3](medium/M3-clinvar-allele-fallback.md) | ClinVar allele check falls back to rsID-only | `scorer/mod.rs:255` | Missing variants table → no allele verification |

### Test Results

| Doc | Description |
|-----|-------------|
| [huAE4518 Analysis](test-results/huAE4518-analysis.md) | Full live test results with the PGP sample |

## Cross-Reference to Research

| Research Section | Issues |
|-----------------|--------|
| Allele Matching & Strand Orientation | C1, H6, M3 |
| PGx Star-Allele Calling & Phenotyping | C2, C4, H6 |
| ClinVar Pathogenicity Interpretation | C3, H1, H2 |
| GWAS Risk Interpretation & PRS | H5, H7, M1 |
| Responsible Result Presentation | H3, H4, M2 |
| Population Frequencies | (BA1 implemented; popmax correct) |
