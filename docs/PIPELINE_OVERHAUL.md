# GeneSight Pipeline Scientific Overhaul

## Overview

This document captures the complete architecture blueprints for fixing 13 scientific deficiencies
identified by comparing the implemented pipeline against the deep research report
(`docs/research/deep-research-report.md`). Each blueprint was designed by a specialized
architecture agent with full codebase access.

**Date:** 2026-03-15
**Status:** Planning complete, implementation starting

---

## Table of Contents

1. [Deficiency Summary](#deficiency-summary)
2. [Implementation Phases](#implementation-phases)
3. [Blueprint 1: PGx Missing-Data-as-\*1](#blueprint-1-pgx-missing-data-as-1)
4. [Blueprint 2: Palindromic SNP Ambiguity](#blueprint-2-palindromic-snp-ambiguity)
5. [Blueprint 3: GWAS Strand & Effect Context](#blueprint-3-gwas-strand--effect-context)
6. [Blueprint 4: ClinVar 0-Star & Conflicting Classifications](#blueprint-4-clinvar-0-star--conflicting-classifications)
7. [Blueprint 5: ClinVar Germline/Somatic Classification](#blueprint-5-clinvar-germlinesomatic-classification)
8. [Blueprint 6: Assembly/Build Tracking](#blueprint-6-assemblybuild-tracking)
9. [Blueprint 7: Universal Clinical Confirmation Language](#blueprint-7-universal-clinical-confirmation-language)
10. [Blueprint 8: VCF Normalization](#blueprint-8-vcf-normalization)
11. [Blueprint 9: PGx CPIC Data Infrastructure](#blueprint-9-pgx-cpic-data-infrastructure)
12. [Blueprint 10: PGx Phasing & Ambiguity](#blueprint-10-pgx-phasing--ambiguity)
13. [Total Impact Summary](#total-impact-summary)

---

## Deficiency Summary

| # | Deficiency | Severity | Blueprint |
|---|-----------|----------|-----------|
| 1 | PGx missing data silently treated as \*1 (wildtype) | CRITICAL | 1 |
| 2 | Palindromic SNPs (A/T, C/G) not flagged as strand-ambiguous | CRITICAL | 2 |
| 3 | GWAS risk allele not strand-normalized, no effect context | HIGH | 3 |
| 4 | ClinVar 0-star entries treated as clinically valid | HIGH | 4 |
| 5 | ClinVar conflicting classifications silently dropped | HIGH | 4 |
| 6 | ClinVar germline/somatic/oncogenicity not distinguished | HIGH | 5 |
| 7 | No genome assembly (GRCh37/38) tracking | MEDIUM | 6 |
| 8 | No universal DTC confirmation language | HIGH | 7 |
| 9 | VCF multi-allelic records not split | MEDIUM | 8 |
| 10 | VCF alleles not trimmed/left-aligned | MEDIUM | 8 |
| 11 | PGx covers only 4 genes, needs CPIC expansion | HIGH | 9 |
| 12 | PGx phasing ambiguity not detected (TPMT\*3A) | HIGH | 10 |
| 13 | BA1 threshold not disease-specific | MEDIUM | (deferred) |

---

## Implementation Phases

```
Phase 1 (can start immediately, in parallel):
  +-- Blueprint 2  (Palindromic SNP fix)
  +-- Blueprint 6  (Assembly tracking)
  +-- Blueprint 8  (VCF normalization)
  +-- Blueprint 9  (CPIC data expansion)

Phase 2 (depends on Blueprint 2 for allele matching):
  +-- Blueprint 1  (PGx missing-data)
  +-- Blueprint 3  (GWAS strand)
  +-- Blueprint 4  (ClinVar quality gates)
  +-- Blueprint 5  (ClinVar germline/somatic)

Phase 3 (depends on Blueprint 1 + 9):
  +-- Blueprint 10 (PGx phasing)

Phase 4 (depends on all prior):
  +-- Blueprint 7  (Clinical confirmation language)
```

---

## Blueprint 1: PGx Missing-Data-as-\*1

**Severity:** CRITICAL
**Problem:** When a defining SNP for a star allele is absent from the user's array data,
the pipeline treats this as "reference allele observed" and defaults to \*1/\*1. This can
produce dangerous false Normal/Rapid/Ultrarapid metabolizer calls. CPIC's #1 pitfall.

### New Types

```rust
pub enum CoverageStatus {
    Complete,                    // All defining positions observed
    Partial { missing: Vec<String>, coverage_pct: f64 },
    Insufficient { missing: Vec<String> },
}

// DEFAULT_MINIMUM_COVERAGE = 1.0 (all positions required)
```

### Changes

| File | Change |
|------|--------|
| `pgx/diplotype.rs` | Replace `is_complete: bool` with `CoverageStatus` on `DiplotypeCall` |
| `pgx/phenotype.rs` | Guard in `call_phenotype()`: return "Indeterminate" when `Insufficient` |
| `pgx/mod.rs` | Pass coverage status through `PgxResult` |
| `report/html.rs` | Render coverage warning badge |
| `report/markdown.rs` | Render coverage warning |

### Test Cases (7)
1. All positions present -> `CoverageStatus::Complete`, normal phenotype call
2. One of two positions missing -> `Partial`, phenotype still called with limitation
3. All positions missing -> `Insufficient`, phenotype = "Indeterminate"
4. Single-SNP gene, position present -> `Complete`
5. Single-SNP gene, position missing -> `Insufficient`
6. CYP2D6 with \*4 defining position absent -> verify limitation mentions CYP2D6 CNV
7. Report renders "Indeterminate" with explanation when `Insufficient`

---

## Blueprint 2: Palindromic SNP Ambiguity

**Severity:** CRITICAL
**Problem:** `match_alleles()` at `allele/mod.rs:64-77` returns a direct match for palindromic
SNPs (A/T, C/G) instead of `StrandAmbiguous`. The frequency-aware resolver
`match_alleles_with_frequency()` exists but is dead code (never called).

### Changes

| File | Change |
|------|--------|
| `allele/mod.rs` | Fix lines 64-77: return `StrandAmbiguous` unconditionally for A/T, C/G |
| `scorer/mod.rs` | Wire `match_alleles_with_frequency()` at lines 62-67 using `db_af` from `av.frequency` |

### Algorithm

```
1. Check if SNP is palindromic: is_palindromic(ref, alt)
   - A/T or T/A -> true
   - C/G or G/C -> true
   - All others -> false

2. If palindromic:
   a. Try frequency-based resolution using try_resolve_palindromic(user_af, db_af)
   b. If resolution succeeds -> return resolved AlleleMatch
   c. If resolution fails -> return StrandAmbiguous

3. If not palindromic: existing logic (complement if needed, match directly)
```

### Test Cases (6)
1. Non-palindromic A/G SNP -> direct match, not ambiguous
2. Palindromic A/T SNP without frequency data -> `StrandAmbiguous`
3. Palindromic C/G SNP with frequency data resolving strand -> correct match
4. Palindromic A/T SNP with frequency data but ambiguous AF -> `StrandAmbiguous`
5. Scorer skips `StrandAmbiguous` variants with caveat in report
6. Integration: palindromic ClinVar variant flagged correctly

---

## Blueprint 3: GWAS Strand & Effect Context

**Severity:** HIGH
**Problem:** `count_risk_allele_copies()` at `scorer/mod.rs:463-482` does naive char comparison
without strand normalization. No effect context (OR classification, ancestry caveat, historical
OR inversion note) is provided.

### New Types

```rust
pub enum RiskAlleleCopies {
    Known { copies: u8, via_complement: bool },
    Ambiguous,       // Palindromic SNP, can't resolve
    Indeterminate,   // No risk allele specified in GWAS
}
```

### Changes

| File | Change |
|------|--------|
| `scorer/mod.rs` | Replace `count_risk_allele_copies` with strand-aware version using ref/alt and AlleleMatch |
| `scorer/mod.rs` | New `gwas_effect_context()` adding OR/beta classification, ancestry caveat |
| `db/gwas.rs` | Add `study_date` and `ancestry_cohort` fields (future schema) |

### Effect Context Rules
- OR > 2.0: "Substantially elevated risk"
- OR 1.5-2.0: "Moderately elevated risk"
- OR 1.2-1.5: "Slightly elevated risk"
- OR < 1.0: Historical inversion note ("Pre-2021 GWAS may report OR<1 for protective alleles")
- Single-SNP caveat: "This association is based on a single variant"
- Ancestry caveat: "Effect size may vary across populations"

### Test Cases (11)
1. Risk allele matches user's alt allele directly -> `Known { copies: 1 }`
2. Risk allele matches complement of user's allele -> `Known { copies: 1, via_complement: true }`
3. Palindromic risk allele -> `Ambiguous`
4. No risk allele in GWAS -> `Indeterminate`
5. User is homozygous ref -> `Known { copies: 0 }` -> SKIP
6. User has 2 copies of risk allele -> `Known { copies: 2 }`
7. OR > 2.0 gets "Substantially elevated" context
8. OR < 1.0 gets historical inversion note
9. Beta coefficient > 0.1 with p < 5e-8 -> Tier2
10. Single-SNP association gets caveat
11. Integration: full GWAS scoring with strand awareness

---

## Blueprint 4: ClinVar 0-Star & Conflicting Classifications

**Severity:** HIGH
**Problem:** 0-star ClinVar entries (no assertion criteria provided) are currently presented
at Tier2, which is too generous. Conflicting classifications are silently dropped at
`scorer/mod.rs:139`.

### New Types

```rust
// New variant for ResultCategory enum
ResultCategory::ClinVarConflicting

// Conflict structure on ClinVarAnnotation
pub struct ClinVarConflictDetail {
    pub classifications: Vec<(String, u8)>,  // (classification, star_count)
    pub summary: String,
}
```

### Control Flow (revised scorer)

```
1. Allele match check (HomozygousRef -> SKIP)
2. Zero-star gate: review_stars == 0 -> Tier3 with caveat OR skip
3. Conflicting gate: significance contains "conflicting" ->
   Tier3 + ResultCategory::ClinVarConflicting with structured conflict data
4. BA1 filter: max_population_af > 0.05 -> downgrade or skip
5. Normal scoring: pathogenic/likely_pathogenic + stars >= 2 -> Tier1
```

### Changes

| File | Change |
|------|--------|
| `scorer/mod.rs` | Reorder control flow: allele-match -> zero-star -> conflicting -> BA1 -> scoring |
| `models/mod.rs` | Add `ClinVarConflicting` to `ResultCategory` enum |
| `models/annotation.rs` | Add `ClinVarConflictDetail` struct |
| `report/html.rs` | Render conflicting classifications with structure |
| `report/markdown.rs` | Same |

### Test Cases (7)
1. 0-star pathogenic -> Tier3 with "no assertion criteria" caveat
2. 0-star benign -> SKIP entirely
3. 2-star pathogenic -> Tier1 (unchanged)
4. Conflicting interpretations -> Tier3 + ClinVarConflicting category
5. VUS (uncertain significance) -> Tier3 with appropriate language
6. 1-star pathogenic -> Tier2 (existing behavior preserved)
7. Integration: mixed bag of ClinVar quality levels produces correct tier distribution

---

## Blueprint 5: ClinVar Germline/Somatic Classification

**Severity:** HIGH
**Problem:** Since 2024, ClinVar separates germline/somatic/oncogenicity classifications.
Current code treats all as germline, which can produce misleading results for somatic-only
variants.

### New Types

```rust
pub enum ClinVarClassificationType {
    Germline,
    Somatic,
    Oncogenicity,
}
```

### Schema Change

```sql
-- Add to clinvar table
ALTER TABLE clinvar ADD COLUMN classification_type TEXT DEFAULT 'germline';
-- Values: 'germline', 'somatic', 'oncogenicity'
```

### Changes

| File | Change |
|------|--------|
| `data/schema/schema.sql` | Add `classification_type` column to clinvar table |
| `data/seed/build_seed_db.py` | Populate classification_type for seed entries |
| `db/clinvar.rs` | Read classification_type, add to `ClinVarResult` |
| `models/annotation.rs` | Add `ClinVarClassificationType` to `ClinVarAnnotation` |
| `scorer/mod.rs` | Germline-preference: when multiple entries exist for same rsID, prefer germline |
| `scorer/mod.rs` | Somatic-only entries -> informational Tier3 notice |
| `report/html.rs` | Show classification type badge |

### Test Cases (6)
1. Germline pathogenic -> normal scoring (unchanged)
2. Somatic pathogenic -> Tier3 informational notice
3. Both germline + somatic for same rsID -> germline takes precedence
4. Oncogenicity classification -> Tier3 informational
5. Missing classification_type (old DB) -> defaults to germline
6. Schema migration: old DB without column still works

---

## Blueprint 6: Assembly/Build Tracking

**Severity:** MEDIUM
**Problem:** No concept of genome assembly (GRCh37 vs GRCh38) in the pipeline. Consumer
arrays use GRCh37, modern databases use GRCh38. rsID-based matching masks the problem
but position-based matching would fail silently.

### New Types

```rust
pub enum GenomeAssembly {
    GRCh36,   // Old 23andMe chips
    GRCh37,   // All current consumer arrays (hg19)
    GRCh38,   // Modern WGS VCFs, newer databases (hg38)
    Unknown,
}

pub struct ParsedFile {
    pub variants: Vec<Variant>,
    pub assembly: GenomeAssembly,
}
```

### Assembly Detection Rules

| Format | Assembly |
|--------|----------|
| 23andMe (all versions) | GRCh37 (always) |
| AncestryDNA | GRCh37 (always) |
| VCF `##reference=GRCh38` | GRCh38 |
| VCF `##reference=hg19` | GRCh37 |
| VCF `##contig=<...,assembly=GRCh38>` | GRCh38 |
| VCF (no reference header) | Unknown |

### Schema Change

```sql
CREATE TABLE IF NOT EXISTS db_metadata (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
-- INSERT OR REPLACE INTO db_metadata VALUES ('assembly', 'GRCh38');
```

### Changes

| File | Change |
|------|--------|
| `models/assembly.rs` | **New**: `GenomeAssembly` enum with Display, FromStr |
| `parser/liftover.rs` | **New**: `LiftoverAdapter` trait (stub, no implementation) |
| `models/mod.rs` | Re-export GenomeAssembly |
| `parser/mod.rs` | Add `ParsedFile`, `parse_auto_with_metadata()` |
| `parser/twentythreeandme.rs` | `detect_assembly()` scanning `# build 37` header lines |
| `parser/ancestry.rs` | `detect_assembly()` |
| `parser/vcf.rs` | `detect_assembly()` scanning `##reference` and `##contig` headers |
| `db/mod.rs` | `query_db_assembly()` reading `db_metadata` table |
| `lib.rs` | `check_assembly_compatibility()`, populate Report fields |
| `models/report.rs` | Add `input_assembly`, `db_assembly`, `assembly_warnings` |
| `report/markdown.rs` | Assembly rows in summary, warning section |
| `report/html.rs` | Same |
| `main.rs` | Use `parse_auto_with_metadata`, print assembly to stderr |
| `data/schema/schema.sql` | Add `db_metadata` table, bump schema to v3 |

### What is Deferred
- **LiftOver**: No implementation. `LiftoverAdapter` trait is the placeholder.
- **Position-based matching guard**: Exists as `AssemblyMismatch` error type but nothing calls it.
- **GRCh36 detection**: Documented but not actively handled.

### Test Cases (17)
- `from_str_header` tests for all aliases (hg19, hg38, b37, b38, GRCh37, GRCh38)
- Assembly compatibility tests (same known, different known, unknown + known)
- Parser assembly detection tests (23andMe, AncestryDNA, VCF with/without headers)
- Pipeline assembly mismatch warning test
- Report rendering tests

---

## Blueprint 7: Universal Clinical Confirmation Language

**Severity:** HIGH
**Problem:** Only ACMG SF genes get "confirm clinically" language. The research report
emphasizes ALL DTC findings need confirmation language, with urgency tiering.

### New Types

```rust
pub enum ConfirmationUrgency {
    HighImpact,                      // ACMG SF v3.2 genes
    ClinicalConfirmationRecommended, // ClinVar P/LP, PGx actionable
    InformationalOnly,               // GWAS, SNPedia, low-magnitude
}
```

### Changes

| File | Change |
|------|--------|
| `models/mod.rs` | Add `ConfirmationUrgency` enum |
| `models/report.rs` | Add `confirmation_urgency: ConfirmationUrgency` to `ScoredResult` |
| `models/report.rs` | Add `dtc_context: String` to `Report` |
| `scorer/mod.rs` | Assign `ConfirmationUrgency` to every `ScoredResult` |
| `scorer/mod.rs` | Push `dtc_raw_data_caveat()` into all `ScoredResult.limitations` |
| `lib.rs` | Strengthen PGx disclaimer with FDA language |
| `report/html.rs` | "Understanding This Report" section, urgency badges |
| `report/markdown.rs` | Same |

### Universal DTC Caveat Text

> "This result is derived from direct-to-consumer (DTC) microarray genotyping data,
> which has not been validated in a clinical laboratory setting. DTC genotyping has
> known limitations including strand ambiguity, limited coverage, and potential
> genotyping errors. Any clinically relevant finding should be confirmed through
> clinical-grade testing before making medical decisions."

### PGx FDA Disclaimer Addition

> "Pharmacogenomic results from consumer genotyping arrays have NOT been reviewed
> or approved by the FDA for clinical use. Do not alter any medication regimen
> based solely on these results. Consult a healthcare provider or clinical
> pharmacogenomics service for validated testing."

### Test Cases (12)
1. BRCA1 pathogenic -> `HighImpact`
2. ClinVar P/LP non-ACMG gene -> `ClinicalConfirmationRecommended`
3. PGx actionable finding -> `ClinicalConfirmationRecommended`
4. GWAS association -> `InformationalOnly`
5. SNPedia magnitude < 3 -> `InformationalOnly`
6. SNPedia magnitude >= 3 -> `ClinicalConfirmationRecommended`
7. Every ScoredResult has DTC caveat in limitations
8. PGx disclaimer contains FDA language
9. HTML "Understanding This Report" section renders
10. Markdown "Understanding This Report" section renders
11. HighImpact badge is visually distinct (red/urgent)
12. InformationalOnly has no scary language

---

## Blueprint 8: VCF Normalization

**Severity:** MEDIUM
**Problem:** VCF parser does not handle multi-allelic records, allele trimming, or
left-alignment. Multi-allelic GT index `"2"` falls through to `NoCall`.

### New Module: `normalizer/mod.rs`

```rust
pub enum NormalizationStatus {
    Original,
    Trimmed { bases_trimmed: u8 },
    MultiAllelicSplit { alt_index: usize },
    LeftAligned,
    LeftAlignedAndTrimmed,
    NormalizationFailed(String),
}

pub struct NormalizationStats {
    pub total_records: usize,
    pub multiallelic_split: usize,
    pub trimmed: usize,
    pub left_aligned: usize,
    pub failed: usize,
    pub issues: Vec<String>,
}
```

### Key Functions

- `normalize_vcf_record(chrom, pos, ref, alt_alleles, gt, rsid)` -> `Vec<(Variant, NormalizationStatus)>`
- `split_multiallelic(record)` -> `Vec<record>` (one per ALT allele)
- `trim_alleles(ref, alt, pos)` -> `(trimmed_ref, trimmed_alt, new_pos)`
- `left_align_approx(ref, alt, pos)` -> `(ref, alt, pos)` (string-level only, no FASTA)
- `resolve_gt_for_split(gt, alt_index, total_alts)` -> `Genotype`

### Multi-Allelic GT Resolution

For GT `"1/2"` with `ALT=G,T`:
- Split record 0 (ALT=G): GT index 1 = this ALT, index 2 = other -> `Heterozygous` for this record
- Split record 1 (ALT=T): GT index 2 = this ALT, index 1 = other -> `Heterozygous` for this record

### Trimming Algorithm

1. **Right-trim** (suffix, no pos change): strip matching trailing bases, min 1 base each
2. **Left-trim** (prefix, increment pos): strip matching leading bases, min 1 base each

### Changes

| File | Change |
|------|--------|
| `normalizer/mod.rs` | **New**: normalization module |
| `models/variant.rs` | Add `normalization_status: NormalizationStatus` |
| `models/report.rs` | Add `normalized_count: usize`, `normalization_issues: Vec<String>` |
| `parser/vcf.rs` | Rewrite to two-stage: raw parse + normalize. Delete `parse_gt` |
| `parser/mod.rs` | Add `pub mod normalize` |
| `lib.rs` | Add `pub mod normalizer`, zero-init Report normalization fields |
| `main.rs` | VCF format detection -> `parse_with_stats`, inject stats into report |
| `report/markdown.rs` | Normalization stats row, notes section |
| `report/html.rs` | Same |
| All test Variant literals | Add `normalization_status: NormalizationStatus::Original` |

### Priority Order
1. `split_multiallelic` + GT resolution (fixes actively wrong behavior)
2. `trim_alleles` (enables canonical form for position-based lookup)
3. `left_align_approx` (nice-to-have, stub with limitation doc)

### Test Cases (16)
- SNP unchanged, multi-allelic split (biallelic, triallelic)
- GT remapping for `0/1`, `1/1`, `1/2` with multiple ALTs
- Trimming: shared prefix, suffix, both, deletion, insertion
- Symbolic alleles (`<DEL>`) -> `NormalizationFailed`
- Stats accumulation
- Integration with VCF parser

---

## Blueprint 9: PGx CPIC Data Infrastructure

**Severity:** HIGH
**Problem:** Only 4 genes (CYP2C19, CYP2D6, CYP2C9, SLCO1B1) with 7 hardcoded allele
definitions. CPIC covers 24+ genes.

### New Gene Coverage

**Tier 1 (CPIC Level A, implement now):**
CYP2C19, CYP2D6, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, SLCO1B1, UGT1A1, VKORC1, CYP4F2, IFNL3

**Tier 2 (CPIC Level B, implement later):**
CYP2B6, NAT2, G6PD, HLA-A, HLA-B, RYR1, CACNA1S

### New Data Pipeline

```
CPIC REST API (api.cpicpgx.org/v1/)
    |  fetch_cpic.sh (curl per gene)
    v
data/raw/cpic/<GENE>.json
    |  import_cpic.py (parse JSON, batch insert)
    v
genesight.db (pgx_allele_definitions, pgx_diplotype_phenotypes,
              pgx_drug_recommendations, pgx_gene_metadata)
```

### Schema Additions

```sql
CREATE TABLE IF NOT EXISTS pgx_gene_metadata (
    gene TEXT PRIMARY KEY,
    cpic_level TEXT NOT NULL,          -- "A", "B", "C", "D"
    tier INTEGER NOT NULL,
    phenotype_classification_method TEXT,
    activity_score_ranges TEXT,        -- JSON
    has_cnv_risk INTEGER DEFAULT 0,
    microarray_coverage_note TEXT,
    cpic_guideline_url TEXT,
    last_updated TEXT
);

CREATE TABLE IF NOT EXISTS pgx_import_log (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene TEXT NOT NULL,
    cpic_version TEXT,
    import_source TEXT NOT NULL,
    allele_count INTEGER,
    diplotype_count INTEGER,
    recommendation_count INTEGER,
    imported_at TEXT DEFAULT (datetime('now'))
);
```

### Gene-Specific Activity Score Classifiers

Added to `phenotype.rs::classify_by_activity_score`:

| Gene | Poor | Intermediate | Normal/Extensive |
|------|------|-------------|-----------------|
| CYP3A5 | score <= 0.0 | 0.0 < score < 2.0 | score >= 2.0 |
| DPYD | score <= 0.0 (DPD Deficient) | 0.0 < score < 2.0 (HAS) | score >= 2.0 (Normal DPD) |
| TPMT/NUDT15 | score <= 0.0 | 0.0 < score < 2.0 | score >= 2.0 |
| VKORC1 | score < 1.5 (Highly Increased) | 1.5 <= score < 2.0 (Increased) | score >= 2.0 (Normal) |

### Changes

| File | Change |
|------|--------|
| `data/fetch/fetch_cpic.sh` | **New**: CPIC API download script |
| `data/import/import_cpic.py` | **New**: JSON -> SQLite importer |
| `data/schema/schema.sql` | Add `pgx_gene_metadata`, `pgx_import_log` tables |
| `data/seed/build_seed_db.py` | Expand with 5 new genes, seed `pgx_gene_metadata` |
| `pgx/phenotype.rs` | Gene-specific activity score classifiers |
| `pgx/definitions.rs` | `GeneMetadata` struct, `load_gene_metadata()` |
| `pgx/mod.rs` | Load gene metadata in `StarAlleleCaller::load()` |

### New VARIANTS in Seed DB

| rsID | Gene / Allele |
|------|---------------|
| rs776746 | CYP3A5*3 |
| rs3918290 | DPYD*2A |
| rs55886062 | DPYD*13 |
| rs67376798 | DPYD c.2846A>T |
| rs75017182 | DPYD HapB3 |
| rs1800462 | TPMT*2 |
| rs1800460 | TPMT*3B/*3A |
| rs1142345 | TPMT*3C/*3A |
| rs116855232 | NUDT15*3 |

### Test Cases (15+)
- Each new gene: score=0.0, score=1.0, score=2.0 boundary tests
- DPYD produces "DPD Activity" labels, not "Metabolizer"
- VKORC1 produces "Warfarin Sensitivity" labels
- TPMT includes phasing caveat
- Version tracking works with empty/populated pgx_data_version table

---

## Blueprint 10: PGx Phasing & Ambiguity

**Severity:** HIGH
**Problem:** Unphased array data cannot distinguish cis vs trans configurations for
multi-SNP star alleles. TPMT\*3A (rs1800460 + rs1142345 in cis) is the canonical example:
- Cis (\*3A/\*1): Intermediate Metabolizer
- Trans (\*3B/\*3C): Poor Metabolizer
These produce DIFFERENT clinical recommendations for thiopurine dosing.

### New Types

```rust
pub struct PhaseAmbiguity {
    pub ambiguous_positions: Vec<String>,
    pub alternative_diplotypes: Vec<AlternativeDiplotype>,
    pub clinical_impact: ClinicalPhaseImpact,
}

pub enum ClinicalPhaseImpact {
    Uniform,             // All alternatives -> same phenotype
    DifferentPhenotypes, // Alternatives span different metabolizer categories
}

pub struct AlternativeDiplotype {
    pub allele1: String,
    pub allele2: String,
    pub activity_score: f64,
    pub phenotype: String,
    pub population_prior: f64,
    pub is_primary: bool,
}

pub struct DiplotypePair {
    pub allele1: String,
    pub allele2: String,
    pub activity_score: f64,
}
```

### Algorithm: Consistency-Based Diplotype Enumeration

```
1. Build observed_counts: HashMap<rsid, u8> from user AlleleMatch values
   (HomozygousRef=0, Heterozygous=1, HomozygousAlt=2)

2. For each star allele, build alt_copies_required: HashMap<rsid, u8>
   (*1 = all zeros, *3A = {rs1800460: 1, rs1142345: 1})

3. Enumerate all pairs (A, B) where A <= B lexicographically:
   Check: for every defining position, reqs_A[pos] + reqs_B[pos] == observed[pos]

4. Collect all consistent pairs, sort by activity score ascending (most conservative first)

5. Primary = pair with lowest activity score (worst case = patient safety)
   Alternative = all other consistent pairs
```

### TPMT\*3A Walkthrough

Input: user heterozygous at rs1800460 AND rs1142345.

| Pair | Required at rs1800460 | Required at rs1142345 | Sum matches observed? | Consistent? |
|------|-----------------------|-----------------------|----------------------|-------------|
| \*3A/\*1 | 1+0=1 | 1+0=1 | Yes | Yes |
| \*3B/\*3C | 1+0=1 | 0+1=1 | Yes | Yes |
| \*3A/\*3A | 1+1=2 | 1+1=2 | No (observed=1) | No |

Result: Two consistent pairs. \*3B/\*3C has lower activity (0.0) vs \*3A/\*1 (1.0).
Conservative call = \*3B/\*3C (Poor Metabolizer). Ambiguity flagged.

### Changes

| File | Change |
|------|--------|
| `pgx/phasing.rs` | **New**: `detect_phase_ambiguity()`, consistency enumeration |
| `pgx/diplotype.rs` | Add `DiplotypePair`, `is_ambiguous`, `alternative_diplotypes` to `DiplotypeCall` |
| `pgx/phenotype.rs` | Add `AlternativePhenotype`, propagate ambiguity, set `clinical_confirmation_recommended` |
| `pgx/mod.rs` | Add `pub mod phasing` |
| `report/html.rs` | Amber "AMBIGUOUS CALL" badge, alternatives table, CSS |
| `report/markdown.rs` | `> **AMBIGUOUS CALL**` blockquote with alternatives |

### Test Cases (8)
1. TPMT het/het -> ambiguous, two consistent pairs, conservative = Poor
2. TPMT het at one site only -> unambiguous (\*3B/\*1 or \*3C/\*1)
3. TPMT hom-alt at one site -> unambiguous (\*3C/\*3C)
4. Single-SNP allele gene (CYP2C19) -> never ambiguous
5. All ref -> \*1/\*1, not ambiguous
6. Ambiguous with same phenotype -> `Uniform`, no clinical escalation
7. Ambiguous with different phenotypes -> `DifferentPhenotypes`, confirmation recommended
8. Report renders alternatives table with conservative label

---

## Total Impact Summary

| Metric | Count |
|--------|-------|
| New files | ~8 (assembly.rs, liftover.rs, normalizer/mod.rs, phasing.rs, import_cpic.py, fetch_cpic.sh, pgx_gene_metadata, docs) |
| Modified files | ~25 |
| New test cases | ~105 |
| Schema version | 2 -> 3 |
| PGx gene coverage | 4 -> 9+ genes |
| New enums/types | ~15 |

### Risk Assessment

| Change | Risk | Mitigation |
|--------|------|------------|
| Palindromic SNP fix | Low | Simple conditional, well-tested |
| ClinVar quality gates | Medium | May reduce reported findings; users expect results |
| PGx coverage expansion | Low | Additive data, existing engine handles it |
| PGx phasing | Medium | Complex algorithm, needs thorough testing |
| VCF normalization | Low | Only affects VCF users, existing formats unchanged |
| Assembly tracking | Low | Informational only in Phase 1 |
| Clinical confirmation | Low | Additive language, no behavior change |
