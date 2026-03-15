# C2: New PGx Pipeline (Diplotype/Phasing/Coverage) Is Dead Code

**Severity:** CRITICAL
**Status:** Correct implementation exists but is unreachable from the live pipeline

## Problem

The codebase contains two completely separate PGx code paths:

**Path A (LIVE — wired into `lib.rs::analyze()`):**
```
lib.rs:208-209 → annotator::annotate_variants_with_config()
    → db::pharmgkb::query_batch()  — raw per-rsID PharmGKB records
    → scorer::score_variants() → score_pharma()  — no allele check, no diplotype
```

**Path B (DEAD — never called from `analyze()`):**
```
pgx/definitions.rs::load_allele_definitions()
    → pgx/diplotype.rs::call_diplotype()     — coverage-aware, enumerated
    → pgx/phasing.rs::detect_phase_ambiguity() — TPMT*3A cis/trans detection
    → pgx/phenotype.rs::call_phenotype_with_coverage() — Indeterminate for missing data
```

Path B implements all research requirements correctly:
- CoverageStatus tracking (Complete/Partial/Insufficient)
- Diplotype enumeration from unphased data
- TPMT*3A phasing ambiguity detection
- "Indeterminate" phenotype when data is insufficient
- Gene-specific activity score classifiers for 9 genes

Path A has none of this. Path A is what runs.

## Evidence

Grep for any call to the new pipeline functions outside their own module:
- `call_diplotype` → 0 matches outside `pgx/diplotype.rs`
- `detect_phase_ambiguity` → 0 matches outside `pgx/phasing.rs`
- `call_phenotype_with_coverage` → 0 matches outside `pgx/phenotype.rs`
- `load_allele_definitions` → 0 matches outside `pgx/definitions.rs`

`lib.rs::analyze_with_config_and_assembly()` (lines 189-254) calls:
1. `annotator::annotate_variants_with_config()` — line 208-209
2. `scorer::score_variants()` — line 213

Neither of these calls any function from `pgx/diplotype.rs`, `pgx/phasing.rs`, or
`pgx/definitions.rs`.

## Modules Affected

| Module | Status | What It Does Right |
|--------|--------|--------------------|
| `pgx/diplotype.rs` | Dead code | CoverageStatus, diplotype enumeration, missing-data != *1 |
| `pgx/phasing.rs` | Dead code | TPMT*3A ambiguity, conservative call (worst-case phenotype) |
| `pgx/phenotype.rs` | Partially used | `call_phenotype()` called by `mod.rs`, but `call_phenotype_with_coverage()` is dead |
| `pgx/definitions.rs` | Dead code | Structured allele definition loader (but has schema mismatch — see C4) |
| `pgx/mod.rs` | Partially used | `StarAlleleCaller::call_gene()` exists but is also not called from `analyze()` |

## Scientific Requirements Not Met (Due to Dead Code)

From the research report:

1. **Missing data treated as wildtype** (Section: PGx Star-Allele Calling):
   > "'Not tested' is not 'wildtype'. When defining sites are missing, a 'Normal'
   > result is often not derivable and must be returned as 'Unclear/Indeterminate'."

2. **Phase ambiguity** (Section: PGx Star-Allele Calling):
   > "Forcing TPMT*3A without phasing → incorrect diplotype (cis/trans)."

3. **Diplotype inference** (Section: PGx Star-Allele Calling):
   > "Enumerate the possibility space of all haplotype pairs that explain the (unphased)
   > genotype counts; in case of ambiguity, return 'ambiguous call' rather than a forced result."

All three are correctly implemented in Path B but unreachable.

## Fix Requirements

1. Wire `StarAlleleCaller` (or the newer diplotype pipeline) into `lib.rs::analyze()`
2. Replace or supplement the PharmGKB annotation path with proper star-allele calling
3. Ensure the new pipeline's results flow into `ScoredResult` for report rendering
4. Fix C4 (schema mismatch) first if using `definitions.rs`

## Relationship to Other Issues

- **Blocks C4**: Schema mismatch must be fixed before `definitions.rs` can query real DB
- **Subsumes C1**: Wiring the new pipeline replaces `score_pharma()` with allele-aware logic
- **Enables Blueprint 10**: Phasing ambiguity detection depends on this being wired in
