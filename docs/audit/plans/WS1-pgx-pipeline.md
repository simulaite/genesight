# Work Stream 1: PGx Pipeline Fixes

**Issues:** #1 (C1), #2 (C2), #4 (C4), #8 (H4), #10 (H6)
**Can run in parallel with:** WS2, WS3
**Cross-stream dependency:** FDA disclaimer constant (Step 5) is also used by WS3 renderers

---

## Dependency Order

```
Step 1: C4 — Fix definitions.rs schema mismatch     [no deps]
Step 2: C1 — Add allele check to score_pharma()      [no deps, parallel with Step 1]
Step 3: H6 — Fix PGx strand normalization             [after Step 1]
Step 4: H4 — Add FDA PGx disclaimer constant          [no deps]
Step 5: C2 — Wire PGx pipeline into analyze()         [after Steps 1, 3, 4]
```

---

## Step 1: Fix `definitions.rs` Schema Mismatch (C4, Issue #4)

**File:** `crates/genesight-core/src/pgx/definitions.rs`

### Problem
SQL query at line 63-66 uses wrong column names vs actual schema:

| Actual Column | definitions.rs (WRONG) |
|---------------|----------------------|
| `allele_name` | `star_allele` |
| `rsid` | `defining_rsid` |
| `alt_allele` | `variant_allele` |
| `activity_score` | `activity_value` |

### Changes

1. **Rename `AlleleDefiningVariant` struct fields** (lines 15-21):
   ```rust
   pub struct AlleleDefiningVariant {
       pub allele_name: String,    // was: star_allele
       pub rsid: String,           // was: defining_rsid
       pub alt_allele: String,     // was: variant_allele
       pub function: String,
       pub activity_score: f64,    // was: activity_value
   }
   ```

2. **Fix SQL query** (lines 63-66):
   ```rust
   "SELECT gene, allele_name, rsid, alt_allele, function, activity_score \
    FROM pgx_allele_definitions \
    ORDER BY gene, allele_name"
   ```

3. **Fix row mapping** (lines 74-80): Update field names in the `AlleleDefiningVariant` constructor.

4. **Fix `defining_rsids` push** (line 93): `variant.defining_rsid` → `variant.rsid`

5. **Fix `load_drug_recommendations()` SQL** (line 161): Remove `cpic_guideline_url` from SELECT (column doesn't exist in schema). Set `DrugRecommendation.cpic_guideline_url` to `None`.

6. **Fix test fixture** (lines 201-215): Change CREATE TABLE and INSERT to use correct column names.

**Then fix downstream consumers:**

7. **`pgx/diplotype.rs`**: All refs to `AlleleDefiningVariant` fields:
   - `def_var.defining_rsid` → `def_var.rsid`
   - `def_var.variant_allele` → `def_var.alt_allele`
   - Also fix test fixtures that construct `AlleleDefiningVariant`

8. **`pgx/phasing.rs`**: Same field renames:
   - `def_var.defining_rsid` → `def_var.rsid`
   - `def_var.variant_allele` → `def_var.alt_allele`
   - Fix test fixtures (lines 223-263)

### Tests
- All existing `definitions.rs`, `diplotype.rs`, `phasing.rs` tests must pass after fixture updates
- Run: `cargo test -p genesight-core`

---

## Step 2: Add Allele Check to `score_pharma()` (C1, Issue #1)

**File:** `crates/genesight-core/src/scorer/mod.rs`

### Problem
`score_pharma()` at lines 469-510 has NO genotype check. Reports phenotype for every rsID match. 6 of 8 PGx results for huAE4518 are false positives.

### Change

Add allele check at top of `score_pharma()` (line 474), reusing the existing `clinvar_allele_check()` function which already handles:
- Direct match
- Complement match
- Palindromic detection
- Indel detection
- NoCall handling

```rust
fn score_pharma(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    pharma: &PharmaAnnotation,
) -> Option<ScoredResult> {
    // NEW: Gate on allele match — skip if user doesn't carry the variant allele
    let allele_result = clinvar_allele_check(av);
    match allele_result {
        ClinvarAlleleResult::Copies(0) => return None,  // homozygous ref → skip
        ClinvarAlleleResult::IndelNotDetectable => return None,
        ClinvarAlleleResult::NoCallGenotype => return None,
        _ => {}  // Copies(1|2), Palindromic, NoAlleleData → proceed
    }

    // ... rest of existing function unchanged ...
}
```

### Tests
- `pharma_homozygous_ref_skipped`: genotype GG, ref=G, alt=A → returns None
- `pharma_heterozygous_passes`: genotype GA, ref=G, alt=A → returns Some(result)
- `pharma_no_allele_data_passes`: no ref/alt in DB → returns Some (legacy)
- Existing `analyze_full_pipeline` in lib.rs must still pass

---

## Step 3: Fix PGx Strand Normalization (H6, Issue #10)

**File:** `crates/genesight-core/src/pgx/mod.rs`

### Problem
Line 150 does raw char comparison: `observed.chars().filter(|&c| c == alt_char).count()`
No complement matching. Opposite-strand data → wrong counts.

### Change

Replace lines 149-150:
```rust
// OLD:
let alt_char = def.alt_allele.chars().next()?;
let count = observed.chars().filter(|&c| c == alt_char).count() as u8;

// NEW:
let alt_char = def.alt_allele.chars().next()?;
let count = observed.chars().filter(|&c| {
    let m = crate::allele::match_single_allele(alt_char, c);
    matches!(m, crate::allele::AlleleMatch::DirectMatch | crate::allele::AlleleMatch::ComplementMatch)
}).count() as u8;
```

Add import at top of file:
```rust
use crate::allele::{match_single_allele, AlleleMatch};
```

**Also fix `pgx/diplotype.rs`** — same raw char comparison pattern at lines 170-181.

### Tests
- `call_gene_complement_strand`: DB alt=A, user genotype=CT (complement of GA) → count=1
- Existing tests: `call_gene_homozygous_alt` (AA), `call_gene_heterozygous` (GA) still pass

---

## Step 4: Add FDA PGx Disclaimer (H4, Issue #8)

**File:** `crates/genesight-core/src/scorer/mod.rs`

### Change

Add constant after line 47 (after `ONCOGENICITY_LIMITATION`):
```rust
const PGX_FDA_DISCLAIMER: &str = "Pharmacogenomic results from consumer genotyping \
    arrays have NOT been reviewed or approved by the U.S. Food and Drug Administration \
    (FDA) for clinical use. Do not alter any medication regimen based solely on these \
    results. Consult a healthcare provider or clinical pharmacogenomics service for \
    validated testing.";
```

In `score_pharma()` at line 508, change:
```rust
// OLD:
limitations: vec![DTC_RAW_DATA_CAVEAT.to_string()],
// NEW:
limitations: vec![PGX_FDA_DISCLAIMER.to_string(), DTC_RAW_DATA_CAVEAT.to_string()],
```

Also make `PGX_FDA_DISCLAIMER` pub(crate) so lib.rs can use it in Step 5.

### Tests
- `pharma_result_has_fda_disclaimer`: verify limitations contains "FDA"

---

## Step 5: Wire PGx Pipeline into `analyze()` (C2, Issue #2)

**File:** `crates/genesight-core/src/lib.rs`

### Problem
`analyze_with_config_and_assembly()` (lines 189-254) never calls `StarAlleleCaller`,
`call_diplotype`, `detect_phase_ambiguity`, or `call_phenotype_with_coverage`.

### Changes

1. **Add imports:**
   ```rust
   use pgx::StarAlleleCaller;
   use std::collections::HashMap;
   ```

2. **Add `run_pgx_pipeline()` function** (~50 lines):
   ```rust
   fn run_pgx_pipeline(
       main_db: &Connection,
       annotated: &[models::AnnotatedVariant],
   ) -> Vec<models::report::ScoredResult> {
       // Try to create StarAlleleCaller; return empty if table missing
       let caller = match StarAlleleCaller::from_db(main_db) {
           Ok(c) if !c.supported_genes().is_empty() => c,
           _ => return Vec::new(),
       };

       // Build rsID → genotype string map
       let mut genotype_map: HashMap<String, String> = HashMap::new();
       for av in annotated {
           if let Some(rsid) = &av.variant.rsid {
               let gt_str = format!("{}", av.variant.genotype);
               if gt_str != "--" && gt_str != "NoCall" {
                   genotype_map.insert(rsid.clone(), gt_str);
               }
           }
       }

       let mut results = Vec::new();
       for gene in caller.supported_genes() {
           if let Some(call) = caller.call_gene(gene, &genotype_map) {
               results.push(pgx_call_to_scored(call));
           }
       }
       results
   }
   ```

3. **Add `pgx_call_to_scored()` converter** (~35 lines):
   Creates a `ScoredResult` from a `StarAlleleCall`. Uses a synthetic `AnnotatedVariant`
   with rsid = `"PGx:{gene}"` as a gene-level marker. Sets:
   - `tier: Tier1Reliable`
   - `category: Pharmacogenomics`
   - `confirmation_urgency: ClinicalConfirmationRecommended`
   - `limitations`: call.limitations + PGX_FDA_DISCLAIMER + DTC_RAW_DATA_CAVEAT

4. **Wire into `analyze_with_config_and_assembly()`** after line 213:
   ```rust
   // Step 2: Score annotated variants
   let scored = scorer::score_variants(&annotated);

   // Step 2b: PGx gene-level analysis (star allele calling)
   let pgx_results = run_pgx_pipeline(main_db, &annotated);

   // Merge
   let mut all_scored = scored;
   all_scored.extend(pgx_results);
   ```
   Then use `all_scored` in the tier filter (line 216).

### Tests
- `analyze_with_pgx_pipeline`: DB with pgx_allele_definitions → gene-level PGx results appear
- `analyze_pgx_table_missing_graceful`: DB without table → no crash, no PGx results
- `pgx_result_has_fda_disclaimer`: gene-level results include FDA text
- All existing `analyze_*` tests pass (they don't have pgx_allele_definitions table)

---

## Verification

After all 5 steps:
```bash
cargo test -p genesight-core
cargo test
cargo clippy
cargo fmt -- --check
cargo run -p genesight-cli -- analyze tests/fixtures/pgp/huAE4518_23andme_v4.txt --format json 2>&1 | head -100
```

Verify huAE4518 PGx results:
- rs12248560 CYP2C19 CC (ref) → should NOT appear (was false positive)
- rs3892097 CYP2D6 CT (het) → should still appear with correct phenotype
- New gene-level PGx results should appear (CYP2D6 *1/*4, etc.)
