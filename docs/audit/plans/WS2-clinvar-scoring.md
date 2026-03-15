# Work Stream 2: ClinVar Scoring & Import Fixes

**Issues:** #3 (C3), #5 (H1), #6 (H2), #15 (M3)
**Can run in parallel with:** WS1, WS3
**No cross-stream dependencies**

---

## Dependency Order

```
Step 1: C3 + H2 — Fix import_clinvar.py (star mapping + classification_type) [no deps]
Step 2: Seed DB — Fix build_seed_db.py and rebuild                            [after Step 1]
Step 3: H1 + M3 — Fix scorer (AR carrier detection + NoAlleleData caveat)     [no deps, parallel]
```

---

## Step 1: Fix `import_clinvar.py` (C3 + H2, Issues #3, #6)

**File:** `data/import/import_clinvar.py`

### 1a. Fix REVIEW_STARS (C3)

Replace lines 38-47 with correct ClinVar star mapping per https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/:

```python
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations of pathogenicity": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
}
```

Specific corrections:
- `"reviewed by expert panel"`: 4 → **3**
- `"criteria provided, multiple submitters, no conflicts"`: 3 → **2**
- `"criteria provided, conflicting classifications"`: 2 → **1**
- `"criteria provided, conflicting interpretations of pathogenicity"`: 2 → **1**
- `"criteria provided, single submitter"`: 2 → **1**
- `"no assertion criteria provided"`: 1 → **0**

### 1b. Add classification_type Parsing (H2)

Add column constants:
```python
COL_ORIGIN = 14
COL_SOMATIC_SIGNIFICANCE = 34  # SomaticClinicalImpact column
COL_ONCOGENICITY = 37          # Oncogenicity column
```

Add function:
```python
def determine_classification_type(fields):
    """Determine germline/somatic/oncogenicity from ClinVar variant_summary row."""
    origin = fields[COL_ORIGIN].strip().lower() if len(fields) > COL_ORIGIN else ""
    somatic_sig = fields[COL_SOMATIC_SIGNIFICANCE].strip() if len(fields) > COL_SOMATIC_SIGNIFICANCE else "-"
    onco_sig = fields[COL_ONCOGENICITY].strip() if len(fields) > COL_ONCOGENICITY else "-"

    origins = {o.strip() for o in origin.split(";")} if origin else set()
    has_germline = bool(origins & {"germline", "inherited", "de novo"})
    has_somatic = "somatic" in origins

    if has_somatic and not has_germline:
        return "somatic"
    if somatic_sig not in ("-", "", "na") and not has_germline:
        return "somatic"
    if onco_sig not in ("-", "", "na") and not has_germline:
        return "oncogenicity"
    return "germline"
```

Update `SCHEMA_SQL` to include `classification_type TEXT DEFAULT 'germline'` in clinvar table.

Update the INSERT statement to include classification_type.

Update the `best` dict tracking to include classification_type.

---

## Step 2: Fix Seed Database (Issues #3, #6)

**File:** `data/seed/build_seed_db.py`

### 2a. Add classification_type to clinvar table

In the clinvar CREATE TABLE (around line 634), add:
```sql
classification_type TEXT DEFAULT 'germline'
```

### 2b. Update CLINVAR_ENTRIES tuples

Add `"germline"` as new element to each entry tuple. All seed entries are germline.

### 2c. Fix star ratings in seed data

Check each entry against corrected mapping:
- rs1800553 (GJB2): if currently 3 → change to 2 ("multiple submitters" = 2 stars)
- Verify all others are correct

### 2d. Add db_metadata table with assembly

```python
cur.execute("""
    CREATE TABLE IF NOT EXISTS db_metadata (
        key TEXT PRIMARY KEY,
        value TEXT NOT NULL
    )
""")
cur.execute("INSERT INTO db_metadata (key, value) VALUES ('assembly', 'GRCh38')")
```

### 2e. Rebuild

```bash
python3 data/seed/build_seed_db.py
```

---

## Step 3: Fix Scorer — AR Carrier Detection + NoAlleleData Caveat (H1 + M3, Issues #5, #15)

**File:** `crates/genesight-core/src/scorer/mod.rs`

### 3a. Add AR_GENES Constant (H1)

Add after the ACMG_SF_GENES constant (around line 32):

```rust
/// Curated autosomal recessive (AR) genes. Heterozygous pathogenic = carrier, not affected.
/// Sources: OMIM, ClinGen Gene-Disease Validity.
/// Excludes genes with both AD and AR modes where het is clinically significant (e.g., BRCA1).
const AR_GENES: &[&str] = &[
    "CFTR",      // Cystic Fibrosis
    "HBB",       // Sickle Cell / Beta-thalassemia
    "HEXA",      // Tay-Sachs
    "GJB2",      // Hearing Loss (AR)
    "GJB6",      // Hearing Loss (AR)
    "SLC26A4",   // Pendred Syndrome
    "PAH",       // Phenylketonuria
    "SMN1",      // Spinal Muscular Atrophy
    "GBA1",      // Gaucher Disease
    "ASPA",      // Canavan Disease
    "MEFV",      // Familial Mediterranean Fever
    "ATP7B",     // Wilson Disease
    "HFE",       // Hemochromatosis
    "BLM",       // Bloom Syndrome
    "FANCA",     // Fanconi Anemia
    "FANCC",     // Fanconi Anemia
    "MYO7A",     // Usher Syndrome
    "USH2A",     // Usher Syndrome
    "GALT",      // Galactosemia
    "ACADM",     // MCAD Deficiency
    "SMPD1",     // Niemann-Pick
    "BCKDHA",    // Maple Syrup Urine Disease
    "BCKDHB",    // Maple Syrup Urine Disease
    "SERPINA1",  // Alpha-1 Antitrypsin Deficiency
    "CYP21A2",  // Congenital Adrenal Hyperplasia
];
```

### 3b. Add Carrier Detection Logic

In `score_clinvar()`, after allele check determines `copies` (around line 237), extract:

```rust
let alt_copies: Option<u8> = match allele_result {
    ClinvarAlleleResult::Copies(c) => Some(c),
    _ => None,
};

let is_ar_carrier = alt_copies == Some(1)
    && clinvar.gene_symbol.as_deref().is_some_and(|g| AR_GENES.contains(&g));
```

Then in the pathogenic scoring paths (lines 420-461), when `is_ar_carrier`:
- Use `ResultCategory::CarrierStatus` instead of `MonogenicDisease`
- Use `ConfirmationUrgency::ClinicalConfirmationRecommended` (not HighImpact)
- Add carrier-specific limitation text:
  "You carry one copy of a pathogenic variant in an autosomal recessive gene. Carriers
  typically do not develop symptoms but may pass the variant to offspring."
- Add "(carrier)" suffix to summary

Apply same logic to all three pathogenic paths: 0-star, 1-star, >=2-star.

### 3c. Add NoAlleleData Caveat (M3)

Change `score_clinvar()` line 255 from:
```rust
ClinvarAlleleResult::NoAlleleData => None,
```
to:
```rust
ClinvarAlleleResult::NoAlleleData => {
    tracing::warn!(rsid = rsid, "ClinVar allele verification not possible — no ref/alt data");
    Some("Allele verification was not possible for this variant (no reference/alternate \
         allele data in the database). This result is based on rsID matching only and \
         may include false positives.".to_string())
},
```

Add tier downgrade flag:
```rust
let no_allele_verification = matches!(allele_result, ClinvarAlleleResult::NoAlleleData);
```

In pathogenic scoring, when `no_allele_verification && base_tier == Tier1Reliable`:
downgrade to `Tier2Probable`.

### Tests

**Carrier detection (H1):**
- `cftr_het_pathogenic_is_carrier`: CFTR, 1 copy, 4-star → CarrierStatus, Tier1
- `cftr_hom_pathogenic_is_monogenic`: CFTR, 2 copies, 4-star → MonogenicDisease
- `brca1_het_is_still_monogenic`: BRCA1 not in AR_GENES → MonogenicDisease
- `hbb_het_pathogenic_is_carrier`: HBB, 1 copy → CarrierStatus

**NoAlleleData caveat (M3):**
- `no_allele_data_adds_limitation`: ref/alt=None → limitation contains "rsID matching"
- `no_allele_data_downgrades_tier1`: 3-star pathogenic, no allele data → Tier2 (not Tier1)
- `with_allele_data_no_rsid_limitation`: ref/alt present → no "rsID matching" text

---

## Verification

```bash
cargo test -p genesight-core
cargo clippy
cargo fmt -- --check
```

Verify with huAE4518: ClinVar results should have correct star assignments after seed rebuild.
