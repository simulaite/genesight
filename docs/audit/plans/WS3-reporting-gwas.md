# Work Stream 3: Report Rendering & GWAS Scoring Fixes

**Issues:** #7 (H3), #8 (H4), #9 (H5), #11 (H7), #12 (H8), #13 (M1), #14 (M2)
**Can run in parallel with:** WS1, WS2
**Cross-stream dependency:** Uses `PGX_FDA_DISCLAIMER` from WS1 Step 4 (scorer/mod.rs)

---

## Dependency Order

```
Step 1: H7 — Rename PolygenicRiskScore → GwasAssociation      [no deps, model change]
Step 2: H5 + M1 — Add GWAS caveats (OR explanation, inversion) [after Step 1]
Step 3: H8 — Add assembly mismatch per-result limitations       [no deps]
Step 4: H3 + M2 + H4 — Render urgency, DTC context, FDA PGx   [after Steps 1-3]
```

---

## Step 1: Rename `PolygenicRiskScore` → `GwasAssociation` (H7, Issue #11)

### Model Change

**File:** `crates/genesight-core/src/models/report.rs`

At line 56, rename the enum variant and add serde alias for backward compat:
```rust
/// GWAS single-variant associations (genome-wide significant hits)
#[serde(alias = "PolygenicRiskScore")]
GwasAssociation,
```

### Mechanical Rename (all match arms)

These files reference `PolygenicRiskScore` and must be updated to `GwasAssociation`:

| File | Lines | What to change |
|------|-------|---------------|
| `scorer/mod.rs` | ~518, 599, 668 | `ResultCategory::PolygenicRiskScore` → `ResultCategory::GwasAssociation` |
| `scorer/mod.rs` | Display impl (~875) | `"Polygenic Risk Score"` → `"GWAS Association"` |
| `report/html.rs` | ~679, 693 | `PolygenicRiskScore` → `GwasAssociation`, label `"PRS"` → `"GWAS"` |
| `report/markdown.rs` | ~222 | `PolygenicRiskScore` → `GwasAssociation` (keep sort_key = 3) |
| `tui/app.rs` (CLI) | ~286 | match arm rename |
| `tui/ui.rs` (CLI) | ~683, 1425 | match arm rename, label `"PRS"` → `"GWAS"` |
| `gui/state.rs` | ~201 | match arm rename |
| `gui/theme.rs` | ~197, 212 | match arm rename (keep colors) |

**Note:** The Rust compiler's exhaustive matching will catch any missed arms — the code won't compile until all are updated.

### Tests
```rust
#[test]
fn gwas_association_json_backward_compat() {
    let old = r#""PolygenicRiskScore""#;
    let cat: ResultCategory = serde_json::from_str(old).expect("deser");
    assert_eq!(cat, ResultCategory::GwasAssociation);
}

#[test]
fn gwas_association_serializes_new_name() {
    let json = serde_json::to_string(&ResultCategory::GwasAssociation).unwrap();
    assert!(json.contains("GwasAssociation"));
}
```

---

## Step 2: Add GWAS Caveats (H5 + M1, Issues #9, #13)

**File:** `crates/genesight-core/src/scorer/mod.rs`

### Add Constants

After existing caveat constants (around line 47):

```rust
/// Caveat: OR is relative, not absolute risk.
const GWAS_OR_RELATIVE_CAVEAT: &str = "The odds ratio (OR) is a relative measure \
    comparing your odds to someone without this allele. It is NOT an absolute \
    probability. The actual risk depends on baseline prevalence in your population. \
    An OR of 1.3 means approximately 1.3x the odds, not a 30% chance.";

/// Caveat: Beta coefficient context.
const GWAS_BETA_CAVEAT: &str = "The beta coefficient is the estimated effect size \
    per copy of the effect allele in the study's units. Without population mean and \
    standard deviation, clinical significance cannot be determined from DTC data alone.";

/// Caveat: Historical GWAS OR inversion.
const GWAS_OR_INVERSION_CAVEAT: &str = "The GWAS Catalog changed curation conventions \
    around January 2021. Earlier entries may have inverted odds ratios. Risk allele \
    directionality should be verified against the original publication.";
```

### Wire into Scoring

In `score_gwas_hit()` (around line 583), after the "Population-specific" limitation:
```rust
// Add OR/Beta explanation
if hit.odds_ratio.is_some() {
    limitations.push(GWAS_OR_RELATIVE_CAVEAT.to_string());
}
if hit.beta.is_some() && hit.odds_ratio.is_none() {
    limitations.push(GWAS_BETA_CAVEAT.to_string());
}
// Add historical inversion caveat (unconditional — no study_date field available)
limitations.push(GWAS_OR_INVERSION_CAVEAT.to_string());
```

Apply same to `score_gwas_hit_fallback()` (around line 698).

### Tests
- `gwas_result_has_or_caveat`: hit with odds_ratio → limitations contains "relative measure"
- `gwas_result_has_beta_caveat`: hit with beta only → limitations contains "beta coefficient"
- `gwas_result_has_inversion_caveat`: any hit → limitations contains "January 2021"

---

## Step 3: Assembly Mismatch Per-Result Limitations (H8, Issue #12)

**File:** `crates/genesight-core/src/scorer/mod.rs`

### Add New Function

```rust
/// Score with assembly context. If mismatched, append warning to every result.
pub fn score_variants_with_assembly(
    annotated: &[AnnotatedVariant],
    input_assembly: GenomeAssembly,
    db_assembly: GenomeAssembly,
) -> Vec<ScoredResult> {
    let mut results = score_variants(annotated);

    if !input_assembly.is_compatible_with(db_assembly)
        && input_assembly != GenomeAssembly::Unknown
        && db_assembly != GenomeAssembly::Unknown
    {
        let warning = format!(
            "Assembly mismatch: input uses {input_assembly} but database uses {db_assembly}. \
             Position-based lookups may be incorrect for this variant."
        );
        for result in &mut results {
            result.limitations.push(warning.clone());
        }
    }

    results
}
```

Add import: `use crate::models::GenomeAssembly;`

**File:** `crates/genesight-core/src/lib.rs`

Change line 213:
```rust
// OLD:
let scored = scorer::score_variants(&annotated);
// NEW:
let scored = scorer::score_variants_with_assembly(&annotated, input_assembly, db_assembly);
```

**File:** `data/seed/build_seed_db.py`

Add `db_metadata` table with assembly='GRCh38' (see WS2 Step 2d — whoever gets there first).

### Tests
- `assembly_mismatch_adds_limitation`: GRCh37 vs GRCh38 → every result has "Assembly mismatch"
- `assembly_match_no_limitation`: same assembly → no extra limitation
- `assembly_unknown_no_limitation`: Unknown assembly → no extra limitation

---

## Step 4: Render Urgency, DTC Context, FDA PGx (H3 + M2 + H4, Issues #7, #14, #8)

### 4a. HTML Renderer

**File:** `crates/genesight-core/src/report/html.rs`

**Add CSS** (in `write_html_head`, before `</style>`):
```css
.urgency-high {
    background: #fef2f2; border: 2px solid #dc2626; border-radius: 8px;
    padding: 0.75rem 1rem; margin: 0.5rem 0; color: #991b1b; font-weight: 600;
}
.urgency-clinical {
    background: #fff7ed; border: 2px solid #f59e0b; border-radius: 8px;
    padding: 0.5rem 1rem; margin: 0.5rem 0; color: #92400e;
}
.urgency-info { color: #6b7280; font-size: 0.8rem; font-style: italic; margin: 0.25rem 0; }
.fda-disclaimer {
    background: #fef2f2; border-left: 4px solid #dc2626; padding: 0.75rem 1rem;
    margin: 1rem 0; border-radius: 0 6px 6px 0; font-size: 0.85rem; color: #991b1b;
}
.dtc-context {
    background: #eff6ff; border-left: 4px solid #3b82f6; padding: 1rem 1.25rem;
    margin: 1.5rem 0; border-radius: 0 6px 6px 0; font-size: 0.9rem;
}
```

**Add import:** `use crate::models::report::ConfirmationUrgency;`

**Add helper:**
```rust
fn urgency_html(urgency: ConfirmationUrgency) -> &'static str {
    match urgency {
        ConfirmationUrgency::HighImpact =>
            "<div class=\"urgency-high\">&#9888; CLINICAL CONFIRMATION STRONGLY RECOMMENDED</div>",
        ConfirmationUrgency::ClinicalConfirmationRecommended =>
            "<div class=\"urgency-clinical\">Clinical confirmation recommended before medical decisions.</div>",
        ConfirmationUrgency::InformationalOnly =>
            "<div class=\"urgency-info\">Informational only — no clinical action warranted from DTC data alone.</div>",
    }
}
```

**Render urgency** in `write_results()` after summary div (around line 742):
```rust
out.push_str(urgency_html(result.confirmation_urgency));
```

Also render in Tier 1 detail section (around line 462).

**Add `write_dtc_context()`:**
```rust
fn write_dtc_context(out: &mut String, dtc_context: &str) {
    if dtc_context.is_empty() { return; }
    out.push_str("<div class=\"dtc-context\"><strong>About Your Data Source</strong>");
    out.push_str(&html_escape(dtc_context));
    out.push_str("</div>\n");
}
```

**Add `write_fda_pgx_disclaimer()`:**
```rust
fn write_fda_pgx_disclaimer(out: &mut String, results: &[ScoredResult]) {
    if !results.iter().any(|r| r.category == ResultCategory::Pharmacogenomics) { return; }
    out.push_str("<div class=\"fda-disclaimer\">");
    out.push_str("<strong>FDA Notice — Pharmacogenomic Results</strong>");
    out.push_str("<p>Pharmacogenomic results from consumer genotyping arrays have NOT been \
        reviewed or approved by the FDA for clinical use. Do not alter any medication regimen \
        based solely on these results.</p>");
    out.push_str("</div>\n");
}
```

**Call in `render()`** between disclaimer and results:
```rust
write_disclaimer(&mut out, &report.disclaimer);
write_dtc_context(&mut out, &report.dtc_context);
write_fda_pgx_disclaimer(&mut out, &report.results);
// ... write_summary, write_results ...
```

### 4b. Markdown Renderer

**File:** `crates/genesight-core/src/report/markdown.rs`

**Add import:** `use crate::models::report::ConfirmationUrgency;`

**Render urgency** in `write_results()` after details line:
```rust
match result.confirmation_urgency {
    ConfirmationUrgency::HighImpact => {
        writeln!(out, "> **!! CLINICAL CONFIRMATION STRONGLY RECOMMENDED !!**").ok();
    }
    ConfirmationUrgency::ClinicalConfirmationRecommended => {
        writeln!(out, "> *Clinical confirmation recommended.*").ok();
    }
    ConfirmationUrgency::InformationalOnly => {}
}
```

**Add `write_dtc_context()`:**
```rust
fn write_dtc_context(out: &mut String, dtc_context: &str) {
    if dtc_context.is_empty() { return; }
    writeln!(out, "### About Your Data Source\n").ok();
    writeln!(out, "> {dtc_context}\n").ok();
}
```

**Add `write_fda_pgx_disclaimer()`:**
```rust
fn write_fda_pgx_disclaimer(out: &mut String, results: &[ScoredResult]) {
    if !results.iter().any(|r| r.category == ResultCategory::Pharmacogenomics) { return; }
    writeln!(out, "### FDA Notice — Pharmacogenomic Results\n").ok();
    writeln!(out, "> **Pharmacogenomic results from consumer genotyping arrays have NOT been \
        reviewed or approved by the FDA for clinical use.** Do not alter any medication \
        based solely on these results.\n").ok();
}
```

Call in `render()` between disclaimer and summary.

### Tests

**HTML:**
- `render_urgency_high_impact`: HighImpact → HTML contains "urgency-high" + "CLINICAL CONFIRMATION"
- `render_urgency_clinical`: → HTML contains "urgency-clinical"
- `render_dtc_context`: non-empty dtc_context → HTML contains "About Your Data Source"
- `render_fda_when_pgx`: PGx result → HTML contains "fda-disclaimer"
- `render_no_fda_without_pgx`: no PGx → no "fda-disclaimer"

**Markdown:**
- `render_urgency_high_impact_md`: HighImpact → contains "!! CLINICAL CONFIRMATION"
- `render_dtc_context_md`: → contains "About Your Data Source"
- `render_fda_when_pgx_md`: → contains "FDA Notice"

---

## Verification

```bash
cargo test -p genesight-core
cargo test  # full workspace including CLI and GUI
cargo clippy
cargo fmt -- --check
```

Generate test reports and visually inspect:
```bash
cargo run -p genesight-cli -- analyze tests/fixtures/pgp/huAE4518_23andme_v4.txt --format html > /tmp/report.html
open /tmp/report.html
```

Check for:
- Red urgency banner on BRCA1 findings
- Orange urgency on ClinVar P/LP and PGx
- FDA disclaimer section (only if PGx results present)
- "About Your Data Source" section
- GWAS results labeled "GWAS Association" (not "Polygenic Risk Score")
- OR caveat in GWAS limitations
