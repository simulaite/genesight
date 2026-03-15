# C1: PGx Scorer Has No Allele Check

**Severity:** CRITICAL
**Status:** Producing false results in live test
**Proven by:** 6 of 8 PGx results for huAE4518 are false positives

## Problem

The `score_pharma()` function at `crates/genesight-core/src/scorer/mod.rs:469-510` scores
PharmGKB annotations **without checking whether the user actually carries the variant allele**.
Every rsID match in the `pharmacogenomics` table produces a scored result regardless of
the user's genotype.

This is the exact bug the research report identifies as the #1 consumer-genomics pitfall:
> "rsID-based Anzeige von 'Pathogenic/Risk' ohne Allelvergleich produziert garantierte
> False Positives"

## Evidence (Live Test)

```
rs12248560  CYP2C19  User=CC (ref)  Alt=T  → Reported: "Ultrarapid Metabolizer (*17)"  → WRONG
rs4986893   CYP2C19  User=GG (ref)  Alt=A  → Reported: "Poor Metabolizer (*3)"         → WRONG
rs4244285   CYP2C19  User=GG (ref)  Alt=A  → Reported: "Poor Metabolizer (*2)"         → WRONG
rs1799853   CYP2C9   User=CC (ref)  Alt=T  → Reported: "Intermediate Metabolizer (*2)" → WRONG
rs1057910   CYP2C9   User=AA (ref)  Alt=C  → Reported: "Poor Metabolizer (*3)"         → WRONG
rs9923231   VKORC1   User=CC (ref)  Alt=T  → Reported: "Increased Sensitivity"         → WRONG
```

In all 6 cases, the user is homozygous reference — they carry zero copies of the variant
allele. The correct result is "Normal Metabolizer" or no result at all.

## Affected Code

```rust
// scorer/mod.rs:469-510
fn score_pharma(
    av: &AnnotatedVariant,
    rsid: &str,
    genotype: &str,
    pharma: &PharmaAnnotation,
) -> Option<ScoredResult> {
    let level = pharma.evidence_level.trim();
    // ... tier assignment ...
    // NO allele check anywhere — jumps straight to building ScoredResult
    Some(ScoredResult { ... })
}
```

Compare with `score_clinvar()` which correctly calls `clinvar_allele_check()` at line 230
and returns `None` when the user has 0 copies (line 234).

## Root Cause

The ClinVar scoring path was upgraded with allele checking but the PharmGKB path was not.
The `annotator/mod.rs` queries PharmGKB by rsID (line 95-96), and `score_pharma` blindly
trusts that an rsID match means the user is affected.

## Scientific Requirement

From the research report (Section: Allele Matching):
> "Der Vergleich erfolgt auf (chr,pos,REF,ALT) in einer festgelegten Assembly. rsID ist
> nur Lookup-Hilfe."

The pseudocode requires:
1. `normalize_user_alleles_to_plus(user, refdb)` — check alleles match
2. `count_target_allele(alleles, target_allele)` — count 0/1/2 copies
3. Only proceed if count > 0

## Fix Requirements

1. Add allele check gate to `score_pharma()` analogous to `clinvar_allele_check()`
2. When user carries 0 copies of ALT allele → return `None` (skip)
3. Include copy count in the result details
4. Handle palindromic SNPs and complement matching consistently

## Impact of Not Fixing

Every user who has any PGx-relevant rsID in their file — which is essentially every user —
will receive false metabolizer phenotype classifications. For CYP2C19, this could mean a
user with normal metabolism is told they are an "Ultrarapid Metabolizer" or "Poor Metabolizer",
potentially leading to medication errors if they share this with a physician.

The FDA has issued warnings about exactly this class of error:
> "The FDA warns that many PGx test claims have not been reviewed by the FDA and may not
> be scientifically supported."
