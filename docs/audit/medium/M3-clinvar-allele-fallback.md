# M3: ClinVar Allele Check Falls Back to rsID-Only When Variants Table Missing

**Severity:** MEDIUM
**Status:** Silent degradation path exists

## Problem

The `clinvar_allele_check()` function at `scorer/mod.rs:118-184` correctly checks whether
the user carries the ClinVar ALT allele. However, when the allele data is unavailable
(because the `variants` table doesn't exist or has no entry for the rsID), the function
returns `ClinvarAlleleResult::NoAlleleData`.

At `scorer/mod.rs:255`, this result causes the scorer to proceed WITHOUT any allele
verification:

```rust
ClinvarAlleleResult::NoAlleleData => None,  // copies = None, proceed without check
```

When `copies` is `None`, the scorer continues to produce a `ScoredResult` without the
allele copy count — effectively reverting to pure rsID matching.

The `variants` table is populated by `import_clinvar.py` from ClinVar's VCF-style
ref/alt columns. If a database was built before this table was added, or if the import
fails to populate it, all ClinVar results silently degrade to rsID-only matching.

## Impact

This is a safety regression path: the allele check exists but can be silently bypassed.
A database missing the `variants` table (or with incomplete data) will produce false
positives for ClinVar pathogenic variants where the user is homozygous reference.

## Fix Requirements

1. When `NoAlleleData` is returned, add a limitation to the result: "Allele verification
   was not possible for this variant. Result is based on rsID matching only and may be
   a false positive."
2. Consider downgrading tier when allele verification is not possible
3. Log a warning when the `variants` table is missing or empty
4. Add a database integrity check that validates the `variants` table exists and has
   adequate coverage
