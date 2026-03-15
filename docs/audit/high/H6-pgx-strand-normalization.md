# H6: PGx Star Allele Caller Uses Raw Genotype Characters

**Severity:** HIGH
**Status:** Strand normalization bypassed in PGx path

## Problem

The `StarAlleleCaller::call_gene()` function at `pgx/mod.rs:150` counts alt allele
copies using a raw character comparison on the genotype string:

```rust
let alt_char = def.alt_allele.chars().next()?;
let count = observed.chars().filter(|&c| c == alt_char).count() as u8;
```

This operates on the un-normalized genotype string from the user's file. If the user's
data is reported on the opposite strand (complement), the character comparison will fail:
- Database says alt = `T`
- User's file (on complement strand) reports `A` (which IS the T allele on the other strand)
- `observed.chars().filter(|&c| c == 'T')` returns 0 → treated as wildtype

## Contrast with ClinVar/GWAS Path

The ClinVar scoring path (`scorer/mod.rs:118-184`) correctly implements:
1. Direct match check
2. Complement match check (lines 169-180)
3. Palindromic SNP detection with frequency-based resolution

The GWAS path (`allele/mod.rs:205-261`) also uses `match_single_allele()` with
complement checking.

Only the PGx path skips strand awareness entirely.

## Scientific Requirement

From the research report (Section: Allele Matching):

> For robust allele comparisons, the correct rule is: normalize alleles to the same build and
> the same (plus-)orientation, and only then compare.

The PGx pseudocode in the research explicitly requires:
> 1) Harmonize input: bring variants to one assembly, normalize alleles to plus strand,
>    verify REF against reference.
> 2) Load allele definitions
> 3) Named-Allele-Matching: for each star definition, check whether the user's genotype data
>    supports this haplotype definition

Step 1 (normalization) is missing from the PGx path.

## Fix Requirements

1. Use `match_alleles()` or `match_single_allele()` from the `allele` module instead
   of raw character comparison
2. Or: normalize genotype characters to plus-strand before entering PGx logic
3. Handle palindromic PGx sites (rare but possible) with appropriate ambiguity flags
4. Add test cases for complement-strand genotype inputs in PGx calling
