# M1: No Historical GWAS OR Inversion Detection

**Severity:** MEDIUM
**Status:** Not implemented

## Problem

The GWAS Catalog changed its curation policy around January 2021. Before that date,
OR < 1.0 entries were sometimes inverted (OR flipped to > 1.0 and the reported allele
swapped) so that all stored ORs would be > 1. After 2021, ORs are stored as-is.

The scorer at `scorer/mod.rs:736` handles OR < 1.0 with "Protective allele: odds ratio
below 1.0 suggests reduced risk" — which is correct IF the OR is reported accurately.
For pre-2021 entries that were inverted, the OR > 1.0 value is correct but the "risk allele"
may have been swapped, meaning the user might see the wrong allele flagged.

## Scientific Requirement

From the research report (Section: GWAS-Risikointerpretation):

> Für vor Jan 2021 kuratierte Studien wurde OR<1 teils invertiert und das berichtete
> Allel entsprechend gedreht, damit die gespeicherten ORs >1 sind.

> OR aus Top-Hits ohne Beachtung der historischen Inversion/Allelswap interpretieren.

## Fix Requirements

1. Add `study_date` or `catalog_date` column to GWAS schema
2. For entries with pre-2021 dates, add caveat: "Effect allele assignment may reflect
   historical GWAS Catalog conventions. Risk allele directionality should be verified
   against the original publication."
3. Alternatively, use harmonized GWAS summary statistics (GWAS-SSF) which have
   consistent effect_allele/other_allele encoding
