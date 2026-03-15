# Live Test: PGP Individual huAE4518

**Test file:** `tests/fixtures/pgp/huAE4518_23andme_v4.txt`
**Date:** 2026-03-15
**Build:** All 255 tests pass, 0 failures

## Parsing

| Metric | Value |
|--------|-------|
| Format | 23andMe v4 |
| Total variants | 601,783 |
| Annotated | 30 |
| Scored | 17 |

## PGx Results (8 results, 6 false positives)

| rsID | Gene | User Genotype | Ref/Alt | Reported Phenotype | Correct? |
|------|------|---------------|---------|-------------------|----------|
| rs12248560 | CYP2C19 | CC | C/T | Ultrarapid Metabolizer (*17) | **FALSE** — CC is ref, should be Normal |
| rs4986893 | CYP2C19 | GG | G/A | Poor Metabolizer (*3) | **FALSE** — GG is ref |
| rs4244285 | CYP2C19 | GG | G/A | Poor Metabolizer (*2) | **FALSE** — GG is ref |
| rs1799853 | CYP2C9 | CC | C/T | Intermediate Metabolizer (*2) | **FALSE** — CC is ref |
| rs1057910 | CYP2C9 | AA | A/C | Poor Metabolizer (*3) | **FALSE** — AA is ref |
| rs9923231 | VKORC1 | CC | C/T | Increased Sensitivity | **FALSE** — CC is ref |
| rs3892097 | CYP2D6 | CT | G/A | Poor Metabolizer (*4) | **OK** — CT is het (note: strand?) |
| rs1045642 | ABCB1 | GG | A/G | Altered Drug Transport | **OK** — GG is hom alt |

**Root cause:** `score_pharma()` has no allele check (see audit C1).

## ClinVar Results (4 results, 2 suspicious)

| rsID | Gene | Condition | User Genotype | Ref/Alt | Stars | Assessment |
|------|------|-----------|---------------|---------|-------|------------|
| rs63750447 | APP | Alzheimer disease | TT | C/T | 3 | **SUSPICIOUS** — Hom alt for AF=0.0001 variant. Likely strand flip. |
| rs28897696 | BRCA1 | HBOC | GG | A/G | 3 | **SUSPICIOUS** — Hom alt pathogenic BRCA1. DB pos=43093449 (GRCh38?), file pos=41215920 (GRCh37). Likely strand flip. |
| rs11571833 | BRCA2 | HBOC | AA | A/T | 3 | OK — AA is ref, correctly reported as carrier 1 copy (palindromic?) |
| rs80357906 | BRCA2 | HBOC | DD | C/T | 3 | OK — DD (deletion/indel), correctly handled |

**Notes:**
- The two suspicious results (APP, BRCA1) appear to be strand-flip false positives
  where the user's genotype happens to match the ALT on the complement strand
- Both show homozygous pathogenic for extremely rare variants — biologically implausible
- The allele check IS running for ClinVar (unlike PGx), but the strand resolution may
  be incorrectly accepting complement matches as true positives

## GWAS Results (7 results, all correct)

| rsID | Trait | Genotype | Risk Allele | Copies | Palindromic? | Assessment |
|------|-------|----------|-------------|--------|-------------|------------|
| rs4402960 | Type 2 diabetes | GT | T | 1 | No | Correct |
| rs10757274 | Coronary artery disease | AG | G | 1 | No | Correct |
| rs1333049 | Coronary artery disease | GG | C | 2* | C/G palindrome | Flagged |
| rs10811661 | Type 2 diabetes | TT | T | 2* | A/T palindrome | Flagged |
| rs6265 | Cognitive function | TT | T | 2* | A/T palindrome | Flagged |
| rs1800497 | Addiction susceptibility | AG | A | 1 | No | Correct |
| rs9939609 | BMI/Obesity | TT | A | 2* | A/T palindrome | Flagged |

*Palindromic results carry strand ambiguity caveat — correct behavior.

GWAS variants correctly skipped (user has 0 copies of risk allele):
- rs429358 (APOE), rs7412 (LDL), rs7903146 (T2D), rs2187668 (Celiac), rs53576 (Social)

## Report Quality

| Feature | Status |
|---------|--------|
| Medical disclaimer | Present |
| Confidence tiers | Shown for every result |
| Data source attributions | Present (ClinVar, GWAS, gnomAD/dbSNP, PharmGKB) |
| DTC caveat per result | Present on every result |
| Limitations per result | Present (strand ambiguity, population-specific, etc.) |
| Assembly warning | "Could not determine genome assembly of the database" |
| ConfirmationUrgency | Computed but NOT rendered (see audit H3) |
| FDA PGx disclaimer | Missing (see audit H4) |

## Summary

The pipeline is architecturally sound for ClinVar and GWAS paths but has critical
false-positive issues in PGx scoring (no allele check) and potential strand-resolution
issues in ClinVar that produce biologically implausible results for rare pathogenic variants.
