# H1: No Mode-of-Inheritance (AD/AR) in ClinVar Scoring

**Severity:** HIGH
**Status:** Not implemented — `CarrierStatus` category exists but is never assigned

## Problem

The ClinVar scorer treats all pathogenic variants identically regardless of zygosity and
inheritance pattern. A heterozygous carrier of an **autosomal recessive** (AR) pathogenic
variant is reported the same way as a heterozygous carrier of an **autosomal dominant** (AD)
pathogenic variant — both as `MonogenicDisease` with `Tier1Reliable`.

This is scientifically wrong:
- **AD pathogenic + heterozygous** = potentially affected (report as MonogenicDisease)
- **AR pathogenic + heterozygous** = carrier only, NOT affected (report as CarrierStatus)
- **AR pathogenic + homozygous alt** = potentially affected (report as MonogenicDisease)

## Evidence

The `ResultCategory::CarrierStatus` variant exists at `models/report.rs:52`:
```rust
pub enum ResultCategory {
    MonogenicDisease,
    CarrierStatus,      // ← exists but never assigned
    Pharmacogenomics,
    ...
}
```

In `scorer/mod.rs`, the `score_clinvar()` function (lines 200-461):
- Checks copy count (line 230-237) — correctly determines 0/1/2 ALT copies
- Assigns `MonogenicDisease` category for all pathogenic results regardless of copy count
- Never queries or considers mode of inheritance
- Never assigns `CarrierStatus`

There is no MOI field in `ClinVarAnnotation` (`models/annotation.rs`), no MOI lookup
from ClinVar or MedGen, and no inheritance-aware logic anywhere.

## Example of Wrong Behavior

A user heterozygous for **CFTR p.F508del** (the most common CF mutation):
- **Current output**: `MonogenicDisease`, `Tier1Reliable`, `HighImpact`
- **Correct output**: `CarrierStatus`, `Tier1Reliable`, `ClinicalConfirmationRecommended`
  with text: "You are a carrier of this autosomal recessive condition. Carriers typically
  do not develop symptoms but may pass the variant to offspring."

## Scientific Requirement

From the research report (Section: ClinVar-Pathogenität korrekt interpretieren):

> Genotyp (0/1/2 ALT-Kopien) ist für monogene Erkrankungen nur im Kontext der
> **Vererbung (AD/AR/X-linked/mitochondrial)** interpretierbar.

> MOI/zygosity ignorieren → AR-Erkrankungen fälschlich als erkrankt bei heterozygoter
> Trägerschaft.

Programmable MOI data sources identified by the research:
- ClinVar submissions may contain MOI
- ClinVar properties/filters contain MOI categories
- MedGen supports MOI as a property field
- ClinGen Gene-Disease Validity KB has MOI per gene-disease assertion

## Fix Requirements

1. Add MOI data source (simplest: curate a gene→MOI lookup table for common genes,
   or parse from ClinVar XML/properties)
2. Add `mode_of_inheritance` field to `ClinVarAnnotation`
3. In `score_clinvar()`: when pathogenic + het + AR → assign `CarrierStatus` category
4. Adjust `ConfirmationUrgency` for carrier findings
5. Add carrier-specific language to report rendering
