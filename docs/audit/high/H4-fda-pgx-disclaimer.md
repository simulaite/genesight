# H4: No FDA PGx Disclaimer Anywhere in Codebase

**Severity:** HIGH
**Status:** Not implemented — zero matches for "FDA" in entire Rust codebase

## Problem

The research report requires specific FDA disclaimer language for pharmacogenomic results.
No such language exists anywhere in the codebase.

The `score_pharma()` function at `scorer/mod.rs:469-510` only appends the generic
`DTC_RAW_DATA_CAVEAT` to PGx results. The `DISCLAIMER` constant in `lib.rs:58-65`
is a general medical disclaimer with no PGx-specific or FDA-specific content.

## Scientific Requirement

From the research report (Section: Verantwortungsvolle Ergebnisdarstellung):

> Die FDA warnte öffentlich, dass viele PGx-Testclaims (inkl. direkt an Konsumenten
> vermarkteter Tests und Software-Interpretationen) nicht von der FDA geprüft sind und
> wissenschaftlich unzureichend gestützt sein können; Therapieänderungen auf Basis
> solcher Claims können Patientenschaden verursachen.

From the PIPELINE_OVERHAUL.md (Blueprint 7):

> "Pharmacogenomic results from consumer genotyping arrays have NOT been reviewed
> or approved by the FDA for clinical use. Do not alter any medication regimen
> based solely on these results. Consult a healthcare provider or clinical
> pharmacogenomics service for validated testing."

The FDA Table of Pharmacogenetic Associations repeats that genotyping does not replace
clinical vigilance and patient management.

## Fix Requirements

1. Add `PGX_FDA_DISCLAIMER` constant to `scorer/mod.rs` or `lib.rs`:
   ```rust
   const PGX_FDA_DISCLAIMER: &str = "Pharmacogenomic results from consumer genotyping \
       arrays have NOT been reviewed or approved by the U.S. Food and Drug Administration \
       (FDA) for clinical use. Do not alter any medication regimen based solely on these \
       results. Consult a healthcare provider or clinical pharmacogenomics service for \
       validated testing.";
   ```

2. Append this to `limitations` in `score_pharma()` (and any future PGx scoring path)

3. Add a dedicated PGx disclaimer section in both HTML and Markdown report renderers,
   rendered only when the report contains PGx results

4. Reference the FDA Table of Pharmacogenetic Associations URL as a resource link
