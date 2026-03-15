# H7: Polygenic Risk Scores Not Implemented

**Severity:** HIGH
**Status:** `scorer/polygenic.rs` is an empty stub; "PRS" label applied to single SNPs

## Problem

The `scorer/polygenic.rs` module contains only a one-line doc comment and no implementation.
There is no PRS computation anywhere in the codebase.

Despite this, the scorer assigns `ResultCategory::PolygenicRiskScore` to individual GWAS
hits when `p < 5e-8 AND OR > 1.5` (scorer/mod.rs:592-606). This is misleading — a single
SNP is not a polygenic risk score. A PRS is a weighted aggregate across many variants.

## Scientific Requirement

From the research report (Section: GWAS Risk Interpretation):

> Standard PRS is a weighted sum score: PRS = Σ(β_i * G_i), where G_i is the number
> of effect alleles. For interpretable percentiles, a reference distribution
> (Mean/SD) in the matching ancestry cohort is needed.

Requirements for a proper PRS:
1. Select validated weight set (betas) for a specific trait/disease
2. Count effect alleles (0/1/2) per variant after strand normalization
3. Compute weighted sum: PRS = Σ(β_i × G_i)
4. Compare against reference distribution for population-matched percentile
5. Report percentile, not raw score

## Current Mislabeling

A user seeing "Polygenic Risk Score" for a single SNP with OR=1.6 will think they have
received a comprehensive genetic risk assessment, when they have received a single
data point with modest predictive value.

## Fix Options

**Option A (Minimum — rename):** Rename `PolygenicRiskScore` to `GwasAssociation` or
`SingleVariantAssociation`. Accurate labeling without implementing PRS.

**Option B (Full implementation):** Implement PRS computation with:
- Curated weight tables per trait (from PGS Catalog or published GWAS)
- Weighted sum computation
- Reference distribution for at least European ancestry
- Percentile reporting
- Ancestry caveat

Option A is recommended for now; Option B is a Phase 2 feature requiring significant
data infrastructure.
