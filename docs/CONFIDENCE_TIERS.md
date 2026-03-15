# GeneSight — Confidence Tier System

## Overview

Every result in GeneSight is assigned a confidence tier.
This is not optional — it is a core feature of the tool.

The tier system exists because DNA analysis is not a binary yes/no.
A BRCA1 mutation detection has a fundamentally different level of certainty than
a polygenic risk score for depression. Users need to be able to recognize this
at a glance.

---

## Tier 1: Reliable (>95% predictive value)

### What belongs here?

**Monogenic disorders** — A single gene/variant is directly causal.

| Disorder | Gene | Inheritance | ClinVar Status |
|----------|------|-------------|----------------|
| Sickle cell disease | HBB | Autosomal recessive | Pathogenic (4★) |
| Cystic fibrosis | CFTR | Autosomal recessive | Pathogenic (4★) |
| Huntington's disease | HTT | Autosomal dominant | Pathogenic (4★) |
| BRCA1/2 breast cancer risk | BRCA1, BRCA2 | Autosomal dominant | Pathogenic (4★) |
| Hemochromatosis | HFE | Autosomal recessive | Pathogenic (3-4★) |
| Factor V Leiden | F5 | Autosomal dominant | Pathogenic (4★) |

**Carrier status** — Heterozygous carriers of recessive disorders.
Clinically relevant for family planning.

**Pharmacogenomics** — How the body metabolizes medications.

| Gene | Medication(s) | PharmGKB Level | Effect |
|------|---------------|----------------|--------|
| CYP2D6 | Codeine, Tramadol, Tamoxifen | 1A | Poor → no effect; Ultrarapid → overdose risk |
| CYP2C19 | Clopidogrel, Omeprazole | 1A | Poor → reduced efficacy |
| CYP2C9 + VKORC1 | Warfarin | 1A | Dose adjustment required |
| HLA-B*5701 | Abacavir | 1A | Hypersensitivity reaction |
| DPYD | 5-Fluorouracil | 1A | Severe toxicity in poor metabolizers |

**Simple traits** — Few genes, well understood.

| Trait | Gene/SNP | Accuracy |
|-------|----------|----------|
| Lactose tolerance | MCM6 (rs4988235) | >95% |
| Earwax type | ABCC11 (rs17822931) | >95% |
| Asparagus odor detection | Multiple SNPs near OR2M7 | ~90% |
| Bitter taster (PTC) | TAS2R38 | >90% |

### Logic for Tier 1 assignment

```rust
fn is_tier1(annotation: &Annotation) -> bool {
    // ClinVar: pathogenic + review >= 2 stars
    if let Some(cv) = &annotation.clinvar {
        if cv.significance == "Pathogenic"
           && cv.review_stars >= 2 {
            return true;
        }
    }

    // PharmGKB: evidence level 1A or 1B
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("1") {
            return true;
        }
    }

    // SNPedia: magnitude >= 4 AND confirmed in ClinVar
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude >= 4.0 && annotation.clinvar.is_some() {
            return true;
        }
    }

    false
}
```

---

## Tier 2: Probable (60-85% predictive value)

### What belongs here?

**Polygenic risk scores (PRS)** — Hundreds to thousands of variants act together.
No single SNP is causal, but in combination they produce a statistical risk profile.

| Disorder/Trait | Number of SNPs | Explained variance | Predictive power |
|----------------|----------------|-------------------|------------------|
| Coronary artery disease | ~1.7M | ~15% | Top quintile: 3x risk |
| Type 2 diabetes | ~400K | ~10% | Top quintile: 2.5x risk |
| Hypertension | ~900K | ~8% | Moderately predictive |
| BMI/Obesity | ~700K | ~6% | Weakly to moderately predictive |

**Important:** "Top quintile: 3x risk" does NOT mean "you will get this."
It means: "In a population, your genetic risk is higher than 80% of people."
Lifestyle, diet, and environment are often more important.

**Physical traits involving multiple genes:**

| Trait | Accuracy |
|-------|----------|
| Hair color | ~70-85% |
| Freckles | ~70-80% |
| Baldness risk (male) | ~70-80% |
| Eye color | ~75-90% (blue vs. brown good, mixed colors weaker) |

### Logic for Tier 2 assignment

```rust
fn is_tier2(annotation: &Annotation) -> bool {
    // GWAS: significant association with moderate effect size
    if let Some(gwas) = &annotation.gwas_hits {
        let significant = gwas.iter().any(|g|
            g.p_value < 5e-8  // genome-wide significance
            && g.odds_ratio.map_or(false, |or| or > 1.1 && or < 3.0)
        );
        if significant { return true; }
    }

    // SNPedia: magnitude 2-3.9, no ClinVar entry
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude >= 2.0 && snp.magnitude < 4.0
           && annotation.clinvar.is_none() {
            return true;
        }
    }

    // PharmGKB: evidence level 2A or 2B
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("2") {
            return true;
        }
    }

    false
}
```

---

## Tier 3: Speculative (50-65% predictive value)

### What belongs here?

**Complex psychiatric disorders:**

| Disorder | Heritability (twins) | Explained by SNPs | Individual SNP effect |
|----------|---------------------|-------------------|----------------------|
| Schizophrenia | ~80% | ~7% | Minimal |
| Bipolar disorder | ~70% | ~5% | Minimal |
| Major depression | ~40% | ~2% | Tiny |
| Autism | ~50-90% | ~3% | Minimal |
| ADHD | ~70-80% | ~3% | Minimal |

The high heritability in twin studies vs. low explained variance by
known SNPs = "missing heritability." Gene-environment interactions, epigenetics,
and rare variants play major roles.

**Personality and cognition:**

| Trait | SNP effect | Usefulness of testing |
|-------|-----------|----------------------|
| Intelligence (IQ) | ~0.01% per SNP | Practically none |
| Risk-taking | Minimal | Practically none |
| Neuroticism | Minimal | Practically none |

**Athletic aptitude:**

| Gene/SNP | Claim | Reality |
|----------|-------|---------|
| ACTN3 (rs1815739) | "Sprinter gene" | Explains <1% of variance in athletic performance |
| ACE I/D | Endurance vs. strength | Inconsistent study results |

### Logic for Tier 3 assignment

```rust
fn is_tier3(annotation: &Annotation) -> bool {
    // GWAS: weak association
    if let Some(gwas) = &annotation.gwas_hits {
        let weak = gwas.iter().any(|g|
            g.p_value < 5e-8
            && g.odds_ratio.map_or(true, |or| or <= 1.1)
        );
        if weak { return true; }
    }

    // SNPedia: low magnitude
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude > 0.0 && snp.magnitude < 2.0 {
            return true;
        }
    }

    // PharmGKB: evidence level 3 or 4
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("3")
           || pgkb.evidence_level.starts_with("4") {
            return true;
        }
    }

    false
}
```

---

## Report Presentation

### CLI Output (Example)

```
═══════════════════════════════════════════════════════════
  GeneSight Report — 2026-03-14
  File: my_dna.txt (23andMe format)
  Variants analyzed: 637,291
  Annotated variants: 12,847
═══════════════════════════════════════════════════════════

🟢 TIER 1 — Reliable (clinically validated)
───────────────────────────────────────────────────────────

  PHARMACOGENOMICS
  ┌─────────────┬──────────────┬───────────────────────────┐
  │ Gene        │ Medication   │ Status                    │
  ├─────────────┼──────────────┼───────────────────────────┤
  │ CYP2D6      │ Codeine      │ ⚠ Poor Metabolizer        │
  │ CYP2C19     │ Clopidogrel  │ ✓ Normal Metabolizer      │
  └─────────────┴──────────────┴───────────────────────────┘
  Source: PharmGKB (Level 1A) + ClinVar

  CARRIER STATUS
  • Cystic fibrosis (CFTR): Heterozygous — carrier, not affected
    Source: ClinVar (4★ review)

  TRAITS
  • Lactose tolerance: Likely lactose intolerant (CC at rs4988235)
  • Earwax: Dry type

🟡 TIER 2 — Probable (statistical association)
───────────────────────────────────────────────────────────

  POLYGENIC RISK SCORES
  • Coronary artery disease: 72nd percentile (slightly elevated)
    ⚠ Lifestyle has a greater impact than genetics
  • Type 2 diabetes: 45th percentile (average)

  TRAITS
  • Baldness risk: Elevated (68% probability by age 50)
  • Freckles: Likely present

🔴 TIER 3 — Speculative (weak evidence)
───────────────────────────────────────────────────────────

  ⚠ The following results have low predictive value
    and should NOT be used for decision-making.

  • ACTN3 (muscle type): Mixed type (endurance/strength)
  • Caffeine metabolism: Fast (CYP1A2 rs762551 AA)

───────────────────────────────────────────────────────────
⚕ DISCLAIMER: This report is informational, not
  diagnostic. For medical decisions, consult a
  physician or genetic counselor.
───────────────────────────────────────────────────────────
```

---

## What GeneSight deliberately does NOT show

1. **No lifespan predictions**
2. **No IQ scores or "intelligence genes"**
3. **No race-/ethnicity-based risk assessments** without context
4. **No results without a confidence tier**
5. **No results that sound like diagnoses** — always phrased probabilistically
