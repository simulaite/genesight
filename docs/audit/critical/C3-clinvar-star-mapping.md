# C3: ClinVar Import Maps Review Stars Off-by-One

**Severity:** CRITICAL
**Status:** All imported ClinVar data has inflated review stars

## Problem

The ClinVar import script at `data/import/import_clinvar.py:38-47` maps review status
text to star ratings that are uniformly **one level too high**.

### Current (WRONG) Mapping

```python
REVIEW_STARS = {
    "practice guideline": 4,                                          # Correct
    "reviewed by expert panel": 4,                                    # WRONG: should be 3
    "criteria provided, multiple submitters, no conflicts": 3,        # WRONG: should be 2
    "criteria provided, conflicting classifications": 2,              # WRONG: should be 1
    "criteria provided, conflicting interpretations of pathogenicity": 2,  # WRONG: should be 1
    "criteria provided, single submitter": 2,                         # WRONG: should be 1
    "no assertion criteria provided": 1,                              # WRONG: should be 0
    "no assertion provided": 0,                                       # Correct
}
```

### Correct ClinVar Star Rating (per NCBI documentation)

```
4 stars = practice guideline
3 stars = reviewed by expert panel
2 stars = criteria provided, multiple submitters, no conflicts
1 star  = criteria provided, single submitter
1 star  = criteria provided, conflicting classifications
0 stars = no assertion criteria provided
0 stars = no assertion provided
```

Source: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

## Impact

The scorer at `scorer/mod.rs:420` assigns Tier1 when `review_stars >= 2`. With the
off-by-one mapping:

| ClinVar Status | Correct Stars | Imported Stars | Correct Tier | Actual Tier |
|----------------|---------------|----------------|-------------|-------------|
| Single submitter, pathogenic | 1 | **2** | Tier2 | **Tier1** |
| No assertion criteria, pathogenic | 0 | **1** | Tier3/skip | **Tier2** |
| Conflicting classifications | 1 | **2** | Tier3 (conflict) | Text-match saves it* |
| Expert panel | 3 | **4** | Tier1 | Tier1 (no impact) |

*The `is_conflicting` text check at scorer line 214 catches "conflicting" in the
significance string, compensating for the star inflation in this one case.

**Single-submitter pathogenic entries are the most dangerous promotion**: they are the
most common category in ClinVar and have the highest rate of eventual reclassification.
Promoting them to Tier1 (same authority as expert panel) directly contradicts the
research requirement.

## Scientific Requirement

From the research report (Section: Correctly Interpreting ClinVar Pathogenicity):

> **High confidence**: 3-4 stars (Expert Panel/Practice Guideline).
> **Medium**: 2 stars (multiple submitters, criteria provided, no conflicts).
> **Low/Informational**: 1 star (single submitter with criteria) — display, but
> clearly label as non-consensus.
> **0 stars / no criteria**: do not present as a clinical assertion; at most as a raw hint.

## Fix Requirements

```python
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations of pathogenicity": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
}
```

After fixing the mapping, the seed database must be rebuilt. Any existing `genesight.db`
files will have wrong star values until reimported.
