# H5: OR Presented Without Relative-Risk Explanation

**Severity:** HIGH
**Status:** Effect context thresholds implemented; no explanation that OR ≠ absolute risk

## Problem

The `gwas_effect_context()` function at `scorer/mod.rs:725-750` produces qualitative
context strings ("Substantially elevated risk", "Moderately elevated risk", etc.) based
on OR thresholds. This is correct. However:

1. **No explanation that OR is relative, not absolute**: The OR appears in result details
   as `"odds ratio {or:.2}"` without any caveat that this is a relative measure, not a
   probability.

2. **No baseline prevalence**: Without disease prevalence (p0), absolute risk cannot be
   computed. The research provides the formula: `odds = (p0/(1-p0)) * OR^k; p = odds/(1+odds)`.
   No prevalence data exists in the GWAS table or scorer.

3. **Beta coefficients**: For quantitative traits, beta represents additive effect per
   allele on the trait scale. No unit or population reference (mean/SD) is provided to
   contextualize the magnitude.

## Scientific Requirement

From the research report (Section: GWAS-Risikointerpretation):

> OR ist ein Odds-Verhältnis, nicht direkt 'Risiko'. Für eine absolute Risikoapproximation
> brauchst du eine Baseline-Risikoannahme p0. Ohne p0 ist 'absolutes Risiko' nicht
> seriös berechenbar.

> OR als 'absolutes Risiko' ausgeben ohne Baseline p0 und ohne RR/OR-Unterschied zu
> erklären.

## Example of Misleading Output

Current output for a GWAS hit with OR=1.3:
```
Odds ratio: 1.30. Slightly elevated risk.
```

What a user reads: "I have a slightly elevated risk" (interprets as absolute).

What it should say:
```
Odds ratio: 1.30 (relative measure). This means your odds of developing this condition
are approximately 1.3x compared to someone without this allele. The absolute risk
depends on the baseline prevalence in your population, which is not available in this
report. Consult a healthcare provider for personalized risk assessment.
```

## Fix Requirements

1. Add OR/Beta explanation caveat to all GWAS result limitations
2. Add baseline prevalence column to GWAS schema (optional, for future use)
3. Clarify in result details that OR is a relative measure
4. For beta coefficients, note the trait unit when available
5. Consider adding a brief "How to read GWAS results" section in the report
