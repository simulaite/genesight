# H3: `ConfirmationUrgency` Computed but Never Rendered

**Severity:** HIGH
**Status:** Data computed and stored on every result; both renderers ignore it

## Problem

The `ConfirmationUrgency` enum is correctly defined, correctly assigned to every
`ScoredResult`, and serialized to JSON — but neither the HTML nor Markdown renderer
reads or displays it.

### Correctly Implemented

- Enum at `models/report.rs:10-21`: `HighImpact`, `ClinicalConfirmationRecommended`,
  `InformationalOnly`
- Assignment in scorer:
  - ACMG SF genes (BRCA1/2, MLH1, etc.) → `HighImpact` (line 399-400)
  - Non-ACMG ClinVar P/LP and PGx → `ClinicalConfirmationRecommended` (line 401-402, 497)
  - GWAS, SNPedia low-magnitude, VUS, conflicting → `InformationalOnly`

### Not Rendered

- `report/html.rs`: No reference to `confirmation_urgency` in any rendering function.
  The HTML `write_results` function (lines 701-763) renders tier badge, category, summary,
  details, and limitations — but no urgency indicator.
- `report/markdown.rs`: Same — `confirmation_urgency` is never referenced in output.

## Impact

A BRCA1 pathogenic finding (which should have a prominent **"CLINICAL CONFIRMATION
STRONGLY RECOMMENDED"** banner) is visually identical to a low-magnitude SNPedia trivia
hit. The most critical safety signal in the entire system is silently dropped.

## Scientific Requirement

From the research report (Section: Verantwortungsvolle Ergebnisdarstellung):

> Für sehr folgenreiche, medizinisch 'actionable' Gene/Krankheitsbilder existiert in der
> klinischen Genomik die ACMG Secondary Findings (SF) Policy.

> Wenn eine Variante als hochrelevant erscheint, muss sie als 'High impact, confirm
> clinically' ausgegeben werden, nicht als Diagnose.

## Fix Requirements

### HTML Renderer

Add urgency badge rendering in `write_results()`:
- `HighImpact` → Red banner: "CLINICAL CONFIRMATION STRONGLY RECOMMENDED"
- `ClinicalConfirmationRecommended` → Orange badge: "Clinical confirmation recommended"
- `InformationalOnly` → Gray badge: "Informational only"

### Markdown Renderer

Add urgency indicator:
- `HighImpact` → `> **!! CLINICAL CONFIRMATION STRONGLY RECOMMENDED !!**`
- `ClinicalConfirmationRecommended` → `> *Clinical confirmation recommended*`
- `InformationalOnly` → no extra text

### Additionally

`Report.dtc_context` (populated at `lib.rs:247`) is also not rendered by either renderer.
Add a "Understanding This Report" section that renders `dtc_context` text.
