# H2: ClinVar Import Never Populates `classification_type`

**Severity:** HIGH
**Status:** Feature architecturally complete but data never populated

## Problem

The codebase supports germline/somatic/oncogenicity distinction at every layer EXCEPT
the import pipeline:

| Layer | Status |
|-------|--------|
| Schema (`data/schema/schema.sql:40`) | Has `classification_type TEXT DEFAULT 'germline'` |
| DB adapter (`db/clinvar.rs:53-63`) | Reads `classification_type`, prefers germline |
| Model (`models/annotation.rs:36-43`) | Has `ClinVarClassificationType` enum |
| Scorer (`scorer/mod.rs:258-303`) | Demotes somatic/oncogenicity to Tier3 |
| **Import (`import_clinvar.py`)** | **Does NOT read or write `classification_type`** |

The import script inserts rows without specifying `classification_type`, so all entries
get the SQLite default `'germline'`. Somatic and oncogenicity classifications from ClinVar
are imported as germline.

## Impact

Since ClinVar 2024, the same variant can have separate germline and somatic classifications.
A variant classified as "Pathogenic" in a somatic (tumor) context may be benign or VUS in
germline context. Without importing this distinction, all somatic-pathogenic variants are
treated as germline-pathogenic, producing false clinical alerts for consumer germline reports.

## Scientific Requirement

From the research report (Section: ClinVar-Pathogenität korrekt interpretieren):

> Seit 2024 trennt ClinVar klinische Klassifikationstypen (germline, somatic clinical
> impact, oncogenicity) in getrennte Felder; 'clinical_significance' muss daher
> kontextualisiert werden.

> Somatische/onko-Klassifikation als germline 'Pathogenic' ausgeben; ClinVar trennt
> diese Klassifikationstypen.

## Fix Requirements

1. In `import_clinvar.py`: Read the classification type column from `variant_summary.txt`
   (ClinVar's full release has separate columns for germline and somatic significance
   since the 2024 format change)
2. Write `classification_type` to the `clinvar` table during import
3. Rebuild seed database with correct classification types

## ClinVar variant_summary.txt Columns

The variant_summary file includes:
- Column 6 (`ClinicalSignificance`): May contain germline, somatic, or mixed
- Column 7 (`ClinSigSimple`): Simplified version
- The newer format has `ClinicalSignificance` split by classification type

The import script should parse this to determine germline vs somatic provenance.
