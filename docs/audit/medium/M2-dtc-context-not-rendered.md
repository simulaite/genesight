# M2: `Report.dtc_context` Populated but Not Rendered

**Severity:** MEDIUM
**Status:** Data populated in `lib.rs:247`; neither renderer outputs it

## Problem

The `DTC_CONTEXT` constant at `lib.rs:80-85` contains important context about DTC
microarray limitations:

> "This report is based on direct-to-consumer (DTC) microarray genotyping data. DTC
> genotyping has known limitations including strand ambiguity, limited coverage of the
> genome, and potential genotyping errors. Any clinically relevant finding should be
> confirmed through clinical-grade testing (e.g., Sanger sequencing or clinical NGS)
> before making medical decisions."

This is assigned to `Report.dtc_context` at `lib.rs:247` and serialized to JSON output.
However, neither `report/html.rs` nor `report/markdown.rs` renders this field.

The medical `DISCLAIMER` IS rendered by both renderers via `write_disclaimer()`. But
`dtc_context` is a separate, complementary statement specifically about the DTC data
source limitations.

## Fix Requirements

1. Add `write_dtc_context()` function to both renderers
2. Render as a distinct section: "About Your Data Source" or "Understanding DTC Data"
3. Position after the medical disclaimer and before results
