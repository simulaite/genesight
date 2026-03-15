# H8: Assembly Mismatch Warns but Doesn't Block

**Severity:** HIGH
**Status:** Detection implemented; enforcement not implemented

## Problem

The pipeline detects assembly mismatches (GRCh37 vs GRCh38) and generates warnings in
`lib.rs:260-293`, but continues processing and produces results from potentially incorrect
coordinate-based lookups.

Current behavior when GRCh37 input file is analyzed against GRCh38 database:
1. Warning text generated: "Assembly mismatch: input file uses GRCh37 but database uses GRCh38"
2. Warning appended to `Report.assembly_warnings`
3. **Pipeline continues normally and produces results**
4. rsID-based lookups still work (rsIDs are assembly-independent)
5. Position-based lookups would silently return wrong results

## Live Test Evidence

The huAE4518 test showed:
- Input file: GRCh37 (23andMe v4)
- Database: Unknown assembly (no `db_metadata` table populated)
- Warning generated: "Could not determine genome assembly of the database"
- Pipeline produced results anyway

The `variants` table in the seed database appears to use GRCh38 positions (e.g.,
rs28897696 at position 43093449) while the 23andMe file reports GRCh37 position 41215920.
Since lookup is by rsID, this didn't break — but it could if position-based matching
is ever added.

## Scientific Requirement

From the research report (Section: Allele Matching):

> Build mix (GRCh37 vs GRCh38) causes apparent 'allele mismatch' and incorrect
> annotations; liftover often works very well, but can fail for certain variant
> types/regions, therefore 'unmappable' must be handled explicitly.

## Fix Requirements

For Phase 1 (current):
1. Store assembly in `db_metadata` table during import/seed build
2. When assemblies are known and different, add assembly warning to EVERY scored result's
   `limitations` field (not just the report header)
3. Consider blocking analysis entirely when assemblies are incompatible and user hasn't
   passed `--force-assembly` flag

For Phase 2 (future):
4. Implement LiftOver using chain files
5. Track `unmappable` variants that can't be lifted
6. Add per-variant assembly provenance
