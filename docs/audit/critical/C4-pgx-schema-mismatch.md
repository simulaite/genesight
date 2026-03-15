# C4: `definitions.rs` Queries Non-Existent Column Names

**Severity:** CRITICAL
**Status:** New PGx pipeline will crash on real database

## Problem

There are two separate code paths that query the `pgx_allele_definitions` table, and
they use **different column names**:

### Path A: `pgx/mod.rs:81-84` (StarAlleleCaller — matches actual schema)

```rust
"SELECT gene, allele_name, rsid, alt_allele, function, activity_score \
 FROM pgx_allele_definitions ORDER BY gene, allele_name"
```

### Path B: `pgx/definitions.rs:63-66` (load_allele_definitions — WRONG column names)

```rust
"SELECT gene, star_allele, defining_rsid, variant_allele, function, activity_value \
 FROM pgx_allele_definitions ORDER BY gene, star_allele"
```

### Actual Schema (`data/schema/schema.sql` and `data/seed/build_seed_db.py`)

```sql
CREATE TABLE pgx_allele_definitions (
    gene TEXT NOT NULL,
    allele_name TEXT NOT NULL,      -- definitions.rs expects "star_allele"
    rsid TEXT NOT NULL,             -- definitions.rs expects "defining_rsid"
    alt_allele TEXT NOT NULL,       -- definitions.rs expects "variant_allele"
    function TEXT NOT NULL,
    activity_score REAL NOT NULL    -- definitions.rs expects "activity_value"
);
```

### Column Name Mismatches

| Actual Column | `mod.rs` (correct) | `definitions.rs` (wrong) |
|---------------|-------------------|--------------------------|
| `allele_name` | `allele_name` | `star_allele` |
| `rsid` | `rsid` | `defining_rsid` |
| `alt_allele` | `alt_allele` | `variant_allele` |
| `activity_score` | `activity_score` | `activity_value` |

## Impact

If `load_allele_definitions()` is ever called against the real seed database, SQLite will
return an error: `no such column: star_allele`. This means the entire new PGx pipeline
(diplotype calling, phasing detection, coverage tracking) cannot execute.

The function currently only runs in its own unit tests, which create tables with the
*wrong* column names matching the code (see `definitions.rs:203-211`). The tests pass
because they test against a schema that doesn't match production.

## Relationship to Other Issues

- **Blocks C2**: Cannot wire `definitions.rs` into the live pipeline until columns match
- **Root cause**: `definitions.rs` and `mod.rs` were written independently with different
  naming conventions and never integration-tested against the same database

## Fix Requirements

Either:
**Option A**: Fix `definitions.rs` to use the actual column names (`allele_name`, `rsid`,
`alt_allele`, `activity_score`) and update the struct field names to match.

**Option B**: Add column aliases to the schema. Not recommended — creates maintenance burden.

Option A is clearly better. Also fix the unit test table creation in `definitions.rs:203-211`
to use the correct schema.
