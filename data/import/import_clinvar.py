#!/usr/bin/env python3
"""Import ClinVar variant_summary.txt into a GeneSight SQLite database.

Parses the tab-separated ClinVar variant_summary file, filters for GRCh38
rows with valid rsIDs, and populates the `variants` and `clinvar` tables.

Usage:
    python3 import_clinvar.py --input data/raw/clinvar/variant_summary.txt \
                              --output genesight.db
"""

import argparse
import json
import sqlite3
import sys
import time

# ---------------------------------------------------------------------------
# Column indices (0-based) in variant_summary.txt
# ---------------------------------------------------------------------------
COL_GENE_SYMBOL = 4
COL_CLINICAL_SIGNIFICANCE = 6
COL_LAST_EVALUATED = 8
COL_RS_DBSNP = 9
COL_PHENOTYPE_LIST = 13
COL_ASSEMBLY = 16
COL_CHROMOSOME = 18
COL_START = 19
COL_REF_ALLELE = 21
COL_ALT_ALLELE = 22
COL_REVIEW_STATUS = 24

# ---------------------------------------------------------------------------
# Review-status text → star rating
# ---------------------------------------------------------------------------
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 4,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, conflicting classifications": 2,
    "criteria provided, conflicting interpretations of pathogenicity": 2,
    "criteria provided, single submitter": 2,
    "no assertion criteria provided": 1,
    "no assertion provided": 0,
}

# ---------------------------------------------------------------------------
# Schema DDL
# ---------------------------------------------------------------------------
SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS variants (
    rsid TEXT PRIMARY KEY,
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref_allele TEXT,
    alt_allele TEXT
);
CREATE INDEX IF NOT EXISTS idx_variants_chr_pos ON variants(chromosome, position);

CREATE TABLE IF NOT EXISTS clinvar (
    rsid TEXT REFERENCES variants(rsid),
    clinical_significance TEXT,
    review_status INTEGER,
    conditions TEXT,
    gene_symbol TEXT,
    last_updated DATE
);
CREATE INDEX IF NOT EXISTS idx_clinvar_rsid ON clinvar(rsid);
"""


def review_status_to_stars(status_text: str) -> int:
    """Map ClinVar ReviewStatus text to a 0-4 star rating."""
    return REVIEW_STARS.get(status_text.strip().lower(), 0)


def parse_conditions(phenotype_list: str) -> str:
    """Split pipe-separated PhenotypeList into a JSON array string."""
    if not phenotype_list or phenotype_list == "-":
        return "[]"
    conditions = [c.strip() for c in phenotype_list.split("|") if c.strip()]
    return json.dumps(conditions, ensure_ascii=False)


def create_schema(conn: sqlite3.Connection) -> None:
    """Create the variants and clinvar tables (idempotent)."""
    conn.executescript(SCHEMA_SQL)


def import_clinvar(input_path: str, output_path: str) -> None:
    """Read variant_summary.txt and import qualifying rows into SQLite."""
    conn = sqlite3.connect(output_path)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA cache_size = -64000")  # 64 MB cache

    create_schema(conn)

    # We accumulate the best record per rsid in memory (dict keyed by rsid).
    # For ~8.9M rows only a fraction have valid rsIDs on GRCh38, so the dict
    # stays manageable.  We keep the one with the highest star rating.
    best: dict[str, tuple] = {}
    # Each value: (rsid, chromosome, position, ref_allele, alt_allele,
    #              clinical_significance, stars, conditions_json,
    #              gene_symbol, last_evaluated)

    total_rows = 0
    skipped = 0
    kept = 0

    t0 = time.monotonic()

    with open(input_path, "r", encoding="utf-8", errors="replace") as fh:
        # Skip the header line
        header = fh.readline()
        if not header.startswith("#"):
            print(
                "WARNING: Expected header starting with '#', got:",
                header[:80],
                file=sys.stderr,
            )

        for line in fh:
            total_rows += 1

            if total_rows % 500_000 == 0:
                elapsed = time.monotonic() - t0
                print(
                    f"  ... processed {total_rows:,} rows "
                    f"({elapsed:.1f}s, {kept:,} kept so far)"
                )

            fields = line.rstrip("\n").split("\t")

            # Ensure we have enough columns
            if len(fields) <= COL_REVIEW_STATUS:
                skipped += 1
                continue

            # --- Filter: GRCh38 only ---
            assembly = fields[COL_ASSEMBLY].strip()
            if assembly != "GRCh38":
                skipped += 1
                continue

            # --- Filter: valid rsID ---
            rs_raw = fields[COL_RS_DBSNP].strip()
            if rs_raw == "-1" or rs_raw == "" or rs_raw == "-":
                skipped += 1
                continue

            # --- Filter: non-empty clinical significance ---
            clin_sig = fields[COL_CLINICAL_SIGNIFICANCE].strip()
            if not clin_sig:
                skipped += 1
                continue

            rsid = "rs" + rs_raw

            # Parse remaining fields
            chromosome = fields[COL_CHROMOSOME].strip()
            try:
                position = int(fields[COL_START].strip())
            except ValueError:
                skipped += 1
                continue

            ref_allele = fields[COL_REF_ALLELE].strip()
            alt_allele = fields[COL_ALT_ALLELE].strip()
            # Normalise missing alleles
            if ref_allele in ("na", "-", ""):
                ref_allele = None
            if alt_allele in ("na", "-", ""):
                alt_allele = None

            review_text = fields[COL_REVIEW_STATUS].strip()
            stars = review_status_to_stars(review_text)

            conditions_json = parse_conditions(fields[COL_PHENOTYPE_LIST])
            gene_symbol = fields[COL_GENE_SYMBOL].strip() or None
            last_evaluated = fields[COL_LAST_EVALUATED].strip() or None

            # Keep the row with the highest star rating per rsid
            if rsid in best:
                existing_stars = best[rsid][6]
                if stars <= existing_stars:
                    kept += 1  # counted as processed/kept but superseded
                    continue

            best[rsid] = (
                rsid,
                chromosome,
                position,
                ref_allele,
                alt_allele,
                clin_sig,
                stars,
                conditions_json,
                gene_symbol,
                last_evaluated,
            )
            kept += 1

    elapsed_parse = time.monotonic() - t0
    print(
        f"\nParsing complete in {elapsed_parse:.1f}s: "
        f"{total_rows:,} total rows, {skipped:,} skipped, "
        f"{len(best):,} unique rsIDs to import."
    )

    # ------------------------------------------------------------------
    # Batch insert into SQLite
    # ------------------------------------------------------------------
    BATCH_SIZE = 50_000
    records = list(best.values())
    best.clear()  # free memory

    print(f"Inserting {len(records):,} records into database ...")

    conn.execute("BEGIN")
    for i in range(0, len(records), BATCH_SIZE):
        batch = records[i : i + BATCH_SIZE]

        # Variants table
        conn.executemany(
            """
            INSERT OR REPLACE INTO variants
                (rsid, chromosome, position, ref_allele, alt_allele)
            VALUES (?, ?, ?, ?, ?)
            """,
            [(r[0], r[1], r[2], r[3], r[4]) for r in batch],
        )

        # ClinVar table
        conn.executemany(
            """
            INSERT OR REPLACE INTO clinvar
                (rsid, clinical_significance, review_status,
                 conditions, gene_symbol, last_updated)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [(r[0], r[5], r[6], r[7], r[8], r[9]) for r in batch],
        )

        if (i // BATCH_SIZE + 1) % 5 == 0 or i + BATCH_SIZE >= len(records):
            print(f"  ... inserted {min(i + BATCH_SIZE, len(records)):,} / {len(records):,}")

    conn.commit()

    # ------------------------------------------------------------------
    # Final statistics
    # ------------------------------------------------------------------
    (variant_count,) = conn.execute("SELECT COUNT(*) FROM variants").fetchone()
    (clinvar_count,) = conn.execute("SELECT COUNT(*) FROM clinvar").fetchone()

    elapsed_total = time.monotonic() - t0
    print(f"\nDone in {elapsed_total:.1f}s.")
    print(f"  Total rows processed:  {total_rows:,}")
    print(f"  Rows skipped:          {skipped:,}")
    print(f"  Rows imported (clinvar): {clinvar_count:,}")
    print(f"  Unique rsIDs (variants): {variant_count:,}")

    conn.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Import ClinVar variant_summary.txt into GeneSight SQLite DB."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to ClinVar variant_summary.txt",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output SQLite database (e.g. genesight.db)",
    )
    args = parser.parse_args()

    print(f"Importing ClinVar from: {args.input}")
    print(f"Output database:        {args.output}")
    print()

    import_clinvar(args.input, args.output)


if __name__ == "__main__":
    main()
