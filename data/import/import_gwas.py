#!/usr/bin/env python3
"""Import GWAS Catalog associations TSV into the GeneSight SQLite database.

Parses the full GWAS Catalog download (tab-separated) and populates the
`variants` and `gwas` tables. Only standard-library modules are used.

Usage:
    python3 import_gwas.py --input data/raw/gwas/gwas-catalog-download-associations-v1.0-full.tsv \
                           --db genesight.db
"""

import argparse
import csv
import re
import sqlite3
import sys
import time

# Column indices (0-based) from the GWAS Catalog TSV header row.
COL_PUBMEDID = 1
COL_DISEASE_TRAIT = 7
COL_CHR_ID = 11
COL_CHR_POS = 12
COL_MAPPED_GENE = 14
COL_STRONGEST_SNP_RISK_ALLELE = 20
COL_SNPS = 21
COL_RISK_ALLELE_FREQUENCY = 26
COL_P_VALUE = 27
COL_OR_BETA = 30
COL_CI_TEXT = 31

# Minimum number of columns we expect per row.
MIN_COLUMNS = 32

# Regex for a valid rsID: "rs" followed by one or more digits.
RSID_RE = re.compile(r"^rs\d+$")

# Batch size for database inserts.
BATCH_SIZE = 5000

# Progress reporting interval (input rows).
PROGRESS_INTERVAL = 200_000


def parse_float(value: str):
    """Parse a string as float, returning None on failure or empty/NR values."""
    if not value or value.strip().upper() in ("", "NR", "NA", "NAN", "-"):
        return None
    try:
        return float(value.strip())
    except ValueError:
        return None


def parse_position(value: str):
    """Parse chromosome position as integer. Returns None if unparseable.

    Some rows have multiple positions separated by semicolons or 'x' —
    take the first valid integer.
    """
    if not value or not value.strip():
        return None
    # Split on common delimiters used in the GWAS catalog for multiple positions.
    for part in re.split(r"[;x,\s]+", value.strip()):
        try:
            return int(part)
        except ValueError:
            continue
    return None


def parse_risk_allele(snp_risk_allele: str, target_rsid: str):
    """Extract the risk allele letter(s) from the STRONGEST SNP-RISK ALLELE field.

    The format is typically "rs12345-A" or "rs12345-?".  When multiple entries
    are present (semicolon-separated), find the one matching target_rsid.
    Returns None if no allele can be determined.
    """
    if not snp_risk_allele or not snp_risk_allele.strip():
        return None

    entries = [e.strip() for e in snp_risk_allele.split(";")]
    for entry in entries:
        if "-" not in entry:
            continue
        parts = entry.rsplit("-", 1)
        rsid_part = parts[0].strip()
        allele_part = parts[1].strip() if len(parts) > 1 else None

        if rsid_part == target_rsid and allele_part and allele_part != "?":
            return allele_part

    # Fallback: if there is only one entry, use it regardless of rsID match.
    if len(entries) == 1 and "-" in entries[0]:
        allele = entries[0].rsplit("-", 1)[1].strip()
        if allele and allele != "?":
            return allele

    return None


def classify_or_beta(raw_value: str, ci_text: str):
    """Decide whether the OR/BETA column value is an odds ratio or a beta coefficient.

    Returns (odds_ratio, beta) — one will be the parsed float, the other None.

    Heuristic: if the value is in a plausible odds-ratio range [0.5, 20] it is
    stored as odds_ratio; otherwise it is stored as beta.
    """
    val = parse_float(raw_value)
    if val is None:
        return None, None

    if 0.5 <= val <= 20.0:
        return val, None  # odds ratio
    else:
        return None, val  # beta


def import_gwas(input_path: str, db_path: str):
    """Main import routine."""

    print(f"Opening database: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=-64000")  # 64 MB cache

    # Ensure target tables exist (they should already from schema.sql, but be safe).
    conn.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            rsid TEXT PRIMARY KEY,
            chromosome TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref_allele TEXT,
            alt_allele TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS gwas (
            rsid TEXT REFERENCES variants(rsid),
            trait TEXT NOT NULL,
            p_value REAL,
            odds_ratio REAL,
            beta REAL,
            risk_allele TEXT,
            risk_allele_frequency REAL,
            pubmed_id TEXT,
            mapped_gene TEXT
        )
    """)
    conn.commit()

    print(f"Reading GWAS Catalog TSV: {input_path}")
    start_time = time.time()

    total_rows = 0
    skipped_rows = 0
    imported_associations = 0
    unique_rsids: set = set()

    variant_batch: list = []
    gwas_batch: list = []

    def flush_batches():
        """Write accumulated batches to the database inside a transaction."""
        nonlocal variant_batch, gwas_batch
        if not variant_batch and not gwas_batch:
            return
        conn.execute("BEGIN")
        if variant_batch:
            conn.executemany(
                "INSERT OR IGNORE INTO variants (rsid, chromosome, position, ref_allele, alt_allele) "
                "VALUES (?, ?, ?, NULL, NULL)",
                variant_batch,
            )
            variant_batch = []
        if gwas_batch:
            conn.executemany(
                "INSERT INTO gwas (rsid, trait, p_value, odds_ratio, beta, risk_allele, "
                "risk_allele_frequency, pubmed_id, mapped_gene) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                gwas_batch,
            )
            gwas_batch = []
        conn.commit()

    with open(input_path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")

        # Skip the header row.
        try:
            header = next(reader)
        except StopIteration:
            print("ERROR: Input file is empty.", file=sys.stderr)
            conn.close()
            sys.exit(1)

        for row in reader:
            total_rows += 1

            if total_rows % PROGRESS_INTERVAL == 0:
                elapsed = time.time() - start_time
                print(
                    f"  Progress: {total_rows:,} rows read, "
                    f"{imported_associations:,} associations imported "
                    f"({elapsed:.1f}s elapsed)"
                )

            # Ensure the row has enough columns.
            if len(row) < MIN_COLUMNS:
                skipped_rows += 1
                continue

            # --- Extract raw values ---
            snps_raw = row[COL_SNPS].strip()
            disease_trait = row[COL_DISEASE_TRAIT].strip()
            p_value_raw = row[COL_P_VALUE].strip()
            chr_id = row[COL_CHR_ID].strip()
            chr_pos_raw = row[COL_CHR_POS].strip()
            snp_risk_allele_raw = row[COL_STRONGEST_SNP_RISK_ALLELE].strip()
            risk_af_raw = row[COL_RISK_ALLELE_FREQUENCY].strip()
            or_beta_raw = row[COL_OR_BETA].strip()
            ci_text = row[COL_CI_TEXT].strip() if len(row) > COL_CI_TEXT else ""
            pubmed_id = row[COL_PUBMEDID].strip()
            mapped_gene = row[COL_MAPPED_GENE].strip() or None

            # --- Validation ---
            if not disease_trait:
                skipped_rows += 1
                continue

            p_value = parse_float(p_value_raw)
            if p_value is None:
                skipped_rows += 1
                continue

            position = parse_position(chr_pos_raw)
            if position is None:
                skipped_rows += 1
                continue

            if not chr_id:
                skipped_rows += 1
                continue

            # Parse effect size.
            odds_ratio, beta = classify_or_beta(or_beta_raw, ci_text)

            # Parse risk allele frequency.
            risk_af = parse_float(risk_af_raw)

            # --- Handle multiple rsIDs (semicolon-separated) ---
            rsid_candidates = [s.strip() for s in snps_raw.split(";")]
            valid_rsids = [r for r in rsid_candidates if RSID_RE.match(r)]

            if not valid_rsids:
                skipped_rows += 1
                continue

            for rsid in valid_rsids:
                risk_allele = parse_risk_allele(snp_risk_allele_raw, rsid)

                # Add to variant batch.
                variant_batch.append((rsid, chr_id, position))
                unique_rsids.add(rsid)

                # Add to gwas batch.
                gwas_batch.append((
                    rsid,
                    disease_trait,
                    p_value,
                    odds_ratio,
                    beta,
                    risk_allele,
                    risk_af,
                    pubmed_id or None,
                    mapped_gene,
                ))
                imported_associations += 1

            # Flush when batches are large enough.
            if len(gwas_batch) >= BATCH_SIZE:
                flush_batches()

    # Flush remaining records.
    flush_batches()

    elapsed = time.time() - start_time
    conn.close()

    print()
    print("=" * 60)
    print("GWAS Catalog import complete")
    print("=" * 60)
    print(f"  Total input rows:        {total_rows:,}")
    print(f"  Skipped rows:            {skipped_rows:,}")
    print(f"  Imported associations:   {imported_associations:,}")
    print(f"  Unique rsIDs:            {len(unique_rsids):,}")
    print(f"  Time elapsed:            {elapsed:.1f}s")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Import GWAS Catalog TSV into the GeneSight SQLite database."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the GWAS Catalog associations TSV file.",
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Path to the GeneSight SQLite database (will be created if missing).",
    )
    args = parser.parse_args()

    import_gwas(args.input, args.db)


if __name__ == "__main__":
    main()
