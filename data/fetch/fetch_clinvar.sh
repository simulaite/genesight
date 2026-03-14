#!/usr/bin/env bash
# Fetch ClinVar data from NCBI FTP
# Output: data/raw/clinvar/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../raw/clinvar"
mkdir -p "$DATA_DIR"

echo "Downloading ClinVar variant summary..."
curl -L -o "$DATA_DIR/variant_summary.txt.gz" \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

echo "Extracting..."
gunzip -f "$DATA_DIR/variant_summary.txt.gz"

echo "ClinVar data downloaded to $DATA_DIR"
echo "File size: $(du -h "$DATA_DIR/variant_summary.txt" | cut -f1)"
