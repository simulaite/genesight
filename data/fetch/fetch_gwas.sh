#!/usr/bin/env bash
# Fetch GWAS Catalog data from EMBL-EBI
# Output: data/raw/gwas/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../raw/gwas"
mkdir -p "$DATA_DIR"

echo "Downloading GWAS Catalog associations (full)..."
curl -L -o "$DATA_DIR/gwas-catalog-associations-full.zip" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"

echo "Extracting..."
unzip -o "$DATA_DIR/gwas-catalog-associations-full.zip" -d "$DATA_DIR"

echo "GWAS Catalog data downloaded to $DATA_DIR"
ls -lh "$DATA_DIR"
