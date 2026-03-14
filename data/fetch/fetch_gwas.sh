#!/usr/bin/env bash
# Fetch GWAS Catalog data from EMBL-EBI
# Output: data/raw/gwas/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../raw/gwas"
mkdir -p "$DATA_DIR"

echo "Downloading GWAS Catalog associations..."
curl -L -o "$DATA_DIR/gwas_catalog_associations.tsv" \
    "https://www.ebi.ac.uk/gwas/api/search/downloads/full"

echo "GWAS Catalog data downloaded to $DATA_DIR"
echo "File size: $(du -h "$DATA_DIR/gwas_catalog_associations.tsv" | cut -f1)"
