#!/usr/bin/env bash
# Fetch dbSNP common variants subset
# Output: data/raw/dbsnp/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../raw/dbsnp"
mkdir -p "$DATA_DIR"

echo "Downloading dbSNP common variants (GRCh38)..."
# Download the common SNPs VCF (much smaller than full dbSNP)
curl -L -o "$DATA_DIR/common_all.vcf.gz" \
    "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"

echo "dbSNP data downloaded to $DATA_DIR"
echo "File size: $(du -h "$DATA_DIR/common_all.vcf.gz" | cut -f1)"
