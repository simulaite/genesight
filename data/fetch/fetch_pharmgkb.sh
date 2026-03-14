#!/usr/bin/env bash
# Fetch PharmGKB clinical annotations
# Note: Requires accepting PharmGKB license terms at https://www.pharmgkb.org/downloads
# Output: data/raw/pharmgkb/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../raw/pharmgkb"
mkdir -p "$DATA_DIR"

echo "PharmGKB requires manual download due to license agreement."
echo ""
echo "Please:"
echo "  1. Visit https://www.pharmgkb.org/downloads"
echo "  2. Accept the license terms"
echo "  3. Download 'Clinical Annotations' and 'Clinical Annotation Alleles'"
echo "  4. Place the ZIP files in: $DATA_DIR"
echo ""
echo "After downloading, run this script again to extract."

if [ -f "$DATA_DIR/clinical_annotations.zip" ]; then
    echo "Found clinical_annotations.zip, extracting..."
    unzip -o "$DATA_DIR/clinical_annotations.zip" -d "$DATA_DIR"
    echo "PharmGKB data extracted to $DATA_DIR"
fi
