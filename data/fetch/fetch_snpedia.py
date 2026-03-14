#!/usr/bin/env python3
"""Fetch SNPedia data via MediaWiki API.

Respects rate limits: 3-second delay between requests.
Output: data/raw/snpedia/
"""

import json
import os
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / ".." / "raw" / "snpedia"
API_URL = "https://bots.snpedia.com/api.php"
DELAY = 3  # seconds between requests

def fetch_snp_list():
    """Fetch the list of all SNP pages from SNPedia."""
    snps = []
    params = {
        "action": "query",
        "list": "categorymembers",
        "cmtitle": "Category:Is_a_snp",
        "cmlimit": "500",
        "format": "json",
    }

    while True:
        url = f"{API_URL}?{urlencode(params)}"
        req = Request(url, headers={"User-Agent": "GeneSight/0.1 (Open Source DNA Tool)"})
        with urlopen(req) as resp:
            data = json.loads(resp.read())

        members = data.get("query", {}).get("categorymembers", [])
        snps.extend(m["title"] for m in members)
        print(f"  Fetched {len(snps)} SNPs so far...")

        if "continue" in data:
            params["cmcontinue"] = data["continue"]["cmcontinue"]
            time.sleep(DELAY)
        else:
            break

    return snps


def fetch_snp_details(title: str) -> dict | None:
    """Fetch details for a single SNP page."""
    params = {
        "action": "parse",
        "page": title,
        "prop": "wikitext",
        "format": "json",
    }
    url = f"{API_URL}?{urlencode(params)}"
    req = Request(url, headers={"User-Agent": "GeneSight/0.1 (Open Source DNA Tool)"})
    try:
        with urlopen(req) as resp:
            data = json.loads(resp.read())
        return data.get("parse", {})
    except Exception as e:
        print(f"  Error fetching {title}: {e}")
        return None


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("Fetching SNP list from SNPedia...")
    snps = fetch_snp_list()
    print(f"Found {len(snps)} SNPs")

    # Save the list
    list_file = DATA_DIR / "snp_list.json"
    with open(list_file, "w") as f:
        json.dump(snps, f)

    # Fetch details for each SNP
    details_dir = DATA_DIR / "pages"
    details_dir.mkdir(exist_ok=True)

    for i, snp in enumerate(snps):
        out_file = details_dir / f"{snp}.json"
        if out_file.exists():
            continue

        print(f"[{i+1}/{len(snps)}] Fetching {snp}...")
        details = fetch_snp_details(snp)
        if details:
            with open(out_file, "w") as f:
                json.dump(details, f)

        time.sleep(DELAY)

    print(f"SNPedia data saved to {DATA_DIR}")


if __name__ == "__main__":
    main()
