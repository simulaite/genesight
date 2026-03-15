# GeneSight — Data Sources & Access

## Overview

GeneSight exclusively uses publicly available, scientifically curated databases.
All data is stored locally — no external requests are made at runtime.

---

## 1. ClinVar (NCBI)

**What:** Public archive of clinically classified human variants. Over 3 million variants with pathogenicity assessments (pathogenic, likely pathogenic, benign, likely benign, uncertain significance).

**Why:** Gold standard for clinically relevant variants. Each variant has a review status (0-4 stars) indicating the quality of evidence.

**License:** Public Domain (US Government Work)

**Access:**
- FTP: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`
- Relevant files:
  - `tab_delimited/variant_summary.txt.gz` (~100MB) — Main file with all variants
  - `vcf_GRCh38/clinvar.vcf.gz` — VCF format for GRCh38
  - `vcf_GRCh37/clinvar.vcf.gz` — VCF format for GRCh37/hg19
- API: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
- Update cycle: Weekly (every Sunday)

**Fields we use:**
- `RS# (dbSNP)` — rsID for matching
- `ClinicalSignificance` — pathogenic/benign/etc.
- `ReviewStatus` — Evidence level (stars)
- `PhenotypeList` — Associated conditions
- `GeneSymbol` — Gene name
- `Assembly` — GRCh37 or GRCh38

**Attribution:** "ClinVar data provided by the National Center for Biotechnology Information (NCBI), U.S. National Library of Medicine."

---

## 2. SNPedia

**What:** Wiki-based database with ~112,000 SNPs. Each entry links variants to peer-reviewed studies and includes human-readable summaries. Unique feature: Magnitude score (0-10) that estimates the relevance of a variant.

**Why:** Best source for understandable, contextualized descriptions. ClinVar says "pathogenic", SNPedia explains what that means.

**License:** CC-BY-NC-SA 3.0 US
- Open-source use: permitted
- Personal use: permitted
- Commercial use: only with a separate license
- Share-Alike: derivative works must be released under the same license

**Access:**
- MediaWiki API: `https://www.snpedia.com/w/api.php`
- Bulk export via API possible (no official dump available)
- Rate limiting: At least 3 seconds between requests
- Respect robots.txt!

**Scraping strategy:**
```
1. List all SNP pages: api.php?action=query&list=allpages&apnamespace=0&apprefix=Rs
2. For each page: api.php?action=parse&page=Rs1234567&prop=wikitext
3. Genotype pages: api.php?action=parse&page=Rs1234567(A;G)
4. Extract relevant fields: magnitude, repute, summary, genotype-specific text
```

**Existing tools:**
- `TheModernPromethease` (GitHub) — R-based scraper
- `SNPedia-Scraper` (GitHub) — Python scraper with SQLite output (~160MB)

**Fields we use:**
- `rsid` — SNP identifier
- `magnitude` — Relevance (0-10, higher = more important)
- `repute` — "good", "bad", or neutral
- `summary` — Human-readable summary
- Genotype-specific descriptions (e.g., what AA vs AG vs GG means)

**Attribution:** "SNPedia content is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 3.0 United States License. Source: https://www.snpedia.com"

---

## 3. GWAS Catalog (NHGRI-EBI)

**What:** Curated catalog of all published Genome-Wide Association Studies. Links SNPs to traits/diseases including effect size (Odds Ratio, Beta) and p-value.

**Why:** Only comprehensive source for polygenic associations. Required for Tier 2 risk scores.

**License:** Open Access (EMBL-EBI, publicly funded)

**Access:**
- REST API v2: `https://www.ebi.ac.uk/gwas/rest/api/v2/`
- Bulk download: `https://www.ebi.ac.uk/gwas/api/search/downloads/full`
- FTP: `ftp://ftp.ebi.ac.uk/pub/databases/gwas/`
- Relevant file: `gwas-catalog-associations_ontology-annotated.tsv` (~50MB)
- Update cycle: Weekly

**Fields we use:**
- `SNPS` — rsID(s)
- `DISEASE/TRAIT` — Associated trait
- `OR or BETA` — Effect size
- `P-VALUE` — Statistical significance
- `RISK ALLELE FREQUENCY` — Frequency of the risk allele
- `MAPPED_GENE` — Mapped gene
- `STUDY` — PubMed ID of the original study

**Attribution:** "GWAS Catalog data provided by NHGRI-EBI GWAS Catalog. Buniello A, et al. Nucleic Acids Research, 2019."

---

## 4. dbSNP (NCBI)

**What:** Reference database for all known Single Nucleotide Polymorphisms. Each SNP has an rs number that serves as a universal identifier.

**Why:** Provides allele frequencies (how common is my variant across different populations) and is the key to linking all other databases.

**License:** Public Domain (US Government Work)

**Access:**
- FTP: `https://ftp.ncbi.nih.gov/snp/`
- Full dump: ~15GB (we only need a subset)
- Relevant files:
  - `organisms/human_9606/VCF/` — VCF files with allele frequencies
  - Alternative: gnomAD for better frequency data
- API: `https://api.ncbi.nlm.nih.gov/variation/v0/`

**Note:** For allele frequencies, gnomAD is the better source. dbSNP is primarily used for rs number lookups and as a reference.

**Attribution:** "dbSNP data provided by the National Center for Biotechnology Information (NCBI)."

---

## 5. gnomAD (Broad Institute)

**What:** Genome Aggregation Database — allele frequencies from >250,000 exomes and >76,000 genomes, broken down by population (European, African, East Asian, South Asian, etc.).

**Why:** Answers the question "How rare is my variant?" — critical for assessing clinical relevance.

**License:** Open Access (ODC Open Database License for data)

**Access:**
- Download: `https://gnomad.broadinstitute.org/downloads`
- Full dataset: Multi-GB
- For our tool: Only sites VCF with frequencies (~1-2GB for exomes)
- API: GraphQL at `https://gnomad.broadinstitute.org/api`

**Fields we use:**
- `rsid` / Chromosome+Position
- `AF` — Overall allele frequency
- `AF_popmax` — Highest frequency in any population
- Population-specific frequencies (afr, amr, asj, eas, fin, nfe, sas)

**Attribution:** "gnomAD data provided by the Genome Aggregation Database (gnomAD), Broad Institute."

---

## 6. PharmGKB

**What:** Curated knowledge base for pharmacogenetics — which genes influence the effect of which drugs.

**Why:** Pharmacogenetics is one of the most reliable (Tier 1) application areas of DNA analysis.

**License:** CC-BY-SA 4.0 (free for academic/non-commercial use, commercial: license required)

**Access:**
- Download: `https://www.pharmgkb.org/downloads`
- Relevant files:
  - `clinical_annotations.tsv` — Clinical annotations
  - `var_drug_ann.tsv` — Variant-drug associations
  - `clinical_ann_alleles.tsv` — Allele-specific information
- Registration required for bulk download
- API: `https://api.pharmgkb.org/`

**Fields we use:**
- Variant (rsID)
- Drug
- Phenotype category (e.g., "Poor Metabolizer", "Ultrarapid Metabolizer")
- Evidence level (1A, 1B, 2A, 2B, 3, 4)
- Clinical recommendation

**Attribution:** "PharmGKB data © PharmGKB, licensed under CC-BY-SA 4.0. M. Whirl-Carrillo et al. Clinical Pharmacology & Therapeutics (2012)."

---

## Local Database Strategy

### Unified SQLite Schema

All sources are imported into a single SQLite file (`genesight.db`):

```sql
-- Core table: All known variants
CREATE TABLE variants (
    rsid TEXT PRIMARY KEY,          -- rs number (e.g., "rs1234567")
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref_allele TEXT,
    alt_allele TEXT
);

-- ClinVar annotations
CREATE TABLE clinvar (
    rsid TEXT REFERENCES variants(rsid),
    clinical_significance TEXT,     -- pathogenic, benign, etc.
    review_status INTEGER,          -- 0-4 stars
    conditions TEXT,                -- JSON array of conditions
    gene_symbol TEXT,
    last_updated DATE
);

-- SNPedia annotations
CREATE TABLE snpedia (
    rsid TEXT REFERENCES variants(rsid),
    magnitude REAL,                 -- 0-10
    repute TEXT,                    -- good, bad, null
    summary TEXT,                   -- Human-readable summary
    genotype_descriptions TEXT      -- JSON: {"AA": "...", "AG": "...", "GG": "..."}
);

-- GWAS associations (1:N — one SNP can have multiple traits)
CREATE TABLE gwas (
    rsid TEXT REFERENCES variants(rsid),
    trait TEXT NOT NULL,
    p_value REAL,
    odds_ratio REAL,
    beta REAL,
    risk_allele TEXT,
    risk_allele_frequency REAL,
    pubmed_id TEXT,
    mapped_gene TEXT
);

-- Allele frequencies (gnomAD/dbSNP)
CREATE TABLE frequencies (
    rsid TEXT REFERENCES variants(rsid),
    af_total REAL,                  -- Overall allele frequency
    af_afr REAL,                    -- African
    af_amr REAL,                    -- American
    af_eas REAL,                    -- East Asian
    af_eur REAL,                    -- European (non-Finnish)
    af_sas REAL,                    -- South Asian
    source TEXT                     -- "gnomad" or "dbsnp"
);

-- Pharmacogenomics
CREATE TABLE pharmacogenomics (
    rsid TEXT REFERENCES variants(rsid),
    drug TEXT NOT NULL,
    phenotype_category TEXT,        -- Poor/Intermediate/Normal/Rapid/Ultrarapid Metabolizer
    evidence_level TEXT,            -- 1A, 1B, 2A, 2B, 3, 4
    clinical_recommendation TEXT,
    gene_symbol TEXT
);

-- Indexes for fast lookups
CREATE INDEX idx_clinvar_rsid ON clinvar(rsid);
CREATE INDEX idx_snpedia_rsid ON snpedia(rsid);
CREATE INDEX idx_gwas_rsid ON gwas(rsid);
CREATE INDEX idx_freq_rsid ON frequencies(rsid);
CREATE INDEX idx_pharma_rsid ON pharmacogenomics(rsid);
CREATE INDEX idx_variants_chr_pos ON variants(chromosome, position);
```

### Expected Sizes

| Table | Rows (approx.) | Size (approx.) |
|-------|----------------|----------------|
| variants | ~15M (subset) | ~300MB |
| clinvar | ~3M | ~100MB |
| snpedia | ~112K | ~30MB |
| gwas | ~500K | ~20MB |
| frequencies | ~10M (subset) | ~200MB |
| pharmacogenomics | ~50K | ~5MB |
| **Total** | | **~500MB-1GB** |

---

## Update Strategy

- ClinVar: Weekly FTP download, delta import
- SNPedia: Monthly re-scrape (MediaWiki Recent Changes API for deltas)
- GWAS Catalog: Monthly re-download
- gnomAD: Stable releases, update on new version
- PharmGKB: Quarterly check

The CLI tool should have a `genesight update` command that updates all databases.
