-- GeneSight SQLite Database Schema
-- All public genome databases are imported into a single local SQLite file.
-- SNPedia is stored in a separate optional database (snpedia.db) due to CC-BY-NC-SA 3.0 license.

-- ============================================================
-- genesight.db — Main database
-- ============================================================

-- Schema version tracking
CREATE TABLE IF NOT EXISTS schema_version (
    version INTEGER PRIMARY KEY,
    applied_at TEXT NOT NULL DEFAULT (datetime('now'))
);

INSERT INTO schema_version (version) VALUES (1);

-- Key-value metadata (assembly version, build date, etc.)
CREATE TABLE IF NOT EXISTS db_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

-- Core table: All known variants (reference data)
CREATE TABLE IF NOT EXISTS variants (
    rsid TEXT PRIMARY KEY,              -- rs-number (e.g., "rs1234567")
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref_allele TEXT,
    alt_allele TEXT
);
CREATE INDEX IF NOT EXISTS idx_variants_chr_pos ON variants(chromosome, position);

-- ClinVar: Clinically classified variants
CREATE TABLE IF NOT EXISTS clinvar (
    rsid TEXT REFERENCES variants(rsid),
    clinical_significance TEXT,         -- pathogenic, benign, etc.
    review_status INTEGER,              -- 0-4 stars
    conditions TEXT,                    -- JSON array of conditions
    gene_symbol TEXT,
    last_updated DATE
);
CREATE INDEX IF NOT EXISTS idx_clinvar_rsid ON clinvar(rsid);

-- GWAS Catalog: Genome-wide association study results (1:N per rsid)
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
);
CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON gwas(rsid);

-- Allele frequencies (gnomAD/dbSNP)
CREATE TABLE IF NOT EXISTS frequencies (
    rsid TEXT REFERENCES variants(rsid),
    af_total REAL,                      -- Overall allele frequency
    af_afr REAL,                        -- African
    af_amr REAL,                        -- American
    af_eas REAL,                        -- East Asian
    af_eur REAL,                        -- European (non-Finnish)
    af_sas REAL,                        -- South Asian
    source TEXT                         -- "gnomad" or "dbsnp"
);
CREATE INDEX IF NOT EXISTS idx_freq_rsid ON frequencies(rsid);

-- Pharmacogenomics (PharmGKB)
CREATE TABLE IF NOT EXISTS pharmacogenomics (
    rsid TEXT REFERENCES variants(rsid),
    drug TEXT NOT NULL,
    phenotype_category TEXT,            -- Poor/Intermediate/Normal/Rapid/Ultrarapid Metabolizer
    evidence_level TEXT,                -- 1A, 1B, 2A, 2B, 3, 4
    clinical_recommendation TEXT,
    gene_symbol TEXT
);
CREATE INDEX IF NOT EXISTS idx_pharma_rsid ON pharmacogenomics(rsid);

-- PGx star allele definitions (CPIC-style)
CREATE TABLE IF NOT EXISTS pgx_allele_definitions (
    gene TEXT NOT NULL,                 -- Gene symbol (e.g., "CYP2D6")
    allele_name TEXT NOT NULL,          -- Star allele name (e.g., "*4")
    rsid TEXT REFERENCES variants(rsid),
    alt_allele TEXT NOT NULL,           -- Alternate allele defining this star allele
    function TEXT NOT NULL,             -- "No Function", "Decreased Function", "Increased Function", etc.
    activity_score REAL NOT NULL        -- Numeric activity score for phenotype calculation
);
CREATE INDEX IF NOT EXISTS idx_pgx_allele_gene ON pgx_allele_definitions(gene);
CREATE INDEX IF NOT EXISTS idx_pgx_allele_rsid ON pgx_allele_definitions(rsid);

-- PGx diplotype-to-phenotype mapping
CREATE TABLE IF NOT EXISTS pgx_diplotype_phenotypes (
    gene TEXT NOT NULL,                 -- Gene symbol
    diplotype TEXT NOT NULL,            -- e.g., "*1/*2", "*4/*4"
    phenotype TEXT NOT NULL,            -- e.g., "Poor Metabolizer"
    activity_score REAL NOT NULL        -- Combined activity score for this diplotype
);
CREATE INDEX IF NOT EXISTS idx_pgx_diplo_gene ON pgx_diplotype_phenotypes(gene);

-- PGx drug recommendations (CPIC guideline-style)
CREATE TABLE IF NOT EXISTS pgx_drug_recommendations (
    gene TEXT NOT NULL,                 -- Gene symbol
    drug TEXT NOT NULL,                 -- Drug name
    phenotype TEXT NOT NULL,            -- Required phenotype for this recommendation
    recommendation TEXT NOT NULL,       -- Clinical recommendation text
    evidence_level TEXT NOT NULL        -- CPIC evidence level (1A, 1B, 2A, etc.)
);
CREATE INDEX IF NOT EXISTS idx_pgx_drug_gene ON pgx_drug_recommendations(gene);


-- ============================================================
-- snpedia.db — Separate optional database (CC-BY-NC-SA 3.0)
-- Created by: genesight fetch --snpedia
-- ============================================================

-- To be applied to snpedia.db separately:
--
-- CREATE TABLE IF NOT EXISTS snpedia (
--     rsid TEXT PRIMARY KEY,
--     magnitude REAL,                   -- 0-10 importance score
--     repute TEXT,                      -- "good", "bad", or null
--     summary TEXT,                     -- human-readable summary
--     genotype_descriptions TEXT        -- JSON: {"AA": "...", "AG": "...", "GG": "..."}
-- );
-- CREATE INDEX IF NOT EXISTS idx_snpedia_rsid ON snpedia(rsid);
