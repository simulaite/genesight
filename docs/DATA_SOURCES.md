# GeneSight – Datenquellen & Zugang

## Übersicht

GeneSight nutzt ausschließlich öffentlich zugängliche, wissenschaftlich kuratierte Datenbanken.
Alle Daten werden lokal vorgehalten — es werden zur Laufzeit keine externen Anfragen gestellt.

---

## 1. ClinVar (NCBI)

**Was:** Öffentliches Archiv klinisch klassifizierter menschlicher Varianten. Über 3 Millionen Varianten mit Pathogenitäts-Bewertungen (pathogenic, likely pathogenic, benign, likely benign, uncertain significance).

**Warum:** Goldstandard für klinisch relevante Varianten. Jede Variante hat einen Review-Status (0-4 Sterne), der die Evidenzqualität angibt.

**Lizenz:** Public Domain (US Government Work)

**Zugang:**
- FTP: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`
- Relevante Dateien:
  - `tab_delimited/variant_summary.txt.gz` (~100MB) — Hauptdatei mit allen Varianten
  - `vcf_GRCh38/clinvar.vcf.gz` — VCF-Format für GRCh38
  - `vcf_GRCh37/clinvar.vcf.gz` — VCF-Format für GRCh37/hg19
- API: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
- Update-Zyklus: Wöchentlich (jeden Sonntag)

**Felder die wir nutzen:**
- `RS# (dbSNP)` — rsID für Matching
- `ClinicalSignificance` — pathogenic/benign/etc.
- `ReviewStatus` — Evidenz-Level (Sterne)
- `PhenotypeList` — Assoziierte Erkrankungen
- `GeneSymbol` — Gen-Name
- `Assembly` — GRCh37 oder GRCh38

**Attribution:** "ClinVar data provided by the National Center for Biotechnology Information (NCBI), U.S. National Library of Medicine."

---

## 2. SNPedia

**Was:** Wiki-basierte Datenbank mit ~112.000 SNPs. Jeder Eintrag verknüpft Varianten mit Peer-reviewed-Studien und enthält menschenlesbare Zusammenfassungen. Unique: Magnitude-Score (0-10) der die Relevanz einer Variante einschätzt.

**Warum:** Beste Quelle für verständliche, kontextualisierte Beschreibungen. ClinVar sagt "pathogenic", SNPedia erklärt was das bedeutet.

**Lizenz:** CC-BY-NC-SA 3.0 US
- ✅ Open-Source-Nutzung: erlaubt
- ✅ Persönliche Nutzung: erlaubt
- ❌ Kommerzielle Nutzung: nur mit separater Lizenz
- ⚠️ Share-Alike: Abgeleitete Werke müssen unter gleicher Lizenz stehen

**Zugang:**
- MediaWiki API: `https://www.snpedia.com/w/api.php`
- Bulk-Export via API möglich (kein offizieller Dump)
- Rate Limiting: Mindestens 3 Sekunden zwischen Requests
- robots.txt respektieren!

**Scraping-Strategie:**
```
1. Alle SNP-Seiten listen: api.php?action=query&list=allpages&apnamespace=0&apprefix=Rs
2. Für jede Seite: api.php?action=parse&page=Rs1234567&prop=wikitext
3. Genotyp-Seiten: api.php?action=parse&page=Rs1234567(A;G)
4. Relevante Felder extrahieren: magnitude, repute, summary, genotype-specific text
```

**Existierende Tools:**
- `TheModernPromethease` (GitHub) — R-basierter Scraper
- `SNPedia-Scraper` (GitHub) — Python-Scraper mit SQLite-Output (~160MB)

**Felder die wir nutzen:**
- `rsid` — SNP-Identifier
- `magnitude` — Relevanz (0-10, höher = wichtiger)
- `repute` — "good", "bad", oder neutral
- `summary` — Menschenlesbare Zusammenfassung
- Genotyp-spezifische Beschreibungen (z.B. was AA vs AG vs GG bedeutet)

**Attribution:** "SNPedia content is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 3.0 United States License. Source: https://www.snpedia.com"

---

## 3. GWAS Catalog (NHGRI-EBI)

**Was:** Kuratierter Katalog aller veröffentlichten Genome-Wide Association Studies. Verknüpft SNPs mit Traits/Erkrankungen inkl. Effektstärke (Odds Ratio, Beta) und p-Wert.

**Warum:** Einzige umfassende Quelle für polygene Assoziationen. Notwendig für Tier-2-Risikoscores.

**Lizenz:** Open Access (EMBL-EBI, öffentlich finanziert)

**Zugang:**
- REST API v2: `https://www.ebi.ac.uk/gwas/rest/api/v2/`
- Bulk-Download: `https://www.ebi.ac.uk/gwas/api/search/downloads/full`
- FTP: `ftp://ftp.ebi.ac.uk/pub/databases/gwas/`
- Relevante Datei: `gwas-catalog-associations_ontology-annotated.tsv` (~50MB)
- Update-Zyklus: Wöchentlich

**Felder die wir nutzen:**
- `SNPS` — rsID(s)
- `DISEASE/TRAIT` — Assoziierter Trait
- `OR or BETA` — Effektstärke
- `P-VALUE` — Statistische Signifikanz
- `RISK ALLELE FREQUENCY` — Häufigkeit des Risiko-Allels
- `MAPPED_GENE` — Zugeordnetes Gen
- `STUDY` — PubMed-ID der Originalstudie

**Attribution:** "GWAS Catalog data provided by NHGRI-EBI GWAS Catalog. Buniello A, et al. Nucleic Acids Research, 2019."

---

## 4. dbSNP (NCBI)

**Was:** Referenzdatenbank für alle bekannten Single Nucleotide Polymorphisms. Jeder SNP hat eine rs-Nummer, die als universeller Identifier dient.

**Warum:** Liefert Allelfrequenzen (wie häufig ist meine Variante in verschiedenen Populationen) und ist der Schlüssel zum Verknüpfen aller anderen Datenbanken.

**Lizenz:** Public Domain (US Government Work)

**Zugang:**
- FTP: `https://ftp.ncbi.nih.gov/snp/`
- Vollständiger Dump: ~15GB (wir brauchen nur ein Subset)
- Relevante Dateien:
  - `organisms/human_9606/VCF/` — VCF-Files mit Allelfrequenzen
  - Alternativ: gnomAD für bessere Frequenzdaten
- API: `https://api.ncbi.nlm.nih.gov/variation/v0/`

**Hinweis:** Für Allelfrequenzen ist gnomAD die bessere Quelle. dbSNP primär für rs-Nummer-Lookups und als Referenz.

**Attribution:** "dbSNP data provided by the National Center for Biotechnology Information (NCBI)."

---

## 5. gnomAD (Broad Institute)

**Was:** Genome Aggregation Database — Allelfrequenzen aus >250.000 Exomen und >76.000 Genomen, aufgeschlüsselt nach Population (European, African, East Asian, South Asian, etc.).

**Warum:** Beantwortet die Frage "Wie selten ist meine Variante?" — entscheidend für die Einschätzung klinischer Relevanz.

**Lizenz:** Open Access (ODC Open Database License für Daten)

**Zugang:**
- Download: `https://gnomad.broadinstitute.org/downloads`
- Vollständig: Multi-GB
- Für unser Tool: Nur Sites-VCF mit Frequenzen (~1-2GB für Exome)
- API: GraphQL unter `https://gnomad.broadinstitute.org/api`

**Felder die wir nutzen:**
- `rsid` / Chromosom+Position
- `AF` — Gesamte Allelfrequenz
- `AF_popmax` — Höchste Frequenz in irgendeiner Population
- Populations-spezifische Frequenzen (afr, amr, asj, eas, fin, nfe, sas)

**Attribution:** "gnomAD data provided by the Genome Aggregation Database (gnomAD), Broad Institute."

---

## 6. PharmGKB

**Was:** Curated Knowledge Base für Pharmakogenetik — welche Gene beeinflussen die Wirkung welcher Medikamente.

**Warum:** Pharmakogenetik ist einer der zuverlässigsten (Tier 1) Anwendungsbereiche der DNA-Analyse.

**Lizenz:** CC-BY-SA 4.0 (akademisch/nicht-kommerziell frei, kommerziell: Lizenz erforderlich)

**Zugang:**
- Download: `https://www.pharmgkb.org/downloads`
- Relevante Dateien:
  - `clinical_annotations.tsv` — Klinische Annotationen
  - `var_drug_ann.tsv` — Varianten-Medikamenten-Assoziationen
  - `clinical_ann_alleles.tsv` — Allel-spezifische Informationen
- Registrierung erforderlich für Bulk-Download
- API: `https://api.pharmgkb.org/`

**Felder die wir nutzen:**
- Variante (rsID)
- Medikament
- Phänotyp-Kategorie (z.B. "Poor Metabolizer", "Ultrarapid Metabolizer")
- Evidenz-Level (1A, 1B, 2A, 2B, 3, 4)
- Klinische Empfehlung

**Attribution:** "PharmGKB data © PharmGKB, licensed under CC-BY-SA 4.0. M. Whirl-Carrillo et al. Clinical Pharmacology & Therapeutics (2012)."

---

## Lokale Datenbank-Strategie

### Einheitliches SQLite-Schema

Alle Quellen werden in eine einzelne SQLite-Datei importiert (`genesight.db`):

```sql
-- Kern-Tabelle: Alle bekannten Varianten
CREATE TABLE variants (
    rsid TEXT PRIMARY KEY,          -- rs-Nummer (z.B. "rs1234567")
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref_allele TEXT,
    alt_allele TEXT
);

-- ClinVar-Annotationen
CREATE TABLE clinvar (
    rsid TEXT REFERENCES variants(rsid),
    clinical_significance TEXT,     -- pathogenic, benign, etc.
    review_status INTEGER,          -- 0-4 Sterne
    conditions TEXT,                -- JSON-Array von Erkrankungen
    gene_symbol TEXT,
    last_updated DATE
);

-- SNPedia-Annotationen
CREATE TABLE snpedia (
    rsid TEXT REFERENCES variants(rsid),
    magnitude REAL,                 -- 0-10
    repute TEXT,                    -- good, bad, null
    summary TEXT,                   -- Menschenlesbare Zusammenfassung
    genotype_descriptions TEXT      -- JSON: {"AA": "...", "AG": "...", "GG": "..."}
);

-- GWAS-Assoziationen (1:N — ein SNP kann mehrere Traits haben)
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

-- Allelfrequenzen (gnomAD/dbSNP)
CREATE TABLE frequencies (
    rsid TEXT REFERENCES variants(rsid),
    af_total REAL,                  -- Gesamt-Allelfrequenz
    af_afr REAL,                    -- African
    af_amr REAL,                    -- American
    af_eas REAL,                    -- East Asian
    af_eur REAL,                    -- European (non-Finnish)
    af_sas REAL,                    -- South Asian
    source TEXT                     -- "gnomad" oder "dbsnp"
);

-- Pharmakogenetik
CREATE TABLE pharmacogenomics (
    rsid TEXT REFERENCES variants(rsid),
    drug TEXT NOT NULL,
    phenotype_category TEXT,        -- Poor/Intermediate/Normal/Rapid/Ultrarapid Metabolizer
    evidence_level TEXT,            -- 1A, 1B, 2A, 2B, 3, 4
    clinical_recommendation TEXT,
    gene_symbol TEXT
);

-- Indizes für schnelle Lookups
CREATE INDEX idx_clinvar_rsid ON clinvar(rsid);
CREATE INDEX idx_snpedia_rsid ON snpedia(rsid);
CREATE INDEX idx_gwas_rsid ON gwas(rsid);
CREATE INDEX idx_freq_rsid ON frequencies(rsid);
CREATE INDEX idx_pharma_rsid ON pharmacogenomics(rsid);
CREATE INDEX idx_variants_chr_pos ON variants(chromosome, position);
```

### Erwartete Größen

| Tabelle | Zeilen (ca.) | Größe (ca.) |
|---------|-------------|-------------|
| variants | ~15M (Subset) | ~300MB |
| clinvar | ~3M | ~100MB |
| snpedia | ~112K | ~30MB |
| gwas | ~500K | ~20MB |
| frequencies | ~10M (Subset) | ~200MB |
| pharmacogenomics | ~50K | ~5MB |
| **Gesamt** | | **~500MB-1GB** |

---

## Update-Strategie

- ClinVar: Wöchentlicher FTP-Download, Delta-Import
- SNPedia: Monatlicher Re-Scrape (MediaWiki Recent Changes API für Deltas)
- GWAS Catalog: Monatlicher Re-Download
- gnomAD: Stabile Releases, Update bei neuer Version
- PharmGKB: Quartalsweise Check

Das CLI-Tool soll einen `genesight update` Command haben, der alle Datenbanken aktualisiert.
