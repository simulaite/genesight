# GeneSight – Open-Source DNA-Analyse-Tool

## Projekt-Identität

**Name:** GeneSight (Arbeitstitel)
**Sprache:** Deutsch & Englisch (Code und API auf Englisch, Dokumentation bilingual)
**Lizenz:** GPL-3.0-or-later (kompatibel mit CC-BY-NC-SA 3.0 von SNPedia)
**Ziel:** Ein Privacy-first CLI- und Desktop-Tool, das persönliche DNA-Rohdaten (23andMe, AncestryDNA, VCF) gegen öffentliche Genomdatenbanken annotiert und verständliche Reports generiert – ohne dass Daten jemals den Rechner des Nutzers verlassen.

---

## Architektur-Überblick

```
genesight/
├── crates/
│   ├── genesight-core/       # Library Crate: Parser, Annotator, Scorer, Report-Engine
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── parser/       # DNA-Datei-Parser (23andMe, AncestryDNA, VCF)
│   │   │   │   ├── mod.rs
│   │   │   │   ├── twentythreeandme.rs
│   │   │   │   ├── ancestry.rs
│   │   │   │   └── vcf.rs
│   │   │   ├── db/           # Datenbank-Adapter (lokale SQLite)
│   │   │   │   ├── mod.rs
│   │   │   │   ├── clinvar.rs
│   │   │   │   ├── snpedia.rs
│   │   │   │   ├── gwas.rs
│   │   │   │   ├── dbsnp.rs
│   │   │   │   └── pharmgkb.rs
│   │   │   ├── annotator/    # Varianten-Annotation gegen Datenbanken
│   │   │   │   ├── mod.rs
│   │   │   │   ├── clinical.rs    # ClinVar pathogenicity
│   │   │   │   ├── frequency.rs   # gnomAD/dbSNP Allelfrequenzen
│   │   │   │   ├── pharmacogenomics.rs
│   │   │   │   └── traits.rs      # SNPedia traits & magnitude
│   │   │   ├── scorer/       # Risiko-Scoring & Confidence-Tiers
│   │   │   │   ├── mod.rs
│   │   │   │   ├── monogenic.rs   # Einzelgen-Erkrankungen (Tier 1: >95%)
│   │   │   │   ├── pharmaco.rs    # Pharmakogenetik (Tier 1: >95%)
│   │   │   │   ├── polygenic.rs   # Polygene Risikoscores (Tier 2: 60-85%)
│   │   │   │   └── traits.rs      # Merkmale & Lifestyle (Tier 2-3)
│   │   │   ├── report/       # Report-Generierung
│   │   │   │   ├── mod.rs
│   │   │   │   ├── markdown.rs
│   │   │   │   ├── json.rs
│   │   │   │   └── html.rs
│   │   │   └── models/       # Shared Types & Structs
│   │   │       ├── mod.rs
│   │   │       ├── variant.rs
│   │   │       ├── annotation.rs
│   │   │       ├── confidence.rs  # ConfidenceTier enum
│   │   │       └── report.rs
│   │   └── Cargo.toml
│   ├── genesight-cli/        # CLI-Tool (clap)
│   │   ├── src/
│   │   │   └── main.rs
│   │   └── Cargo.toml
│   ├── genesight-server/     # Axum API (optional, für Web-Version)
│   │   ├── src/
│   │   │   └── main.rs
│   │   └── Cargo.toml
│   └── genesight-desktop/    # Tauri App (Phase 2)
│       └── ...
├── data/
│   ├── fetch/                # Scripts zum Herunterladen der Datenbanken
│   │   ├── fetch_clinvar.sh
│   │   ├── fetch_snpedia.py
│   │   ├── fetch_gwas.sh
│   │   ├── fetch_dbsnp.sh
│   │   └── fetch_pharmgkb.sh
│   ├── import/               # Scripts zum Importieren in SQLite
│   │   ├── import_clinvar.rs (oder .py)
│   │   ├── import_snpedia.rs
│   │   └── import_gwas.rs
│   └── schema/               # SQLite Schema-Definitionen
│       └── schema.sql
├── tests/
│   ├── fixtures/             # Test-DNA-Dateien (synthetisch!)
│   │   ├── sample_23andme.txt
│   │   ├── sample_ancestry.txt
│   │   └── sample.vcf
│   └── integration/
├── docs/
│   ├── ARCHITECTURE.md
│   ├── DATA_SOURCES.md
│   ├── LICENSES.md
│   ├── CONFIDENCE_TIERS.md
│   └── CONTRIBUTING.md
├── Cargo.toml               # Workspace
├── CLAUDE.md                 # Diese Datei (Claude Code Kontext)
├── LICENSE                   # GPL-3.0
└── README.md
```

---

## Datenquellen & Lizenzen

### Primäre Datenbanken

| Datenbank | Inhalt | Lizenz | Zugang | Größe (ca.) |
|-----------|--------|--------|--------|-------------|
| **ClinVar** | Klinisch klassifizierte Varianten (pathogenic/benign), >3M Varianten | Public Domain (US Gov) | FTP: `ftp.ncbi.nlm.nih.gov/pub/clinvar/` + REST API | ~100MB (TSV) |
| **SNPedia** | Wiki mit ~112K SNPs, Magnitude-Scores, menschenlesbare Zusammenfassungen | CC-BY-NC-SA 3.0 | MediaWiki API: `snpedia.com/w/api.php` | ~160MB (SQLite dump) |
| **GWAS Catalog** | Genom-weite Assoziationsstudien, polygene Traits | Open Access (EMBL-EBI) | REST API v2: `ebi.ac.uk/gwas/rest/api/v2/` + FTP | ~50MB |
| **dbSNP** | Referenz-SNP-Datenbank (rs-Nummern, Allelfrequenzen) | Public Domain (US Gov) | FTP: `ftp.ncbi.nih.gov/snp/` | ~15GB (vollständig), Subset ~500MB |
| **gnomAD** | Allelfrequenzen aus >250K Genomen | Open Access | Download: `gnomad.broadinstitute.org` | Multi-GB, Subset ~1GB |
| **PharmGKB** | Pharmakogenetik (Medikamenten-Gen-Interaktionen) | CC-BY-SA 4.0 (akademisch frei) | Download + API: `pharmgkb.org` | ~50MB |

### Lizenz-Kompatibilität

- **GPL-3.0** (unser Projekt) ist kompatibel mit:
  - Public Domain (ClinVar, dbSNP) ✅
  - CC-BY-NC-SA 3.0 (SNPedia) ✅ — solange wir nicht-kommerziell bleiben oder SNPedia-Daten als separaten, optional herunterladbaren Datensatz behandeln
  - CC-BY-SA 4.0 (PharmGKB) ✅
  - Open Access (GWAS Catalog, gnomAD) ✅

- **Wichtig:** CC-BY-NC-SA 3.0 von SNPedia bedeutet:
  - ✅ Open-Source-Projekt: kein Problem
  - ✅ Persönliche/akademische Nutzung: kein Problem
  - ⚠️ Falls jemand das Projekt kommerziell forken will: SNPedia-Daten müssen entfernt oder separat lizenziert werden
  - → **Lösung:** SNPedia-Daten als optionalen Download behandeln, nicht im Repo bündeln

### Attributions-Pflichten

Jede Nutzung muss korrekt attribuieren:
- ClinVar: "ClinVar data provided by NCBI (National Center for Biotechnology Information)"
- SNPedia: "SNPedia content is licensed under CC-BY-NC-SA 3.0 by SNPedia.com"
- GWAS Catalog: "GWAS Catalog provided by NHGRI-EBI"
- PharmGKB: "PharmGKB data © PharmGKB, CC-BY-SA 4.0"

---

## Confidence-Tier-System

Alle Ergebnisse werden in drei Zuverlässigkeitsstufen eingeteilt:

### Tier 1: Zuverlässig (>95% Genauigkeit)
- **Monogenetische Erkrankungen** — Einzelne Variante ist direkt kausal (z.B. BRCA1/2, CFTR, Huntington)
- **Carrier Status** — Trägerstatus für rezessive Erkrankungen
- **Pharmakogenetik** — Medikamenten-Metabolismus (CYP2D6, CYP2C19, etc.)
- **Einfache Merkmale** — Laktosetoleranz, Ohrenschmalz-Typ, etc.
- Quelle: Primär ClinVar (review status ≥ 2 Sterne), PharmGKB (Level 1-2)

### Tier 2: Wahrscheinlich (60-85% Genauigkeit)
- **Polygene Risikoscores** — Diabetes, Herzkrankheiten, Bluthochdruck
- **Körperliche Merkmale** — Haarfarbe, Sommersprossen, Glatzenrisiko
- Quelle: GWAS Catalog, SNPedia (Magnitude ≥ 3)

### Tier 3: Spekulativ (50-65% Genauigkeit)
- **Komplexe Erkrankungen** — Depression, Schizophrenie, Autismus
- **Persönlichkeitsmerkmale** — Intelligenz, Risikobereitschaft
- **Sportliche Eignung** — ACTN3, VO2max-Prädisposition
- Quelle: GWAS Catalog (niedrige Effektstärke), SNPedia (Magnitude < 3)

**Regel:** Jedes Ergebnis MUSS ein `ConfidenceTier` zugewiesen bekommen. Der Report zeigt dies prominent an.

---

## Unterstützte Eingabeformate

### 23andMe Raw Data
```
# rsid  chromosome  position  genotype
rs4477212  1  82154  AA
rs3094315  1  752566  AG
```
- Tab-separated, Kommentarzeilen beginnen mit `#`
- Header-Zeile: `rsid  chromosome  position  genotype`
- Genotyp: 2 Buchstaben (z.B. AA, AG, CT), `--` für nicht-aufgerufen, `I` oder `D` für Indels

### AncestryDNA Raw Data
```
rsid  chromosome  position  allele1  allele2
rs4477212  1  82154  A  A
```
- Tab-separated, Kommentarzeilen beginnen mit `#`
- Allele sind getrennt in zwei Spalten

### VCF (Variant Call Format)
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
1  82154  rs4477212  G  A  .  PASS  .  GT  0/1
```
- Standard-Bioinformatik-Format
- Komplexer zu parsen, aber am vollständigsten

---

## Entwicklungs-Phasen

### Phase 1: Daten & CLI (AKTUELL)
1. **Daten-Fetching-Scripts** — Alle Datenbanken herunterladen
2. **SQLite-Import** — Daten in lokale, abfragbare Datenbank importieren
3. **DNA-Parser** — 23andMe, AncestryDNA, VCF Dateien einlesen
4. **Annotation-Engine** — Varianten gegen lokale DB matchen
5. **CLI-Tool** — `genesight analyze my_dna.txt --format markdown`
6. **Report-Generator** — Markdown/JSON/HTML Output mit Confidence-Tiers

### Phase 2: Desktop App
7. **Tauri-Integration** — GUI um den Core
8. **Auto-Update** — Datenbank-Updates im Hintergrund
9. **LLM-Integration** — Optional: Ergebnisse per LLM verständlich zusammenfassen

### Phase 3: Web & Community
10. **Axum API** — Für Web-Version (mit explizitem Privacy-Disclaimer)
11. **Community-Reports** — Anonymisierte, aggregierte Statistiken

---

## Coding-Konventionen

### Rust
- **Edition:** 2021
- **MSRV:** 1.75+
- **Error Handling:** `thiserror` für Library-Errors, `anyhow` für CLI/App
- **Async:** `tokio` (für Daten-Fetching und Server), sync für Core-Logik
- **Serialization:** `serde` + `serde_json`
- **CLI:** `clap` v4 (derive API)
- **Database:** `rusqlite` (mit bundled SQLite)
- **HTTP Client:** `reqwest` (für Daten-Fetching)
- **Testing:** Unit-Tests in jedem Modul, Integration-Tests in `tests/`

### Code-Stil
- `cargo fmt` und `cargo clippy` müssen sauber durchlaufen
- Alle öffentlichen Funktionen haben Doc-Comments
- Keine `unwrap()` in Library-Code — nur in Tests und CLI mit Kontext
- Englische Variablen- und Funktionsnamen
- Deutsche Kommentare sind OK, Doc-Comments auf Englisch

### Git
- Conventional Commits: `feat:`, `fix:`, `docs:`, `data:`, `refactor:`
- Branch-Schema: `feat/parser-23andme`, `data/clinvar-import`
- Keine DNA-Rohdaten im Repo — nur synthetische Testdaten

---

## Wichtige Regeln

1. **Keine echten DNA-Daten im Repository.** Testdaten müssen synthetisch generiert werden.
2. **Keine medizinischen Diagnosen.** Der Report ist informativ, nicht diagnostisch. Jeder Report enthält einen Disclaimer.
3. **Privacy first.** Keine Telemetrie, keine Daten-Uploads, keine Analytics. Alles lokal.
4. **Confidence-Tiers sind Pflicht.** Kein Ergebnis ohne zugewiesene Zuverlässigkeitsstufe.
5. **Attributions sind Pflicht.** Jede Datenquelle muss im Report korrekt attribuiert werden.
6. **Offline-fähig.** Nach initialem Datenbank-Download muss das Tool komplett offline funktionieren.

---

## Aktueller Fokus: Phase 1 – Daten holen & CLI

### Aufgabe 1: Daten-Fetching-Scripts erstellen
- `data/fetch/fetch_clinvar.sh` — ClinVar VCF + variant_summary.txt von NCBI FTP
- `data/fetch/fetch_snpedia.py` — SNPedia via MediaWiki API scrapen (respektiere Rate Limits: 3s Delay)
- `data/fetch/fetch_gwas.sh` — GWAS Catalog TSV-Download
- `data/fetch/fetch_dbsnp.sh` — dbSNP relevante Subset-Daten

### Aufgabe 2: SQLite-Schema & Import
- Einheitliches Schema in `data/schema/schema.sql`
- Import-Scripts die heruntergeladene Daten in SQLite transformieren
- Ziel: Eine einzelne `genesight.db` Datei (~500MB-1GB)

### Aufgabe 3: DNA-Parser
- `genesight-core` Parser für 23andMe, AncestryDNA, VCF
- Ergebnis: `Vec<Variant>` mit rsID, Chromosom, Position, Genotyp

### Aufgabe 4: Annotation & CLI
- Lookup jeder Variante gegen die lokale SQLite
- CLI-Interface: `genesight analyze <file> [--format json|md|html] [--tier 1|2|3]`
- Report-Output mit Confidence-Tiers
