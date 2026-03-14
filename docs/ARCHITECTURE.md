# GeneSight – Architektur

## Design-Prinzipien

1. **Privacy first** — DNA-Daten verlassen nie den Rechner des Nutzers
2. **Offline-fähig** — Nach initialem Datenbank-Download keine Netzwerkverbindung nötig
3. **Modular** — Core-Library ist unabhängig von CLI/Desktop/Web
4. **Transparent** — Jedes Ergebnis zeigt Quelle und Konfidenz

---

## Crate-Struktur

```
genesight (workspace)
│
├── genesight-core          # Library — kein IO, kein Netzwerk, pure Logik
│   ├── parser::*           # DNA-Dateien → Vec<Variant>
│   ├── db::*               # SQLite-Queries (nimmt Connection als Parameter)
│   ├── annotator::*        # Variant + DB → AnnotatedVariant
│   ├── scorer::*           # AnnotatedVariant → ScoredResult + ConfidenceTier
│   ├── report::*           # Vec<ScoredResult> → Report (MD/JSON/HTML)
│   └── models::*           # Shared Types
│
├── genesight-cli           # Binary — IO, Argparse, DB-Öffnung
│   └── main.rs             # clap, rusqlite::Connection, ruft core auf
│
├── genesight-server        # Binary — Axum API (Phase 3)
│   └── main.rs             # Upload-Endpoint, ruft core auf
│
└── genesight-desktop       # Binary — Tauri App (Phase 2)
    └── ...                 # Tauri commands wrappen core
```

### Warum diese Trennung?

`genesight-core` hat **keine** Abhängigkeiten auf:
- Dateisystem-IO (bekommt `&[u8]` oder `&str`, nicht Pfade)
- Netzwerk (bekommt `rusqlite::Connection`, öffnet sie nicht selbst)
- CLI-Framework (kein clap)
- Async Runtime (alles synchron)

Das macht die Library:
- **Testbar** — Unit-Tests brauchen kein Dateisystem
- **Wiederverwendbar** — CLI, Tauri, Axum, WASM können alle dieselbe Library nutzen
- **Kompilierbar für WASM** — Zukunftsoption: im Browser laufen lassen

---

## Datenfluss

```
┌─────────────┐     ┌──────────┐     ┌───────────┐     ┌──────────┐     ┌──────────┐
│  DNA-Datei  │ ──► │  Parser  │ ──► │ Annotator │ ──► │  Scorer  │ ──► │  Report  │
│ (23andMe/   │     │          │     │           │     │          │     │ (MD/JSON │
│  VCF/etc.)  │     │ → Vec<   │     │ → Vec<    │     │ → Vec<   │     │  /HTML)  │
│             │     │ Variant> │     │ Annotated │     │ Scored   │     │          │
│             │     │          │     │ Variant>  │     │ Result>  │     │          │
└─────────────┘     └──────────┘     └─────┬─────┘     └──────────┘     └──────────┘
                                           │
                                    ┌──────┴──────┐
                                    │  SQLite DB  │
                                    │ (genesight  │
                                    │    .db)     │
                                    └─────────────┘
```

### Schritt für Schritt

1. **Parser** liest die DNA-Datei und produziert `Vec<Variant>`
   - Erkennt Format automatisch (23andMe vs AncestryDNA vs VCF)
   - Normalisiert auf einheitliches `Variant`-Struct
   - Filtert ungültige Zeilen (--calls, no-calls)

2. **Annotator** nimmt jede `Variant` und schlägt sie in der lokalen DB nach
   - Batch-Queries für Performance (nicht 600K Einzel-Queries)
   - Produziert `AnnotatedVariant` mit allen gefundenen Annotationen

3. **Scorer** bewertet jede annotierte Variante
   - Weist `ConfidenceTier` zu (Tier 1/2/3)
   - Berechnet polygene Risikoscores (summiert über viele Varianten)
   - Bestimmt Pharmakogenetik-Phänotypen (Metabolizer-Status)

4. **Report** formatiert die Ergebnisse
   - Gruppiert nach Tier und Kategorie
   - Fügt Disclaimer und Attributions hinzu
   - Output: Markdown, JSON, oder HTML

---

## Kern-Datenmodelle

```rust
/// Eine einzelne DNA-Variante aus der Nutzerdatei
pub struct Variant {
    pub rsid: Option<String>,       // z.B. "rs1234567"
    pub chromosome: String,          // "1"-"22", "X", "Y", "MT"
    pub position: u64,               // Genomische Position
    pub genotype: Genotype,          // z.B. Genotype::Heterozygous('A', 'G')
    pub source_format: SourceFormat, // TwentyThreeAndMe, AncestryDNA, VCF
}

pub enum Genotype {
    Homozygous(char),           // AA, GG, etc.
    Heterozygous(char, char),   // AG, CT, etc.
    NoCall,                      // --
    Indel(String),              // Insertion/Deletion
}

pub enum SourceFormat {
    TwentyThreeAndMe,
    AncestryDNA,
    Vcf,
}

/// Eine Variante mit allen Datenbank-Annotationen
pub struct AnnotatedVariant {
    pub variant: Variant,
    pub clinvar: Option<ClinVarAnnotation>,
    pub snpedia: Option<SnpediaAnnotation>,
    pub gwas_hits: Vec<GwasHit>,
    pub frequency: Option<AlleleFrequency>,
    pub pharmacogenomics: Option<PharmaAnnotation>,
}

/// Bewertetes Ergebnis mit Konfidenz
pub struct ScoredResult {
    pub variant: AnnotatedVariant,
    pub tier: ConfidenceTier,
    pub category: ResultCategory,
    pub summary: String,            // Menschenlesbare Zusammenfassung
    pub details: String,            // Ausführlichere Erklärung
}

pub enum ConfidenceTier {
    Tier1Reliable,    // >95% — klinisch validiert
    Tier2Probable,    // 60-85% — statistische Assoziation
    Tier3Speculative, // 50-65% — schwache Evidenz
}

pub enum ResultCategory {
    MonogenicDisease,
    CarrierStatus,
    Pharmacogenomics,
    PolygenicRiskScore,
    PhysicalTrait,
    ComplexTrait,
    Ancestry,
}
```

---

## Performance-Überlegungen

### DNA-Datei-Größen
- 23andMe: ~600K Varianten, ~15MB Textdatei
- AncestryDNA: ~700K Varianten, ~20MB
- WGS VCF: ~4-5M Varianten, ~1-5GB

### SQLite-Query-Strategie

**Nicht so:**
```rust
// ❌ 600.000 einzelne Queries
for variant in variants {
    db.query("SELECT * FROM clinvar WHERE rsid = ?", &[&variant.rsid]);
}
```

**Sondern so:**
```rust
// ✅ Batch-Queries mit temporärer Tabelle
db.execute("CREATE TEMP TABLE user_variants (rsid TEXT PRIMARY KEY)");
// Bulk-Insert der User-Varianten
// Dann JOINs gegen die Annotationstabellen
db.query("
    SELECT uv.rsid, c.*
    FROM user_variants uv
    LEFT JOIN clinvar c ON uv.rsid = c.rsid
");
```

Erwartete Performance:
- Parsing: <2 Sekunden für 600K Varianten
- Annotation (Batch): <10 Sekunden gegen alle Datenbanken
- Report-Generierung: <1 Sekunde
- **Gesamt: <15 Sekunden** für einen kompletten Report

---

## CLI-Interface (Phase 1)

```
genesight 0.1.0
Open-source DNA analysis tool

USAGE:
    genesight <COMMAND>

COMMANDS:
    analyze     Analyze a DNA file and generate a report
    fetch       Download/update reference databases
    info        Show database status and statistics
    help        Print help information

ANALYZE:
    genesight analyze <FILE> [OPTIONS]

    OPTIONS:
        -f, --format <FORMAT>    Output format [default: markdown] [possible: markdown, json, html]
        -o, --output <FILE>      Output file [default: stdout]
        -t, --tiers <TIERS>      Which tiers to include [default: 1,2,3]
        --db <PATH>              Path to genesight.db [default: ~/.genesight/genesight.db]
        --snpedia-db <PATH>      Path to snpedia.db (optional)
        --no-disclaimer          Omit disclaimer (for piping/scripting)
        -v, --verbose            Show all annotated variants, not just notable ones

FETCH:
    genesight fetch [OPTIONS]

    OPTIONS:
        --all                    Download all databases
        --clinvar                Download ClinVar
        --gwas                   Download GWAS Catalog
        --snpedia                Download SNPedia (CC-BY-NC-SA 3.0)
        --gnomad                 Download gnomAD frequencies
        --pharmgkb               Download PharmGKB
        --db-dir <PATH>          Database directory [default: ~/.genesight/]

INFO:
    genesight info [OPTIONS]

    OPTIONS:
        --db <PATH>              Path to genesight.db
```

---

## Zukunft: Tauri Desktop App (Phase 2)

```
┌─────────────────────────────────────────────────┐
│  Tauri Window (WebView)                         │
│  ┌───────────────────────────────────────────┐  │
│  │  SvelteKit / React Frontend               │  │
│  │  - Datei-Auswahl (Drag & Drop)           │  │
│  │  - Report-Anzeige (Tier-basiert)         │  │
│  │  - Datenbank-Management                   │  │
│  └──────────────────┬────────────────────────┘  │
│                     │ Tauri Commands (IPC)       │
│  ┌──────────────────┴────────────────────────┐  │
│  │  Rust Backend (genesight-core)            │  │
│  │  - Parser, Annotator, Scorer, Report      │  │
│  │  - SQLite DB Management                   │  │
│  │  - Optional: lokales LLM (Ollama)         │  │
│  └───────────────────────────────────────────┘  │
└─────────────────────────────────────────────────┘
```

Tauri-Commands sind dünne Wrapper:

```rust
#[tauri::command]
fn analyze_dna(file_path: String, db_path: String) -> Result<Report, String> {
    let data = std::fs::read_to_string(&file_path).map_err(|e| e.to_string())?;
    let conn = Connection::open(&db_path).map_err(|e| e.to_string())?;

    let variants = genesight_core::parser::parse_auto(&data)?;
    let annotated = genesight_core::annotator::annotate_batch(&conn, &variants)?;
    let scored = genesight_core::scorer::score_all(&annotated)?;
    let report = genesight_core::report::generate(&scored, Format::Json)?;

    Ok(report)
}
```
