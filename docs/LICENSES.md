# GeneSight – Lizenz-Übersicht

## Projekt-Lizenz

**GeneSight** steht unter der **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

### Warum GPL-3.0?

1. **Copyleft:** Stellt sicher, dass Forks und Ableitungen ebenfalls Open Source bleiben
2. **Kompatibilität:** GPL-3.0 ist kompatibel mit CC-BY-NC-SA 3.0 (SNPedia) im nicht-kommerziellen Kontext
3. **Community:** Fördert Beiträge zurück ins Projekt
4. **Privacy-Signal:** Unterstreicht, dass der Code inspizierbar ist — wichtig bei DNA-Daten

---

## Datenbank-Lizenzen im Detail

### Public Domain (keine Einschränkungen)

| Datenbank | Begründung |
|-----------|-----------|
| **ClinVar** | US Government Work — per 17 U.S.C. § 105 nicht urheberrechtsfähig |
| **dbSNP** | US Government Work — gleiche Begründung |

Diese Daten können ohne jede Einschränkung genutzt, verteilt und kommerziell verwendet werden.
Attribution ist nicht rechtlich erforderlich, aber wissenschaftlich gute Praxis.

### CC-BY-NC-SA 3.0 US (SNPedia)

**Erlaubt:**
- Teilen — Kopieren und Weiterverbreiten in jedem Format
- Bearbeiten — Remixen, Verändern und Aufbauen auf dem Material

**Bedingungen:**
- **BY (Namensnennung):** SNPedia muss als Quelle genannt werden
- **NC (Nicht-kommerziell):** Keine kommerzielle Nutzung ohne separate Lizenz
- **SA (Share-Alike):** Abgeleitete Werke müssen unter der gleichen Lizenz stehen

**Für GeneSight bedeutet das:**
- ✅ Das Open-Source-Projekt kann SNPedia-Daten frei nutzen
- ✅ Nutzer können das Tool für persönliche DNA-Analyse verwenden
- ✅ Akademische Forschung ist erlaubt
- ❌ Ein kommerzieller Fork müsste SNPedia-Daten entfernen
- → **Architektur-Entscheidung:** SNPedia-Daten werden als separater, optionaler Download behandelt, nicht im Repo gebündelt

### CC-BY-SA 4.0 (PharmGKB)

**Erlaubt:**
- Teilen und Bearbeiten, auch kommerziell

**Bedingungen:**
- **BY (Namensnennung):** PharmGKB muss als Quelle genannt werden
- **SA (Share-Alike):** Abgeleitete Werke unter gleicher oder kompatibler Lizenz

**Hinweis:** PharmGKB hat zusätzliche Nutzungsbedingungen für kommerzielle Nutzung.
Für akademische und nicht-kommerzielle Open-Source-Nutzung: frei verfügbar.

### Open Access (GWAS Catalog, gnomAD)

| Datenbank | Lizenz | Details |
|-----------|--------|---------|
| **GWAS Catalog** | EMBL-EBI Terms of Use | Frei für alle Zwecke, Attribution erbeten |
| **gnomAD** | ODC Open Database License | Frei für alle Zwecke inkl. kommerziell |

---

## Lizenz-Kompatibilitäts-Matrix

```
GPL-3.0 (GeneSight Code)
├── ✅ Public Domain (ClinVar, dbSNP) — kein Konflikt
├── ✅ CC-BY-SA 4.0 (PharmGKB) — kompatibel ab GPL-3.0
├── ✅ ODC-ODbL (gnomAD) — kompatibel
├── ✅ Open Access (GWAS Catalog) — kompatibel
└── ⚠️ CC-BY-NC-SA 3.0 (SNPedia) — kompatibel NUR wenn:
    - Das Gesamtprojekt nicht-kommerziell genutzt wird, ODER
    - SNPedia-Daten als separater, optionaler Download behandelt werden
```

### Architektur-Lösung für CC-BY-NC-SA

```
genesight.db (Haupt-Datenbank)
├── clinvar     → Public Domain ✅
├── dbsnp       → Public Domain ✅
├── gwas        → Open Access ✅
├── gnomad      → ODC-ODbL ✅
└── pharmgkb    → CC-BY-SA 4.0 ✅

snpedia.db (Separate, optionale Datenbank)
└── snpedia     → CC-BY-NC-SA 3.0 ⚠️
```

Das CLI-Tool funktioniert ohne `snpedia.db` — der Nutzer kann sie optional
herunterladen mit `genesight fetch --include-snpedia`.

---

## Attribution im Code

Jeder generierte Report MUSS folgenden Attribution-Block enthalten:

```markdown
---
## Datenquellen

Dieser Report wurde erstellt mit GeneSight (GPL-3.0).
Die folgenden Datenquellen wurden verwendet:

- **ClinVar** — NCBI, National Library of Medicine (Public Domain)
- **SNPedia** — snpedia.com (CC-BY-NC-SA 3.0) [falls genutzt]
- **GWAS Catalog** — NHGRI-EBI (Open Access)
- **gnomAD** — Broad Institute (ODC-ODbL)
- **PharmGKB** — pharmgkb.org (CC-BY-SA 4.0) [falls genutzt]

Dieser Report ist NICHT diagnostisch. Konsultieren Sie einen Arzt oder
genetischen Berater für medizinische Entscheidungen.
---
```

---

## Drittanbieter-Abhängigkeiten (Rust Crates)

Alle Rust-Crates müssen GPL-3.0-kompatibel sein. Erlaubte Lizenzen:

- MIT ✅
- Apache-2.0 ✅
- BSD-2-Clause / BSD-3-Clause ✅
- ISC ✅
- MPL-2.0 ✅ (Copyleft auf Datei-Ebene)
- GPL-3.0 ✅
- Unlicense ✅

**Nicht erlaubt:**
- GPL-2.0-only (ohne "or later") — inkompatibel mit GPL-3.0
- AGPL-3.0 — würde Server-Betreiber zu Source-Disclosure zwingen
- Proprietäre Lizenzen
