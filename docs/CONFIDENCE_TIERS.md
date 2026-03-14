# GeneSight – Confidence-Tier-System

## Überblick

Jedes Ergebnis in GeneSight bekommt eine Zuverlässigkeitsstufe zugewiesen.
Das ist nicht optional — es ist eine Kernfunktion des Tools.

Das Tier-System existiert, weil DNA-Analyse kein binäres Ja/Nein ist.
Ein BRCA1-Mutationsnachweis hat eine fundamental andere Aussagekraft als
ein polygener Risikoscore für Depression. Nutzer müssen das auf den ersten
Blick erkennen können.

---

## Tier 1: Zuverlässig (>95% prädiktiver Wert)

### Was gehört hierhin?

**Monogenetische Erkrankungen** — Ein einzelnes Gen/eine einzelne Variante ist direkt kausal.

| Erkrankung | Gen | Vererbung | ClinVar-Status |
|-----------|-----|-----------|---------------|
| Sichelzellanämie | HBB | Autosomal rezessiv | Pathogenic (4★) |
| Mukoviszidose | CFTR | Autosomal rezessiv | Pathogenic (4★) |
| Huntington | HTT | Autosomal dominant | Pathogenic (4★) |
| BRCA1/2 Brustkrebs-Risiko | BRCA1, BRCA2 | Autosomal dominant | Pathogenic (4★) |
| Hämochromatose | HFE | Autosomal rezessiv | Pathogenic (3-4★) |
| Faktor V Leiden | F5 | Autosomal dominant | Pathogenic (4★) |

**Trägerstatus** — Heterozygote Träger rezessiver Erkrankungen.
Klinisch relevant für Familienplanung.

**Pharmakogenetik** — Wie der Körper Medikamente verstoffwechselt.

| Gen | Medikament(e) | PharmGKB Level | Effekt |
|-----|--------------|----------------|--------|
| CYP2D6 | Codein, Tramadol, Tamoxifen | 1A | Poor → keine Wirkung; Ultrarapid → Überdosis-Risiko |
| CYP2C19 | Clopidogrel, Omeprazol | 1A | Poor → reduzierte Wirkung |
| CYP2C9 + VKORC1 | Warfarin | 1A | Dosis-Anpassung nötig |
| HLA-B*5701 | Abacavir | 1A | Hypersensitivitäts-Reaktion |
| DPYD | 5-Fluorouracil | 1A | Schwere Toxizität bei Poor Metabolizer |

**Einfache Merkmale** — Wenige Gene, gut verstanden.

| Merkmal | Gen/SNP | Genauigkeit |
|---------|---------|-------------|
| Laktosetoleranz | MCM6 (rs4988235) | >95% |
| Ohrenschmalz-Typ | ABCC11 (rs17822931) | >95% |
| Asparagus-Geruch | Mehrere SNPs nahe OR2M7 | ~90% |
| Bitterschmecker (PTC) | TAS2R38 | >90% |

### Logik für Tier-1-Zuweisung

```rust
fn is_tier1(annotation: &Annotation) -> bool {
    // ClinVar: pathogenic + Review ≥ 2 Sterne
    if let Some(cv) = &annotation.clinvar {
        if cv.significance == "Pathogenic"
           && cv.review_stars >= 2 {
            return true;
        }
    }

    // PharmGKB: Evidence Level 1A oder 1B
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("1") {
            return true;
        }
    }

    // SNPedia: Magnitude ≥ 4 UND in ClinVar bestätigt
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude >= 4.0 && annotation.clinvar.is_some() {
            return true;
        }
    }

    false
}
```

---

## Tier 2: Wahrscheinlich (60-85% prädiktiver Wert)

### Was gehört hierhin?

**Polygene Risikoscores (PRS)** — Hunderte bis tausende Varianten wirken zusammen.
Kein einzelner SNP ist kausal, aber in Kombination ergibt sich ein statistisches Risikoprofil.

| Erkrankung/Trait | Anzahl SNPs | Erkl. Varianz | Aussagekraft |
|-----------------|-------------|---------------|-------------|
| Koronare Herzkrankheit | ~1.7M | ~15% | Oberes Quintil: 3x Risiko |
| Typ-2-Diabetes | ~400K | ~10% | Oberes Quintil: 2.5x Risiko |
| Bluthochdruck | ~900K | ~8% | Moderat prädiktiv |
| BMI/Adipositas | ~700K | ~6% | Schwach-moderat prädiktiv |

**Wichtig:** "Oberes Quintil: 3x Risiko" heißt NICHT "du bekommst das".
Es heißt: "In einer Population ist dein genetisches Risiko höher als bei 80% der Menschen."
Lebensstil, Ernährung und Umwelt sind oft wichtiger.

**Körperliche Merkmale mit mehreren Genen:**

| Merkmal | Genauigkeit |
|---------|-------------|
| Haarfarbe | ~70-85% |
| Sommersprossen | ~70-80% |
| Glatzenrisiko (männlich) | ~70-80% |
| Augenfarbe | ~75-90% (blau vs. braun gut, Mischfarben schwächer) |

### Logik für Tier-2-Zuweisung

```rust
fn is_tier2(annotation: &Annotation) -> bool {
    // GWAS: Signifikante Assoziation mit moderater Effektstärke
    if let Some(gwas) = &annotation.gwas_hits {
        let significant = gwas.iter().any(|g|
            g.p_value < 5e-8  // Genomweite Signifikanz
            && g.odds_ratio.map_or(false, |or| or > 1.1 && or < 3.0)
        );
        if significant { return true; }
    }

    // SNPedia: Magnitude 2-3.9, kein ClinVar-Eintrag
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude >= 2.0 && snp.magnitude < 4.0
           && annotation.clinvar.is_none() {
            return true;
        }
    }

    // PharmGKB: Evidence Level 2A oder 2B
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("2") {
            return true;
        }
    }

    false
}
```

---

## Tier 3: Spekulativ (50-65% prädiktiver Wert)

### Was gehört hierhin?

**Komplexe psychiatrische Erkrankungen:**

| Erkrankung | Erblichkeit (Zwillinge) | Erkl. durch SNPs | Einzelner SNP-Effekt |
|-----------|------------------------|-------------------|---------------------|
| Schizophrenie | ~80% | ~7% | Minimal |
| Bipolare Störung | ~70% | ~5% | Minimal |
| Major Depression | ~40% | ~2% | Winzig |
| Autismus | ~50-90% | ~3% | Minimal |
| ADHS | ~70-80% | ~3% | Minimal |

Die hohe Erblichkeit in Zwillingsstudien vs. niedrige erklärte Varianz durch
bekannte SNPs = "Missing Heritability". Gen-Umwelt-Interaktionen, Epigenetik,
und seltene Varianten spielen große Rollen.

**Persönlichkeit & Kognition:**

| Trait | SNP-Effekt | Nutzen eines Tests |
|-------|-----------|-------------------|
| Intelligenz (IQ) | ~0.01% pro SNP | Praktisch keiner |
| Risikobereitschaft | Minimal | Praktisch keiner |
| Neurotizismus | Minimal | Praktisch keiner |

**Sportliche Eignung:**

| Gen/SNP | Behauptung | Realität |
|---------|-----------|---------|
| ACTN3 (rs1815739) | "Sprinter-Gen" | Erklärt <1% der Varianz in sportl. Leistung |
| ACE I/D | Ausdauer vs. Kraft | Inkonsistente Studienlage |

### Logik für Tier-3-Zuweisung

```rust
fn is_tier3(annotation: &Annotation) -> bool {
    // GWAS: Schwache Assoziation
    if let Some(gwas) = &annotation.gwas_hits {
        let weak = gwas.iter().any(|g|
            g.p_value < 5e-8
            && g.odds_ratio.map_or(true, |or| or <= 1.1)
        );
        if weak { return true; }
    }

    // SNPedia: Niedrige Magnitude
    if let Some(snp) = &annotation.snpedia {
        if snp.magnitude > 0.0 && snp.magnitude < 2.0 {
            return true;
        }
    }

    // PharmGKB: Evidence Level 3 oder 4
    if let Some(pgkb) = &annotation.pharmacogenomics {
        if pgkb.evidence_level.starts_with("3")
           || pgkb.evidence_level.starts_with("4") {
            return true;
        }
    }

    false
}
```

---

## Report-Darstellung

### CLI-Output (Beispiel)

```
═══════════════════════════════════════════════════════════
  GeneSight Report — 2026-03-14
  Datei: meine_dna.txt (23andMe Format)
  Varianten analysiert: 637,291
  Annotierte Varianten: 12,847
═══════════════════════════════════════════════════════════

🟢 TIER 1 — Zuverlässig (klinisch validiert)
───────────────────────────────────────────────────────────

  PHARMAKOGENETIK
  ┌─────────────┬──────────────┬───────────────────────────┐
  │ Gen         │ Medikament   │ Status                    │
  ├─────────────┼──────────────┼───────────────────────────┤
  │ CYP2D6      │ Codein       │ ⚠ Poor Metabolizer        │
  │ CYP2C19     │ Clopidogrel  │ ✓ Normal Metabolizer      │
  └─────────────┴──────────────┴───────────────────────────┘
  Quelle: PharmGKB (Level 1A) + ClinVar

  TRÄGERSTATUS
  • Mukoviszidose (CFTR): Heterozygot — Träger, nicht betroffen
    Quelle: ClinVar (4★ Review)

  MERKMALE
  • Laktosetoleranz: Wahrscheinlich laktoseintolerant (CC bei rs4988235)
  • Ohrenschmalz: Trockener Typ

🟡 TIER 2 — Wahrscheinlich (statistische Assoziation)
───────────────────────────────────────────────────────────

  POLYGENE RISIKOSCORES
  • Koronare Herzkrankheit: 72. Perzentil (leicht erhöht)
    ⚠ Lebensstil hat größeren Einfluss als Genetik
  • Typ-2-Diabetes: 45. Perzentil (durchschnittlich)

  MERKMALE
  • Glatzenrisiko: Erhöht (68% Wahrscheinlichkeit bis 50)
  • Sommersprossen: Wahrscheinlich vorhanden

🔴 TIER 3 — Spekulativ (schwache Evidenz)
───────────────────────────────────────────────────────────

  ⚠ Die folgenden Ergebnisse haben geringe prädiktive
    Aussagekraft und sollten NICHT für Entscheidungen
    herangezogen werden.

  • ACTN3 (Muskeltyp): Mischtyp (Ausdauer/Kraft)
  • Koffein-Metabolismus: Schnell (CYP1A2 rs762551 AA)

───────────────────────────────────────────────────────────
⚕ DISCLAIMER: Dieser Report ist informativ, nicht
  diagnostisch. Für medizinische Entscheidungen
  konsultieren Sie einen Arzt oder genetischen Berater.
───────────────────────────────────────────────────────────
```

---

## Was GeneSight bewusst NICHT zeigt

1. **Keine Lebensdauer-Vorhersagen**
2. **Keine IQ-Scores oder "Intelligenz-Gene"**
3. **Keine Rassen-/Ethnizitäts-basierten Risikobewertungen** ohne Kontext
4. **Keine Ergebnisse ohne Confidence-Tier**
5. **Keine Ergebnisse die wie Diagnosen klingen** — immer probabilistisch formuliert
