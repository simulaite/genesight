# Wissenschaftlich korrekte Implementierung einer Consumer-Genomics-Annotation-Pipeline

Zielzustand: Jede Annotation wird erst nach **kanonischer Variantendarstellung** (Build + Koordinaten + REF/ALT + Strang) und **Allele-zu-Allele-Abgleich** erzeugt; rsID ist nur ein Alias und nie der Primärschlüssel für klinische Aussagen. Die sicherheitskritischen Fehler in der aktuellen Pipeline (rsID-Lookup ohne Allelvergleich) erzeugen systematisch Falschzuordnungen bei Pathogenität, Risikoallelen und Pharmakogenetik.

## Allele-Matching und Strang-Orientierung

**Wissenschaftlich korrekter Ansatz (Algorithmus + präzise Regeln)**  
Consumer-Rohdaten sind typischerweise bereits als Nukleotide (A/C/G/T) angegeben; entscheidend ist, **auf welche Referenz (Build) und welchen Strang** sich diese Nukleotide beziehen. Für zwei große DTC-Formate gilt: 23andMe berichtet Genotypen standardmäßig auf dem **Plus-Strang** relativ zu **GRCh37** (und bietet im Raw-Data-Browser auch GRCh38 auf Plus-Strang), und AncestryDNA berichtet Rohdaten auf dem **Forward/Plus-Strang** relativ zu **GRCh37**. citeturn14search0turn14search1 Das reduziert, aber eliminiert nicht, Strangprobleme, weil alle externen Referenzdatenbanken (und/oder historische Exporte) anders kodiert oder anders normalisiert sein können.

Auf der Array-/Manifest-Ebene existieren zusätzliche Strangdefinitionen: Illumina nutzt u.a. **TOP/BOT** (stabil gegenüber Build-Änderungen, nicht identisch mit +/-) sowie Manifest-Spalten wie **IlmnStrand** (TOP/BOT für SNPs bzw. PLUS/MINUS für Indels) und **SourceStrand** (eingereichte Kunden-/Datenbank-Strandangabe). citeturn0search10turn18search0turn21view3 Wenn deine Pipeline nicht nur DTC-Exports, sondern auch (i) generische Illumina-Final-Reports oder (ii) Manifest-basierte Formate unterstützt, muss sie diese Stranglagen explizit auflösen.

dbSNP ist als Primärquelle für rsID→(chr,pos,REF,ALT) geeignet, aber Historie zählt: NCBI beschreibt, dass RefSNPs früher teils als FWD/REV relativ zur Assembly geführt wurden und dass in der neueren dbSNP-Architektur die Allele konsistent „forward to the reported sequence“ (VCF/HGVS/SPDI-konform) berichtet werden; darum existieren auch spezielle VCFs für ehemals „REV“-gemeldete rsIDs. citeturn21view4turn13search8 Für robuste Allelvergleiche ist die korrekte Regel: **Allele auf denselben Build und dieselbe (Plus-)Orientierung normalisieren und erst dann vergleichen.**

ClinVar liefert Varianten u.a. als VCF über FTP; diese VCFs enthalten „simple alleles“ (<10kb) mit präzisen Endpunkten, gemappt auf GRCh37 oder GRCh38 (nicht alle ClinVar-Variantentypen sind im VCF). citeturn13search9turn13search17 Jede ClinVar-Aussage, die du anzeigst, muss daher (a) die richtige Assembly treffen und (b) deine Nutzerallele mit den in ClinVar referenzierten Allelen (REF/ALT bzw. CLNSIG-bezogene Allele) abgleichen.

Die GWAS Catalog „Top Hits“-Dateien enthalten ein Feld „STRONGEST SNP-RISK ALLELE“ (SNP + Risiko-/Effektallel; „?“ falls unbekannt) sowie ein Feld „OR or BETA“ zum Effekt; für ältere Kurationen wurde OR<1 teils invertiert und das berichtete Allel entsprechend umgedreht, damit die gespeicherten ORs >1 sind. citeturn7search0turn16search0turn16search12 Daraus folgt: Du darfst Risikoallele nicht blind übernehmen, sondern musst die jeweilige Katalog-Logik (Datum/Version) und die Allelorientierung prüfen oder auf harmonisierte Summary-Stats ausweichen.

image_group{"layout":"carousel","aspect_ratio":"16:9","query":["Illumina TOP BOT strand diagram","Illumina Infinium manifest IlmnStrand SourceStrand explanation","palindromic SNP A/T C/G strand ambiguity diagram","SNP strand flip plus minus complement A T C G diagram"],"num_per_query":1}

**Strang-/Allelnormalisierung: Referenz-Algorithmus (pseudocode)**  
Kernprinzip: Der Vergleich erfolgt auf **(chr,pos,REF,ALT)** in einer festgelegten Assembly. rsID ist nur Lookup-Hilfe.

```pseudo
# Inputs:
#   user: (assembly, chr, pos, genotype="AG"/"--"/"A", provider_hint)
#   refdb: (assembly, chr, pos, REF, ALTs[])  # aus dbSNP/VCF/FASTA verifiziert
#   target_allele: ein Nukleotid (A/C/G/T) oder bei Indels eine Sequenz (VCF-ALT)

function complement_base(b):
    map = { "A":"T", "T":"A", "C":"G", "G":"C" }
    return map[b]

function is_palindromic_snp(REF, ALT):
    # biallelic SNP palindromic wenn {A,T} oder {C,G}
    s = set([REF, ALT])
    return (s == {"A","T"} or s == {"C","G"})

function normalize_user_alleles_to_plus(user, refdb):
    # 1) Build/Koordinate muss schon refdb-assembly entsprechen (ggf. liftover vorher)
    assert user.assembly == refdb.assembly
    allowed = set([refdb.REF] + refdb.ALTs)

    alleles = parse_genotype(user.genotype)   # "AG"->["A","G"], "--"->[None,None]
    if any(a not in {"A","C","G","T",None} for a in alleles):
        return {status:"unsupported_allele_encoding"}

    if None in alleles:
        return {status:"missing_genotype"}

    # 2) Direktmatch?
    if set(alleles) ⊆ allowed:
        return {status:"ok", oriented:"+", alleles:alleles}

    # 3) Reverse-complement match?
    comp = [complement_base(a) for a in alleles]
    if set(comp) ⊆ allowed:
        # Ambiguität bei palindromischen SNPs: Flip nicht beweisbar ohne Zusatzinfo
        if len(refdb.ALTs)==1 and is_palindromic_snp(refdb.REF, refdb.ALTs[0]):
            return {status:"ambiguous_palindromic_strand"}
        return {status:"ok", oriented:"-", alleles:comp}

    return {status:"allele_mismatch"}  # falsche Koordinate, falscher Build, multi-allelic, rsID-Merge, etc.

function count_target_allele(alleles, target_allele):
    # 0/1/2 Kopien für biallelic SNPs; bei multi-allelic analog
    return sum(1 for a in alleles if a == target_allele)
```

Dieses Schema ist minimal; in der Praxis kommen hinzu: multi-allele rsIDs, Indels, und die Notwendigkeit, **REF gegen Referenz-FASTA zu verifizieren** (oder REF/ALT direkt aus einem verlässlichen VCF zu übernehmen). citeturn12search2turn12search14turn13search8

**Palindromische SNPs (A/T, C/G): sichere Behandlung**  
Bei A/T- oder C/G-SNPs ist „flip vs. nicht flip“ anhand der Allelmenge allein unentscheidbar. Das ist der klassische Fehlerbereich bei GWAS/Array-Datenharmonisierung. Werkzeuge wie **snpflip** markieren reverse/ambiguous SNPs, und **Genotype Harmonizer** kann ambige A/T und G/C SNPs über LD-Patterns gegen ein Referenzpanel ausrichten (in einem DTC-Local-Tool meist zu schwergewichtig; aber als Konzept wichtig). citeturn18search1turn18search7 Für Chip-spezifische Auflösung existieren außerdem kuratierte „strand and build files“ für viele Genotyping-Chips. citeturn18search22

**Datenquellen/Implementierungsartefakte (konkret, maschinenlesbar)**  
```text
# DTC Strand/Build
https://eu.customercare.23andme.com/hc/en-us/articles/115002090907-Raw-Genotype-Data-Technical-Details
https://support.ancestry.com/articles/en_GB/Support_Site/Downloading-DNA-Data

# Illumina Manifest/Strand-Definition
https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565
https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001489
https://www.illumina.com/documents/products/technotes/technote_topbot.pdf

# dbSNP / NCBI Variation Services (rsID↔SPDI/VCF)
https://www.ncbi.nlm.nih.gov/core/assets/snp/docs/RefSNP_orientation_updates.pdf
https://github.com/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/spdi_batch.py
https://pmc.ncbi.nlm.nih.gov/articles/PMC7523648/   # SPDI paper

# ClinVar FTP Primer (VCF scope, GRCh37/38)
https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/

# Strand-Resolution Tools / Chip-Strand-Files
https://github.com/andymckenzie/snpflip
https://pmc.ncbi.nlm.nih.gov/articles/PMC4307387/   # Genotype Harmonizer
https://www.well.ox.ac.uk/~wrayner/strand/
```
citeturn14search0turn14search1turn0search2turn18search0turn21view3turn21view4turn13search0turn13search8turn13search9turn18search1turn18search7turn18search22

**Häufige Fallstricke (consumer-genomics-typisch)**  
(1) rsID-basierte Anzeige von „Pathogenic/Risk“ ohne Allelvergleich produziert garantierte False Positives (dein aktueller Bug ist das Standardmuster). citeturn24view2turn21view4  
(2) Build-Mix (GRCh37 vs GRCh38) verursacht scheinbare „Allele mismatch“ und falsche Annotationen; Liftover funktioniert zwar oft sehr gut, kann aber bei bestimmten Variantentypen/Regionen scheitern, daher muss „unmappable“ explizit behandelt werden. citeturn18search13turn13search8  
(3) Palindromische SNPs ohne Zusatzinfo flippen → silent corruption von Risikoallelen. citeturn18search7turn18search1  
(4) Illumina A/B-Allele oder TOP/BOT ungeprüft als A/C/G/T interpretieren. citeturn21view3turn18search0  
(5) ClinVar-VCF als „vollständige ClinVar-Wahrheit“ verwenden (ist sie nicht; Teile liegen nur im XML/TXT-Full-Release). citeturn13search9turn13search17

**Priorität (Patientensicherheit)**  
Kritisch. Ohne korrekten Strang-/Allelabgleich sind alle nachfolgenden Kategorien (ClinVar, GWAS, PGx) potenziell systematisch falsch.

## PGx-Star-Allele-Calling und Phänotypisierung

**Wissenschaftlich korrekter Ansatz (CPIC/PharmVar-konform)**  
Pharmakogenomische Phänotypen sind i.d.R. nicht „ein SNP → ein Phänotyp“, sondern **Haplotypen (Star Alleles) → Diplotyp → Funktion/Activity Score → Phänotyp → Drug-Guidance**. CPIC beschreibt explizit, dass die Kombination der Allele den Diplotyp (Genotyp) bestimmt und daraus Phänotypklassen abgeleitet werden; für CYP2C19 sind z.B. *1/*17 = Rapid und *17/*17 = Ultrarapid, während *2/*17 trotz *17 als Intermediate klassifiziert wird. citeturn21view0turn27search2

**Warum dein aktueller Ansatz zwangsläufig falsch ist:** rs12248560 ist ein Definierer für CYP2C19*17; ohne Allelvergleich und ohne Diplotyp-Logik wird *1/*1 fälschlich als „Ultrarapid“ etikettiert. PharmVar zeigt *17 als -806C>T (rs12248560) und listet die zugehörigen Kernvarianten. citeturn1search2turn21view0

**Star-Allele-Definitionen: robuste Datenquellen statt Hardcoding**  
Hardcoding von „ein paar bekannten rsIDs“ ist nicht wartbar und wird bei Updates der Star-Nomenklatur schnell falsch (Star-Definitionen sind dynamisch). citeturn15search17turn26search3 CPIC stellt dafür maschinenlesbare Tabellen bereit: **allele_definition**, **allele_functionality_reference**, **diplotype_phenotype**, **frequency**, **gene_cds**. citeturn27search4turn28search15 Das ist der richtige Pfad für ein lokales Tool: Tabellen versionieren, offline shippen, transparent updaten.

**Konkrete, häufig verwendete Definierervarianten (Beispiele, nicht vollständig)**  
Die folgenden Beispiele dienen als Einstieg/Validierung; produktiv sollte die Definition immer aus den CPIC/PharmVar-Tabellen kommen:

- CYP2C19: *17 rs12248560 (Promoter -806C>T), *2 rs4244285, *3 rs4986893. citeturn1search2turn21view0turn27search13  
- CYP2C9: *2 rs1799853, *3 rs1057910. citeturn25search0turn25search5  
- CYP3A5: *3 rs776746 (Splice-Defekt). citeturn25search6turn25search2  
- DPYD: *2A rs3918290 (Splice-Defekt); CPIC arbeitet u.a. mit c.1905+1G>A, c.1679T>G, c.2846A>T, c.1129–5923C>G als zentralen Varianten für DPD-Aktivitätsabschätzung. citeturn25search7turn23view1  
- SLCO1B1: *5 rs4149056 (c.521T>C; V174A). citeturn26search7turn23view2turn26search1  
- TPMT: *3A (rs1800460 + rs1142345 in cis), *3B rs1800460, *3C rs1142345, *2 rs1800462; *3A ist ohne Phasing potentiell mehrdeutig. citeturn26search6turn26search0turn26search10turn26search4  
- NUDT15: *3 rs116855232 (R139C). citeturn25search4turn25search20  
- CYP2D6: sehr viele Allele (PharmVar nennt >130 Kernallele); einzelne Schlüsselvarianten (z.B. *4 rs3892097, *10 rs1065852) sind nur ein kleiner Teil, strukturelle Varianten (Deletion/Duplikation/Hybrid) sind klinisch relevant. citeturn26search3turn26search9turn21view2

**Korrekte End-to-End-Logik (SNP-Genotyp → Star Alleles → Diplotyp → Activity Score → Phänotyp)**  
Für ein lokales Consumer-Tool mit Array/VCF-Input ist die robuste Implementierung:

1) **Input harmonisieren**: Varianten auf eine Assembly bringen (GRCh37 oder GRCh38), Allele auf Plus-Strang normalisieren, REF gegen Referenz prüfen. citeturn14search0turn14search1turn12search2turn13search8  
2) **Allele-Definition laden**: CPIC allele_definition_table für das Gen (Star-Allele als Variant-Kombinationen). citeturn27search0turn28search11  
3) **Named-Alele-Matching**: Für jede Star-Definition prüfen, ob die Nutzer-Genotypdaten diese Haplotypdefinition zulassen (mit Missingness-Tracking). PharmCAT macht das als „Named Allele Matcher“ auf VCF-Input. citeturn1search5turn1search1turn1search13  
4) **Diplotyp-Inferenz ohne Phasing**: Möglichkeitsraum aller Haplotyppaare bilden, die die (ungephaseten) Genotypcounts erklären; bei Mehrdeutigkeit „ambiguous call“ statt erzwungenem Ergebnis. citeturn26search0turn25search1  
5) **Funktion + Activity Score**: CPIC allele_functionality_reference + diplotype_phenotype Tabellen verwenden. Für CYP2D6 wird pro Allel ein Aktivitätswert vergeben, bei Duplikationen multipliziert, dann summiert (Activity Score). citeturn21view2turn28search1turn28search0turn28search2  
6) **Phänotyp & Drug Guidance**: Gen-spezifische CPIC-Empfehlungslogik/Tabellen; für viele Gene existieren CPIC „gene_cds“ Tabellen als strukturierte Entscheidungsunterstützung. citeturn27search6turn28search7turn23view1turn23view2turn31search0

**Pseudocode (diplotype + phenotype, datengetrieben)**  
```pseudo
# datengetrieben: nicht hardcodieren, sondern CPIC-Tabellen laden

function call_gene_pgx(sample_vcf, gene):
    defs  = load_cpic_allele_definitions(gene)        # allele_definition_table.xlsx
    func  = load_cpic_allele_function(gene)           # allele_functionality_reference.xlsx
    d2p   = load_cpic_diplotype_to_phenotype(gene)    # Diplotype_Phenotype_Table.xlsx

    # 1) relevante Varianten extrahieren
    g = subset_variants(sample_vcf, defs.all_required_sites)

    # 2) mögliche haplotypen bestimmen (kompatibel mit beobachteten Genotypen)
    possible_haplotypes = []
    for star in defs.star_alleles:
        if star_is_compatible(star, g):   # berücksichtigt missingness + Konflikte
            possible_haplotypes.append(star)

    # 3) diplotype enumerieren (ungephaset)
    candidates = []
    for a in possible_haplotypes:
        for b in possible_haplotypes:
            if diplotype_explains_genotypes(a, b, g):
                completeness = completeness_score(a, b, g)
                candidates.append((a, b, completeness))

    if candidates empty:
        return {status:"no_call"}

    best = select_best_with_ties(candidates)
    if best.is_ambiguous:
        return {status:"ambiguous", diplotypes:best.list}

    diplotype = best.diplotype  # (a,b)
    phenotype = d2p.lookup(diplotype)  # bevorzugt: CPIC Tabelle, nicht selbst rechnen
    return {status:"ok", diplotype:diplotype, phenotype:phenotype}
```

**Missing SNPs: sichere Default-Regeln**  
„Nicht getestet“ ist nicht „Wildtyp“. CPIC weist in mehreren Guidelines darauf hin, dass Genotyp-Only-Tests seltene/neu entdeckte Varianten nicht erfassen; daraus folgt: wenn definierende Sites fehlen, ist ein „Normal“ oft nicht ableitbar und muss als „Unklar/Indeterminate“ zurückgegeben werden. citeturn21view2turn23view0turn28search0 Praktisch: Für jedes Gen eine **Coverage-Metrik** ausgeben („X% der definierenden Sites für häufige Core-Allele beobachtet“), und Phänotyp nur ausgeben, wenn die CPIC-Tabellen den Call als determiniert erlauben. citeturn27search2turn28search2

**CNVs / strukturelle Varianten (v.a. CYP2D6): harte Grenzen von Consumer-Arrays**  
CYP2D6-Phänotyp hängt stark von Deletionen und Duplikationen („xN“) ab; CPIC beschreibt Deletion (*5) und Duplikationen explizit und dass Aktivitätswerte bei Mehrkopien multipliziert werden. citeturn21view2turn28search1 Viele Consumer-Arrays detektieren diese CNVs nicht zuverlässig; daher muss dein Tool für CYP2D6 standardmäßig einen **„CNV unknown“-Status** führen und keine Ultrarapid-Calls behaupten, wenn keine Kopienzahlinfo vorhanden ist. Tools wie Stargazer können prinzipiell auch mit SNP-Array-Daten arbeiten, aber CNV-Auflösung bleibt je nach Input limitiert; Aldy/Astrolabe sind primär auf Sequenzdaten und strukturelle Modelle ausgelegt. citeturn15search4turn15search3turn15search2

**CPIC-Empfehlungen: hochrelevante Beispiele (datengetrieben referenzieren, nicht paraphrasieren)**  
- CYP2C19–Clopidogrel: CPIC 2022 klassifiziert *1/*17 als Rapid, *17/*17 als Ultrarapid und *2/*17 als Intermediate; Therapieempfehlungen leiten daraus ab. citeturn21view0turn20view0  
- DPYD–Fluoropyrimidine: CPIC nutzt DPYD Activity Score; Tabelle 2 gibt dosisreduzierte Startdosen für partielle DPD-Defizienz (z.B. 50% Reduktion bei bestimmten AS-Konstellationen) und betont Titration/Monitoring. citeturn23view1  
- TPMT/NUDT15–Thiopurine: Table 1 ordnet Diplotypen zu Phänotypen, Empfehlungen reduzieren Startdosen je nach Kombinationsstatus und erlauben Dosisanpassung nach Myelosuppression. citeturn23view0  
- SLCO1B1–Simvastatin: rs4149056 Genotyp wird in Phänotypklassen übersetzt; Empfehlungen vermeiden 80 mg und schlagen niedrigere Dosen oder Alternativstatin vor, besonders bei C-Allelen. citeturn23view2turn17search3  
- CYP2D6–Codein/Tramadol: CPIC empfiehlt, Codein/Tramadol bei Ultrarapid nicht zu nutzen (Toxizitätsrisiko) und Alternativen bei Poor; die 2021-Opioid-Guideline formuliert das mit Activity-Score-Schwellen. citeturn31search0turn31search4

**Datenquellen/Implementierungsartefakte (CPIC „Full Tables“, maschinenlesbar)**  
```text
# CPIC: CYP2C19 (Definition/Funktion/Diplotyp→Phänotyp)
https://files.cpicpgx.org/data/report/current/allele_definition/CYP2C19_allele_definition_table.xlsx
https://files.cpicpgx.org/data/report/current/allele_function_reference/CYP2C19_allele_functionality_reference.xlsx
https://files.cpicpgx.org/data/report/current/diplotype_phenotype/CYP2C19_Diplotype_Phenotype_Table.xlsx
https://cpicpgx.org/gene/cyp2c19/

# CPIC: CYP2D6 (Definition/Funktion/Diplotyp→Phänotyp)
https://files.cpicpgx.org/data/report/current/allele_definition/CYP2D6_allele_definition_table.xlsx
https://files.cpicpgx.org/data/report/current/allele_function_reference/CYP2D6_allele_functionality_reference.xlsx
https://files.cpicpgx.org/data/report/current/diplotype_phenotype/CYP2D6_Diplotype_Phenotype_Table.xlsx
https://files.cpicpgx.org/data/report/current/gene_phenotype/CYP2D6_phenotypes.xlsx

# PharmCAT (Referenzimplementierung)
https://pharmcat.clinpgx.org/using/
https://github.com/PharmGKB/PharmCAT-tutorial
https://pmc.ncbi.nlm.nih.gov/articles/PMC10121724/
```
citeturn27search0turn27search1turn27search2turn27search4turn28search11turn28search1turn28search0turn28search2turn1search1turn1search5turn1search13

**Häufige Fallstricke**  
(1) „Absent in data“ als *1 interpretieren → falsche Normal-/Rapid-/Ultrarapid-Calls. citeturn23view0turn21view2  
(2) TPMT*3A ohne Phasing erzwingen → falscher Diplotyp (cis/trans). citeturn26search0turn26search4  
(3) CYP2D6 ohne CNV-Status als vollständig reporten → Hochrisiko-Fehleinschätzungen. citeturn21view2turn15search9  
(4) Star-Definitions-Updates ignorieren → Drift zwischen Report und aktueller CPIC/PharmVar-Nomenklatur. citeturn15search17turn28search24

**Priorität (Patientensicherheit)**  
Hoch bis kritisch. Falsche PGx-Phänotypen können zu realer Medikamentenfehlsteuerung führen; FDA warnt explizit vor nicht validierten PGx-Claims und auch Software-Interpretationen im DTC-Kontext. citeturn30view0turn9search0

## ClinVar-Pathogenität korrekt interpretieren

**Wissenschaftlich korrekter Ansatz (nicht „contains pathogenic“)**
ClinVar ist ein Archiv eingereichter Interpretationen (SCVs), aggregiert zu Varianten-/Varianten-Konditions-Records (VCV/RCV). citeturn32search12turn24view2 Seit 2024 trennt ClinVar klinische Klassifikationstypen (germline, somatic clinical impact, oncogenicity) in getrennte Felder; „clinical_significance“ muss daher kontextualisiert werden. citeturn24view2

Die korrekte Logik ist dreistufig:

1) **Allele prüfen**: Nur wenn der Nutzer die in ClinVar klassifizierte(n) ALT-Allele trägt, ist die ClinVar-Klassifikation überhaupt genotypebezogen relevant (sonst „homozygous reference“/„kein Träger“). citeturn13search9turn21view4  
2) **Klassifikationstyp + Term interpretieren**: ClinVar-Terms umfassen u.a. Pathogenic/Likely pathogenic/Benign/Likely benign/VUS sowie „risk factor“, „drug response“, „association/protective/other“. citeturn24view1turn24view2  
3) **Evidenzqualität (Review Status/Stars) gewichten**: ClinVar definiert Stars/Review-Status: 4 = practice guideline, 3 = expert panel, 2 = multiple submitters/no conflicts, 1 = criteria provided (single submitter) oder criteria provided (conflicting), 0 = no criteria/no classification usw. citeturn24view0

**Zuverlässigkeitsschwellen (sicherheitsorientiert, ClinVar-konform)**  
Für ein Consumer-Tool ist eine konservative Policy notwendig:

- **Hochvertrauen**: 3–4 Sterne (Expert Panel/Practice Guideline). citeturn24view0turn32search5  
- **Mittel**: 2 Sterne (mehrere Submitter, Kriterien, kein Konflikt). citeturn24view0  
- **Niedrig/Informativ**: 1 Stern (Single submitter mit Kriterien) – anzeigen, aber deutlich als nicht-konsentiert kennzeichnen. citeturn24view0  
- **Konfliktfall**: 1 Stern „criteria provided, conflicting classifications“ – nie als eindeutige Pathogenität ausgeben; stattdessen Konfliktstruktur anzeigen (wer sagt was, welcher Review-Status, welches Datum). citeturn24view0turn24view2  
- **0 Sterne / no criteria**: nicht als klinische Aussage präsentieren; maximal als Rohhinweis. citeturn24view0

**Zygosität + Vererbung: korrekte, programmierbare Einbindung**  
Genotyp (0/1/2 ALT-Kopien) ist für monogene Erkrankungen nur im Kontext der **Vererbung (AD/AR/X-linked/mitochondrial)** interpretierbar. Grundlagen zu Vererbungsmodi sind standardisiert (autosomal dominant/recessive, X-linked usw.). citeturn32search15 Programmtaugliche Quellen für „Mode of Inheritance“ (MOI):

- ClinVar kann MOI auf Submission-Ebene enthalten; wenn Submitter MOI angibt, wird das auf Variantenseiten angezeigt. citeturn32search1turn24view2  
- ClinVar „properties“ und Filterbegriffe enthalten MOI-Kategorien (moi autosomal dominant/recessive usw.). citeturn32search8  
- MedGen unterstützt „mode of inheritance“ als Such-/Property-Feld und verarbeitet u.a. Orphanet/ORDO inklusive MOI. citeturn32search0turn32search4  
- ClinGen Gene-Disease Validity Knowledge Base führt „Mode of Inheritance“ pro gene–disease Assertion und ist öffentlich einsehbar. citeturn2search19turn2search7

**Korrekte ClinVar-Interpretationslogik (pseudocode, risikoarm)**  
```pseudo
# Inputs:
#   variant_call: {user_count_alt, genotype_quality, allele_normalized_ok}
#   clinvar_records: list of RCV-like entries {condition_id, clinsig_term, review_status, last_eval_date, moi?}
# Output:
#   structured interpretation objects (nicht "Diagnose", sondern Variant-Hinweis)

function interpret_clinvar(variant_call, clinvar_records):
    if variant_call.allele_normalized_ok != true:
        return {status:"cannot_interpret_without_allele_match"}

    if variant_call.user_count_alt == 0:
        return {status:"non_carrier"}  # keine ALT-/Risikoallel-Kopie

    # filtere auf germline classification (nicht somatic clinical impact / oncogenicity)
    records = filter_germline_records(clinvar_records)

    # gewichtete Auswahl nach Review Status
    # 4>3>2>1>0; Konfliktstatus separat
    best = pick_by_highest_review_status(records)

    # Wenn best "conflicting classifications": niemals binär entscheiden
    if best.review_status_contains("conflicting"):
        return {status:"conflicting", details: summarize_conflict(records)}

    # clinsig-Terminologie korrekt behandeln
    if best.clinsig in {"Pathogenic","Likely pathogenic"}:
        return {status:"P_or_LP", zygosity:variant_call.user_count_alt, moi:best.moi, caveats:...}
    if best.clinsig in {"Benign","Likely benign"}:
        return {status:"B_or_LB", ...}
    if best.clinsig in {"Uncertain significance"}:
        return {status:"VUS", ...}
    if best.clinsig in {"risk factor","association","protective","drug response","other"}:
        return {status:"non_mendelian_or_pgx_term", ...}
```

**ACMG/AMP-Framework: was es ist und wie es zu ClinVar passt**  
ACMG/AMP (Richards et al., 2015) definiert die Standardterminologie (pathogenic/likely pathogenic/VUS/likely benign/benign) und evidenzbasierte Kriterien inkl. Populationsfrequenz-Regeln (BA1/BS1/PM2 etc.). citeturn2search2turn8search1turn8search5 ClinVar-Terms orientieren sich an dieser Terminologie, aber ClinVar ist kein eigener ACMG-Klassifikator; es aggregiert Submitter-Aussagen und gewichtet Aggregation u.a. nach Review-Status. citeturn24view2turn24view0 Die ClinVar-Star-Systematik ist daher **Evidenz-/Prozessqualität**, nicht ein formales Mapping auf ACMG-Kriterienerfüllung. citeturn24view0turn2search2

**Datenquellen/Implementierungsartefakte**  
```text
# ClinVar Doku: Review Status / Sterne, Klassifikationstypen, Termini
https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
https://www.ncbi.nlm.nih.gov/clinvar/docs/properties/
https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/

# ACMG/AMP 2015 (Primary Source)
https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/

# MedGen MOI
https://www.ncbi.nlm.nih.gov/medgen/docs/search/
https://www.ncbi.nlm.nih.gov/medgen/docs/data/

# ClinGen Gene-Disease Validity KB (MOI sichtbar)
https://search.clinicalgenome.org/
```
citeturn24view0turn24view2turn24view1turn13search9turn2search2turn32search0turn32search4turn2search19

**Häufige Fallstricke**  
(1) „Pathogenic“ anzeigen trotz 0 Kopien ALT; das ist exakt dein rsID-only Bug. citeturn13search9turn24view1  
(2) „Conflicting interpretations“ wie „Pathogenic“ behandeln; ClinVar unterscheidet Konflikte explizit im Review-Status. citeturn24view0turn24view1  
(3) Somatische/onko-Klassifikation als germline „Pathogenic“ ausgeben; ClinVar trennt diese Klassifikationstypen. citeturn24view2  
(4) MOI/zygosity ignorieren → AR-Erkrankungen fälschlich als erkrankt bei heterozygoter Trägerschaft. citeturn32search15turn32search0

**Priorität (Patientensicherheit)**  
Kritisch. Mendelian „Pathogenic“-Labels ohne Allel- und Review-Status-Kontrolle sind der schnellste Weg zu hochschädlichen False Alarms.

## GWAS-Risikointerpretation und PRS

**Wissenschaftlich korrekter Ansatz (Single-Variant und Multi-Variant)**  
GWAS-Top-Hits sind Assoziationen mit meist kleinen Effekten; korrekte Nutzer-Interpretation benötigt: (1) Effekt-/Risikoallel, (2) Effektgröße (OR oder Beta) bezogen auf dieses Allel, (3) Genotyp (0/1/2 Kopien), (4) ggf. Populationsfrequenzen und (5) Baseline-Prävalenz für absolute Risiken. citeturn7search0turn7search20turn37search1

**Welche Allele im GWAS Catalog „Risiko“-Allele sind (und warum das heikel ist)**  
In den kuratierten GWAS Catalog Top-Hits bezeichnet „STRONGEST SNP-RISK ALLELE“ die Variante plus Risiko-/Effektallel; „?“ wenn unbekannt. citeturn7search0turn7search2 Zusätzlich ist die Effektspalte „OR or BETA“ kontextabhängig; für vor Jan 2021 kuratierte Studien wurde OR<1 teils invertiert und das berichtete Allel entsprechend gedreht, damit OR>1 gespeichert wird. citeturn16search0turn16search12 Daraus folgt: Das „Risikoallel“ in Top-Hits ist nicht automatisch ein konsistent auf Plus-Strang harmonisiertes Effektallel in VCF-Sinn.

Für robuste, skalierbare Interpretation sind harmonisierte Summary-Statistics (GWAS-SSF) vorzuziehen, da sie explizit **effect_allele/other_allele, effect size, SE, p** usw. erfassen und durch Pipelines harmonisiert werden können. citeturn7search8turn3search5turn7search5

**Genotyp → Risikoallel-Kopien (0/1/2): korrekte Logik**  
Nach Strang-/Build-Normalisierung ist die GWAS-Seite analog zu ClinVar:

- 0 Kopien Effektallel: Referenz für den (studien-)definierten Baselinevergleich.  
- 1 Kopie: heterozygot.  
- 2 Kopien: homozygot.  

Dieses Zählen ist nur gültig, wenn Effektallel und Nutzerallele auf derselben Referenzorientierung liegen. citeturn18search7turn7search15turn7search5

**OR korrekt in personale (relative/absolute) Risikomaße übersetzen**  
OR ist ein Odds-Verhältnis, nicht direkt „Risiko“. Für eine additive log-Odds-Annahme pro Effektallel gilt typischerweise: **Odds_multiplier = OR^k** (k = 0/1/2 Effektallele). citeturn7search20turn3search17 Für eine absolute Risikoapproximation brauchst du eine Baseline-Risikoannahme p0 (z.B. Prävalenz/Lebenszeitrisiko im Zielkollektiv): odds0 = p0/(1-p0); odds = odds0 * OR^k; p = odds/(1+odds). Die Beziehung OR↔RR hängt stark von p0 ab; bei hohen Baseline-Raten divergieren OR und RR stark. citeturn3search2turn3search10turn3search17 Ohne p0 ist „absolutes Risiko“ nicht seriös berechenbar.

**Beta (quantitative Traits) anders behandeln als OR (binär)**  
Beta ist ein additiver Effekt pro Effektallel auf der Trait-Skala (oder einer transformierten Skala), darum ist die naive Erwartungsverschiebung **ΔTrait ≈ beta * k**; zusätzlich braucht man Einheit/Skalierung und Populationsreferenz (Mittelwert/SD), um Patientenverständnis zu ermöglichen. citeturn3search5turn7search20

**Reporting-Schwellen: was wissenschaftlich vertretbar ist (Consumer-Kontext)**  
GWAS-Literatur betont die Notwendigkeit standardisierter Berichtsbestandteile (inkl. Allele/Strand/Effect sizes) und dass Mindestinformationen vorhanden sein müssen. citeturn7search20turn3search5 Für Consumer-Reporting ist als Minimalfilter vertretbar: nur Assoziationen mit genome-wide significance (klassisch p<5×10^-8) und klarer Effektalleldefinition; alles andere als „explorativ“ abwerten. citeturn7search20turn16search3 Effektgrößen sind oft klein; isolierte Single-SNP-Aussagen sind meist schwach prädiktiv. citeturn3search12turn3search18

**Polygenic Risk Scores (PRS): Standardmethode und notwendige Metadaten**  
Standard-PRS ist ein gewichteter Summenscore: **PRS = Σ (β_i * G_i)**, wobei G_i die Anzahl Effektallele ist. citeturn3search3turn3search18 Für interpretierbare Perzentile braucht man eine Referenzverteilung (Mean/SD) im passenden Ancestry-Kollektiv; ohne Referenz ist der Score ein roher Index ohne klinischen Maßstab. citeturn3search3turn7search7

**Ancestry-Transferabilität: zentrale Limitation**  
GWAS/PRS-Übertragbarkeit ist populationsabhängig; viele GWAS stammen aus europäischen Kohorten, was zu schlechterer Performance in nicht-europäischen Gruppen führen kann. citeturn7search20turn3search18turn3search12

**Datenquellen/Implementierungsartefakte**  
```text
# GWAS Catalog: Top Hits Felddefinitionen (Risikoallel, OR/BETA)
https://www.ebi.ac.uk/gwas/docs/fileheaders

# GWAS Summary Statistics Standards / Harmonisierung
https://pmc.ncbi.nlm.nih.gov/articles/PMC11526975/      # GWAS-SSF in GWAS Catalog
https://www.sciencedirect.com/science/article/pii/S2666979X21000045  # Workshop/Standards
https://github.com/EBISPOT/gwas-sumstats-harmoniser

# PRS Best Practice
https://pmc.ncbi.nlm.nih.gov/articles/PMC7612115/
```
citeturn7search0turn7search8turn3search5turn7search15turn3search3

**Häufige Fallstricke**  
(1) Risikoallel nicht gegen Nutzerallel zählen (dein aktueller Zustand) → alle Treffer werden fälschlich als „relevant“ gemeldet. citeturn7search0turn18search7  
(2) OR aus Top-Hits ohne Beachtung der historischen Inversion/Allelswap interpretieren. citeturn16search0turn16search12  
(3) OR als „absolutes Risiko“ ausgeben ohne Baseline p0 und ohne RR/OR-Unterschied zu erklären. citeturn3search2turn3search10  
(4) GWAS-Effekte ohne Ancestry-Kontext anwenden. citeturn7search20turn3search18

**Priorität (Patientensicherheit)**  
Mittel bis hoch. Direkte medizinische Fehlsteuerung ist seltener als bei ClinVar/PGx, aber Falschinterpretation kann zu riskantem Verhalten und Fehlentscheidungen führen; PRS/Single-SNP sind in der klinischen Utility begrenzt. citeturn3search18turn11search19

## Populationsfrequenzen sinnvoll nutzen

**Wissenschaftlich korrekter Ansatz (ACMG/ClinGen-konform)**  
Populationsallelfrequenzen sind zentral, um Pathogenität plausibel zu machen oder zu widerlegen. ACMG/AMP enthält explizite Frequenzkriterien; BA1 ist „stand-alone benign“ bei hoher Frequenz, BS1 ist „benign strong“ bei Frequenz höher als für die Erkrankung erwartet, während PM2 „absent/very rare in controls“ als pathogenitätsstützend nutzt. citeturn2search2turn8search1turn8search5 Eine aktualisierte BA1-Empfehlung präzisiert, dass BA1 bei AF>0.05 (5%) in einem geeigneten Referenzdatensatz angewendet werden kann. citeturn8search1

**Welche Frequenzzahl ist die richtige? (global vs ancestry-matched, popmax)**  
Für Filtering und Plausibilitätschecks wird empfohlen, **popmax** zu verwenden (Maximum der Kontinentalpopulationen), weil eine Variante, die in einer Population häufig ist, i.d.R. nicht als hochpenetrant-mendeliansch krankheitsverursachend gelten kann. citeturn8search0 gnomAD selbst nutzt/erklärt Populationskategorien und stellt pro Population AF bereit. citeturn8search2turn8search6

**Wie Frequenzen in deine Pipeline gehören (konkret)**  
1) ClinVar P/LP nur dann als „hochrelevant“ klassifizieren, wenn popmax unter krankheitsspezifischen Schwellen liegt; bei sehr hoher Frequenz automatisch „Penetranz niedrig / Klassifikation fraglich / re-evaluate“ markieren. citeturn8search0turn8search5turn2search2  
2) Bei GWAS: Frequenzen sind notwendig, um Baseline-/Genotypverteilungen korrekt zu modellieren (insb. wenn absolute Risiken approximiert werden). citeturn7search20turn3search2  
3) Bei PGx: Frequenz-/Ancestry-Kontext ist relevant, weil relevante Allele je nach Population stark variieren; CPIC stellt hierfür Frequency Tables bereit. citeturn27search3turn28search24

**Ancestry-Inferenz lokal: vorsichtig, aber möglich**  
PCA-/Admixture-basierte Methoden sind Standard zur genetischen Ancestry-Beschreibung, aber es gibt bekannte Missbrauchs- und Fehlinterpretationsrisiken; insbesondere dürfen solche Tools nicht als historische/ethnische Aussagen missverstanden werden. citeturn8search7turn8search3 Für ein privacy-first Tool ist die technisch saubere Variante: rein statistische „gnomAD-superpopulation nearest“ Auswahl zur Frequenzanzeige, ohne identitätsnahe Labels, und immer mit „uncertain“ Option.

**Datenquellen/Implementierungsartefakte**  
```text
# ACMG/AMP + BA1 Update
https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/
https://pmc.ncbi.nlm.nih.gov/articles/PMC6188666/

# Popmax Best Practice
https://pmc.ncbi.nlm.nih.gov/articles/PMC9160216/

# gnomAD Population Labels / Ancestry
https://gnomad.broadinstitute.org/news/2017-02-the-genome-aggregation-database/
https://gnomad-sg.org/help/ancestry
```
citeturn2search2turn8search1turn8search0turn8search2turn8search6

**Häufige Fallstricke**  
(1) Frequenzen nur anzeigen, aber nicht in die Interpretation integrieren → ClinVar-Falschpositive bleiben ungefiltert. citeturn8search0turn2search2  
(2) Globale AF statt popmax/ancestry-matched verwenden → falsche Aussagen bei population-spezifischen Varianten. citeturn8search0turn8search6  
(3) BA1 stumpf als „>5% = benign“ implementieren ohne Krankheitskontext; ClinGen/ACMG betonen krankheitsspezifische Schwellen. citeturn8search5turn8search1

**Priorität (Patientensicherheit)**  
Hoch. Frequenzlogik ist ein primärer Schutz gegen „Pathogenic“-Überinterpretation von häufigen Varianten und against overcalling bei DTC-Rohdaten. citeturn8search1turn33search2

## Verantwortungsvolle Ergebnisdarstellung und Regulatorik

**Wissenschaftlich korrekter Ansatz (Risiko- und Compliance-getrieben)**  
Ein nicht-klinisches Interpretationstool muss strikt zwischen (a) „Informationsdarstellung“ und (b) „medizinischer Empfehlung“ trennen. Dies ist nicht nur ethisch, sondern regulatorisch relevant: Die FDA warnte öffentlich, dass viele PGx-Testclaims (inkl. direkt an Konsumenten vermarkteter Tests und Software-Interpretationen) nicht von der FDA geprüft sind und wissenschaftlich unzureichend gestützt sein können; Therapieänderungen auf Basis solcher Claims können Patientenschaden verursachen. citeturn30view0turn9search0

**FDA-Regelwerk / Labeling-Logik (relevant für DTC-nahe Reports)**  
21 CFR Part 809 beschreibt Labeling-Anforderungen für In-vitro-Diagnostika; 21 CFR 809.10 verlangt u.a. Angaben zu intended use, limitations und performance characteristics. citeturn9search1turn9search9 Eine FDA-Schulung 2024 zu 809.10(b) betont, dass Labeling konsistente Kerninformationen liefern soll (intended use, limitations/warnings, performance) und dass diese Elemente auch in Testreport-Templates abgebildet werden können. citeturn30view1turn9search5

**FDA und Pharmakogenetik: sichere Referenzpunkte**  
Die FDA führt eine „Table of Pharmacogenetic Associations“ (informativ, nicht automatisch ein DTC-Freifahrtschein) und wiederholt darin, dass Genotyping klinische Vigilanz und Patientenmanagement nicht ersetzt. citeturn19search2 Das ist als Template für Safety-Language in einem Tool nützlich: PGx ist kontextabhängig und nicht absolut.

**Disclaimers: wie etablierte Dienste die Grenze ziehen (Belege)**  
- Nebula formuliert explizit „informational/educational only“, „not intended for diagnostic purpose“, „no medical advice“. citeturn10search1turn10search5  
- SelfDecode beschreibt, dass das Produkt nicht zur Diagnose/Behandlung gedacht ist und keine medizinischen Entscheidungen daraus abgeleitet werden sollen. citeturn10search6  
Solche Disclaimers sind notwendig, aber nicht hinreichend: Ohne korrekte Allel- und Qualitätslogik bleibt der Output gefährlich.

**Warum „klinische Bestätigung“ zwingend in die UI gehört**  
Eine Studie zu klinischer Bestätigung von DTC-Rohdatenvarianten fand, dass ein großer Anteil der in DTC-Rohdaten berichteten Varianten in klinischer Bestätigung falsch-positiv war und dass einige als „increased risk“ markierte Varianten klinisch als benign klassifiziert wurden; die Autoren betonen die Notwendigkeit klinischer Bestätigungstests. citeturn33search2turn9search7 Das ist ein direktes Argument für UI-Flags: „Rohdaten sind nicht klinisch validiert; Bestätigung in einem qualifizierten Labor erforderlich.“

**Kommunikation an Nicht-Expert:innen: evidenzbasierte Report-Formate**  
Patientenfreundliche Reports profitieren von prominentem Ergebnis-Summary in verständlicher Sprache und klaren Interpretationsabschnitten; Fachliteratur empfiehlt explizit strukturierte, klare Darstellung, um Fehlinterpretation zu reduzieren. citeturn19search0turn19search17turn19search3 ClinGen stellt zudem Kommunikations-/Consent-Frameworks (z.B. CADRe) bereit, die Disclosure-Strategien strukturieren. citeturn19search12

**High-impact Findings separat flaggen (und warum)**  
Für sehr folgenreiche, medizinisch „actionable“ Gene/Krankheitsbilder existiert in der klinischen Genomik die ACMG Secondary Findings (SF) Policy mit kuratierten Listen (z.B. SF v3.2, v3.3) als Rahmen für verantwortungsvolle Rückgabe in klinischen Sequenzierungen. citeturn33search14turn33search1 Für ein Consumer-Tool folgt daraus eine Safety-Policy: **wenn** (und nur wenn) eine Variante nach (i) Allelmatch, (ii) hoher ClinVar-Review-Qualität, (iii) Frequenzplausibilität, und (iv) plausibler Vererbung/Genotyp-Konstellation als hochrelevant erscheint, muss sie als „High impact, confirm clinically“ ausgegeben werden, nicht als Diagnose. Die DTC-Falschpositiv-Daten stützen diese Härte. citeturn33search2turn8search1turn24view0

**Datenquellen/Implementierungsartefakte**  
```text
# FDA Warnings / Enforcement
https://www.fda.gov/media/125467/download
https://www.fda.gov/news-events/press-announcements/fda-issues-warning-letter-genomics-lab-illegally-marketing-genetic-test-claims-predict-patients

# 21 CFR Part 809
https://www.ecfr.gov/current/title-21/chapter-I/subchapter-H/part-809
https://www.ecfr.gov/current/title-21/chapter-I/subchapter-H/part-809/subpart-B/section-809.10

# FDA Pharmacogenetic Associations
https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations

# DTC false positives (clinical confirmation)
https://www.nature.com/articles/gim201838

# Patient-friendly report design
https://pmc.ncbi.nlm.nih.gov/articles/PMC4254435/
```
citeturn30view0turn9search0turn9search1turn9search9turn19search2turn33search2turn19search0

**Häufige Fallstricke**  
(1) Disclaimers als Ersatz für technische Korrektheit behandeln. citeturn30view0turn33search2  
(2) PGx als „automatische Dosierungsanweisung“ ohne Kontext präsentieren; FDA warnt vor nicht validierten Claims. citeturn30view0turn19search2  
(3) High-impact Varianten wie normale „Info“-Treffer darstellen statt als „confirm clinically“. citeturn33search2turn33search14

**Priorität (Patientensicherheit)**  
Kritisch. Dieser Bereich entscheidet, ob ein technisch korrektes Tool klinisch missbraucht wird.

## Referenzimplementierungen und Best Practices

**Wissenschaftlich korrekter Ansatz: aus Referenztools ableiten, nicht neu erfinden**  
Für PGx ist PharmCAT die Referenzimplementierung: PharmCAT akzeptiert VCF/„outside call“, identifiziert PGx-Genotypen, inferiert Star-Allele und erzeugt Reports mit Guideline-Empfehlungen (CPIC/DPWG); es besteht aus einem Preprocessor und einem Named Allele Matcher. citeturn1search1turn1search13turn1search5 Daraus ist die direkte Architekturableitung: (a) VCF-normalisieren, (b) named-allele matching datengetrieben aus Tabellen, (c) report generation strikt getrennt von matching.

Für generische Variant Annotation (ClinVar, gnomAD, eigene Tabellen) sind etablierte Annotatoren relevant:

- OpenCRAVAT ist modular, lokal betreibbar und pipeline-fähig. citeturn12search5turn12search1  
- Ensembl VEP dokumentiert sehr konkret, wie VCF-Einträge intern normalisiert/trimmed werden und dass Allelvergleich optional deaktiviert werden kann („don’t compare alleles“), was als Warnsignal dient: „colocated rsID“ ist nicht automatisch derselbe Allelzustand. citeturn12search7turn12search11  
- GA4GH VRS definiert Normalisierung als kanonische Repräsentation zur System-übergreifenden Vergleichbarkeit. citeturn12search0turn12search8

**GA4GH/VCF-Normalisierung: was zwingend ist und wann**  
Variation Normalization (VRS) zielt auf kanonische Formen, um „äquivalente“ Varianten eindeutig zu machen. citeturn12search0turn12search4 Für VCFs sind klassische Schritte: (i) Multi-allelic split, (ii) REF-Check gegen FASTA, (iii) left-alignment von Indels in repetitiven Regionen, (iv) trimming gemeinsamer Basen; bcftools norm ist der de-facto Standard-Utility für solche Operationen. citeturn12search2turn12search14turn12search7

Für reine SNP-Array-Genotypen sind left-alignment/trimming selten relevant (weil meist SNVs), aber spätestens wenn du aus Arraydaten VCF erzeugst oder echte VCFs ingestierst, muss Normalisierung Teil der Standardpipeline sein. citeturn1search1turn12search2turn12search0

**ClinVar-Daten: VCF vs Full Release**  
ClinVar-FTP primer stellt klar, dass ClinVar-VCF nur bestimmte Variantentypen (simple alleles, <10kb, präzise) umfasst; eine „vollständige“ ClinVar-Abdeckung kann Full Release (XML/TXT) erfordern. citeturn13search9turn13search17 Für ein lokales Consumer-Tool ist VCF oft ausreichend, aber die Limitations müssen UI-seitig sichtbar sein („coverage limitations“).

**Best-Practice-Pipeline (kompakt, implementierbar)**  
```pseudo
ingest(file):
  detect_provider_and_build(file)             # DTC header / heuristics
  parse_genotypes(file)
  normalize_to_internal_assembly(GRCh37|38)   # optional liftover + unmappable tracking
  resolve_strand_and_alleles()                # plus-strand canonical
  represent_variants_as (chr,pos,REF,ALT)     # canonical key, rsID as secondary

annotate_clinvar():
  load_clinvar_vcf_for_assembly()
  allele_match_then_interpret_with_review_status()

annotate_gwas():
  prefer_summary_stats_harmonised()
  allele_match_then_effect_allele_count()
  compute_OR_or_beta_effects()
  (optional) compute_PRS_with_reference_distribution()

annotate_pgx():
  build_vcf_subset_for_pgx_sites
  run_named_allele_matching_datengetrieben (PharmCAT-like)
  output_diplotype + phenotype + guideline pointers
  carry_forward "missing sites" + "CNV unknown" flags
```

**Datenquellen/Implementierungsartefakte**  
```text
# PharmCAT
https://pharmcat.clinpgx.org/using/
https://github.com/PharmGKB/PharmCAT-tutorial
https://pmc.ncbi.nlm.nih.gov/articles/PMC10121724/

# OpenCRAVAT
https://docs.opencravat.org/

# Ensembl VEP (VCF handling)
https://www.ensembl.org/info/docs/tools/vep/vep_formats.html

# GA4GH VRS Normalization
https://vrs.ga4gh.org/en/1.2/impl-guide/normalization.html

# bcftools norm
https://samtools.github.io/bcftools/bcftools.html
```
citeturn1search1turn1search5turn1search13turn12search5turn12search7turn12search0turn12search2

**Häufige Fallstricke**  
(1) Annotator-Output (rsID co-location) als identische Allele behandeln; VEP dokumentiert explizit, dass Allelvergleich eine eigene Entscheidung ist, nicht implizit sicher. citeturn12search11turn12search7  
(2) PGx als „regelbasiertes Hardcoding“ statt datengetriebener Tabellen (CPIC xlsx) implementieren → Update-Drift. citeturn27search0turn28search24  
(3) ClinVar-VCF als vollständig ansehen → Coverage-Lücken. citeturn13search9turn13search17  
(4) Normalisierung/REF-Check überspringen → stille Mismatches bei Indels/Repeats. citeturn12search2turn12search0

**Priorität (Patientensicherheit)**  
Hoch. Referenzarchitekturen (PharmCAT, VEP/VRS-Normalisierung) liefern die sichersten Blaupausen, um die Hauptklassen von Silent-Failure (Allel-/Strang-/Repräsentationsmismatch) zu vermeiden. citeturn1search13turn12search0turn12search11