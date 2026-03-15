# Scientifically Correct Implementation of a Consumer Genomics Annotation Pipeline

Target state: Every annotation is produced only after **canonical variant representation** (build + coordinates + REF/ALT + strand) and **allele-to-allele matching**; rsID is merely an alias and never the primary key for clinical assertions. The safety-critical bugs in the current pipeline (rsID lookup without allele comparison) systematically produce false assignments for pathogenicity, risk alleles, and pharmacogenetics.

## Allele Matching and Strand Orientation

**Scientifically correct approach (algorithm + precise rules)**
Consumer raw data are typically already provided as nucleotides (A/C/G/T); the critical factor is **which reference (build) and which strand** these nucleotides refer to. For two major DTC formats: 23andMe reports genotypes by default on the **plus strand** relative to **GRCh37** (and also offers GRCh38 on the plus strand in the raw data browser), and AncestryDNA reports raw data on the **forward/plus strand** relative to **GRCh37**. citeturn14search0turn14search1 This reduces but does not eliminate strand issues, because all external reference databases (and/or historical exports) may be differently encoded or differently normalized.

At the array/manifest level, additional strand definitions exist: Illumina uses, among others, **TOP/BOT** (stable across build changes, not identical to +/-) as well as manifest columns such as **IlmnStrand** (TOP/BOT for SNPs and PLUS/MINUS for indels) and **SourceStrand** (submitted customer/database strand designation). citeturn0search10turn18search0turn21view3 If your pipeline supports not only DTC exports but also (i) generic Illumina final reports or (ii) manifest-based formats, it must explicitly resolve these strand orientations.

dbSNP is suitable as the primary source for rsID→(chr,pos,REF,ALT), but history matters: NCBI describes that RefSNPs were formerly sometimes reported as FWD/REV relative to the assembly and that in the newer dbSNP architecture, alleles are consistently reported "forward to the reported sequence" (VCF/HGVS/SPDI-compliant); therefore, special VCFs also exist for formerly "REV"-reported rsIDs. citeturn21view4turn13search8 For robust allele comparisons, the correct rule is: **Normalize alleles to the same build and the same (plus) orientation and only then compare.**

ClinVar delivers variants as VCF via FTP, among other formats; these VCFs contain "simple alleles" (<10kb) with precise endpoints, mapped to GRCh37 or GRCh38 (not all ClinVar variant types are included in the VCF). citeturn13search9turn13search17 Every ClinVar assertion you display must therefore (a) match the correct assembly and (b) compare your user alleles with the alleles referenced in ClinVar (REF/ALT or CLNSIG-related alleles).

The GWAS Catalog "Top Hits" files contain a field "STRONGEST SNP-RISK ALLELE" (SNP + risk/effect allele; "?" if unknown) as well as a field "OR or BETA" for the effect; for older curations, OR<1 was sometimes inverted and the reported allele flipped accordingly, so that stored ORs are >1. citeturn7search0turn16search0turn16search12 This means: You must not blindly adopt risk alleles, but instead must check the respective catalog logic (date/version) and allele orientation, or switch to harmonized summary stats.

image_group{"layout":"carousel","aspect_ratio":"16:9","query":["Illumina TOP BOT strand diagram","Illumina Infinium manifest IlmnStrand SourceStrand explanation","palindromic SNP A/T C/G strand ambiguity diagram","SNP strand flip plus minus complement A T C G diagram"],"num_per_query":1}

**Strand/allele normalization: reference algorithm (pseudocode)**
Core principle: The comparison is performed on **(chr,pos,REF,ALT)** in a defined assembly. rsID is only a lookup aid.

```pseudo
# Inputs:
#   user: (assembly, chr, pos, genotype="AG"/"--"/"A", provider_hint)
#   refdb: (assembly, chr, pos, REF, ALTs[])  # verified from dbSNP/VCF/FASTA
#   target_allele: a nucleotide (A/C/G/T) or for indels a sequence (VCF-ALT)

function complement_base(b):
    map = { "A":"T", "T":"A", "C":"G", "G":"C" }
    return map[b]

function is_palindromic_snp(REF, ALT):
    # biallelic SNP is palindromic when {A,T} or {C,G}
    s = set([REF, ALT])
    return (s == {"A","T"} or s == {"C","G"})

function normalize_user_alleles_to_plus(user, refdb):
    # 1) Build/coordinate must already match refdb assembly (liftover beforehand if needed)
    assert user.assembly == refdb.assembly
    allowed = set([refdb.REF] + refdb.ALTs)

    alleles = parse_genotype(user.genotype)   # "AG"->["A","G"], "--"->[None,None]
    if any(a not in {"A","C","G","T",None} for a in alleles):
        return {status:"unsupported_allele_encoding"}

    if None in alleles:
        return {status:"missing_genotype"}

    # 2) Direct match?
    if set(alleles) ⊆ allowed:
        return {status:"ok", oriented:"+", alleles:alleles}

    # 3) Reverse-complement match?
    comp = [complement_base(a) for a in alleles]
    if set(comp) ⊆ allowed:
        # Ambiguity for palindromic SNPs: Flip not provable without additional info
        if len(refdb.ALTs)==1 and is_palindromic_snp(refdb.REF, refdb.ALTs[0]):
            return {status:"ambiguous_palindromic_strand"}
        return {status:"ok", oriented:"-", alleles:comp}

    return {status:"allele_mismatch"}  # wrong coordinate, wrong build, multi-allelic, rsID merge, etc.

function count_target_allele(alleles, target_allele):
    # 0/1/2 copies for biallelic SNPs; analogous for multi-allelic
    return sum(1 for a in alleles if a == target_allele)
```

This schema is minimal; in practice, additional considerations include: multi-allele rsIDs, indels, and the need to **verify REF against a reference FASTA** (or to take REF/ALT directly from a reliable VCF). citeturn12search2turn12search14turn13search8

**Palindromic SNPs (A/T, C/G): safe handling**
For A/T or C/G SNPs, "flip vs. no flip" is undecidable based on the allele set alone. This is the classic error domain in GWAS/array data harmonization. Tools such as **snpflip** flag reverse/ambiguous SNPs, and **Genotype Harmonizer** can align ambiguous A/T and G/C SNPs via LD patterns against a reference panel (in a DTC local tool usually too heavyweight; but important as a concept). citeturn18search1turn18search7 For chip-specific resolution, curated "strand and build files" also exist for many genotyping chips. citeturn18search22

**Data sources/implementation artifacts (concrete, machine-readable)**
```text
# DTC Strand/Build
https://eu.customercare.23andme.com/hc/en-us/articles/115002090907-Raw-Genotype-Data-Technical-Details
https://support.ancestry.com/articles/en_GB/Support_Site/Downloading-DNA-Data

# Illumina Manifest/Strand Definition
https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565
https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001489
https://www.illumina.com/documents/products/technotes/technote_topbot.pdf

# dbSNP / NCBI Variation Services (rsID↔SPDI/VCF)
https://www.ncbi.nlm.nih.gov/core/assets/snp/docs/RefSNP_orientation_updates.pdf
https://github.com/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/spdi_batch.py
https://pmc.ncbi.nlm.nih.gov/articles/PMC7523648/   # SPDI paper

# ClinVar FTP Primer (VCF scope, GRCh37/38)
https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/

# Strand Resolution Tools / Chip Strand Files
https://github.com/andymckenzie/snpflip
https://pmc.ncbi.nlm.nih.gov/articles/PMC4307387/   # Genotype Harmonizer
https://www.well.ox.ac.uk/~wrayner/strand/
```
citeturn14search0turn14search1turn0search2turn18search0turn21view3turn21view4turn13search0turn13search8turn13search9turn18search1turn18search7turn18search22

**Common pitfalls (typical for consumer genomics)**
(1) rsID-based display of "Pathogenic/Risk" without allele comparison produces guaranteed false positives (your current bug is the standard pattern). citeturn24view2turn21view4
(2) Build mix (GRCh37 vs GRCh38) causes apparent "allele mismatch" and incorrect annotations; liftover often works very well, but can fail for certain variant types/regions, so "unmappable" must be explicitly handled. citeturn18search13turn13search8
(3) Palindromic SNPs flipped without additional info → silent corruption of risk alleles. citeturn18search7turn18search1
(4) Illumina A/B alleles or TOP/BOT interpreted as A/C/G/T without validation. citeturn21view3turn18search0
(5) Using ClinVar VCF as "complete ClinVar truth" (it is not; parts exist only in the XML/TXT full release). citeturn13search9turn13search17

**Priority (patient safety)**
Critical. Without correct strand/allele matching, all subsequent categories (ClinVar, GWAS, PGx) are potentially systematically wrong.

## PGx Star Allele Calling and Phenotyping

**Scientifically correct approach (CPIC/PharmVar-compliant)**
Pharmacogenomic phenotypes are generally not "one SNP → one phenotype", but rather **haplotypes (star alleles) → diplotype → function/activity score → phenotype → drug guidance**. CPIC explicitly describes that the combination of alleles determines the diplotype (genotype) and that phenotype classes are derived from it; for CYP2C19, for example, *1/*17 = Rapid and *17/*17 = Ultrarapid, while *2/*17 is classified as Intermediate despite *17. citeturn21view0turn27search2

**Why your current approach is inevitably wrong:** rs12248560 is a definer for CYP2C19*17; without allele comparison and without diplotype logic, *1/*1 is falsely labeled as "Ultrarapid". PharmVar shows *17 as -806C>T (rs12248560) and lists the associated core variants. citeturn1search2turn21view0

**Star allele definitions: robust data sources instead of hardcoding**
Hardcoding "a few known rsIDs" is not maintainable and quickly becomes incorrect when star nomenclature is updated (star definitions are dynamic). citeturn15search17turn26search3 CPIC provides machine-readable tables for this purpose: **allele_definition**, **allele_functionality_reference**, **diplotype_phenotype**, **frequency**, **gene_cds**. citeturn27search4turn28search15 This is the correct path for a local tool: version the tables, ship them offline, update transparently.

**Concrete, commonly used definer variants (examples, not exhaustive)**
The following examples serve as an entry point/validation; in production, the definition should always come from the CPIC/PharmVar tables:

- CYP2C19: *17 rs12248560 (Promoter -806C>T), *2 rs4244285, *3 rs4986893. citeturn1search2turn21view0turn27search13
- CYP2C9: *2 rs1799853, *3 rs1057910. citeturn25search0turn25search5
- CYP3A5: *3 rs776746 (splice defect). citeturn25search6turn25search2
- DPYD: *2A rs3918290 (splice defect); CPIC works with c.1905+1G>A, c.1679T>G, c.2846A>T, c.1129–5923C>G among others as central variants for DPD activity estimation. citeturn25search7turn23view1
- SLCO1B1: *5 rs4149056 (c.521T>C; V174A). citeturn26search7turn23view2turn26search1
- TPMT: *3A (rs1800460 + rs1142345 in cis), *3B rs1800460, *3C rs1142345, *2 rs1800462; *3A is potentially ambiguous without phasing. citeturn26search6turn26search0turn26search10turn26search4
- NUDT15: *3 rs116855232 (R139C). citeturn25search4turn25search20
- CYP2D6: very many alleles (PharmVar lists >130 core alleles); individual key variants (e.g., *4 rs3892097, *10 rs1065852) represent only a small part, structural variants (deletion/duplication/hybrid) are clinically relevant. citeturn26search3turn26search9turn21view2

**Correct end-to-end logic (SNP genotype → star alleles → diplotype → activity score → phenotype)**
For a local consumer tool with array/VCF input, the robust implementation is:

1) **Harmonize input**: Bring variants to one assembly (GRCh37 or GRCh38), normalize alleles to plus strand, verify REF against reference. citeturn14search0turn14search1turn12search2turn13search8
2) **Load allele definitions**: CPIC allele_definition_table for the gene (star alleles as variant combinations). citeturn27search0turn28search11
3) **Named allele matching**: For each star definition, check whether the user genotype data permit this haplotype definition (with missingness tracking). PharmCAT does this as a "Named Allele Matcher" on VCF input. citeturn1search5turn1search1turn1search13
4) **Diplotype inference without phasing**: Form the possibility space of all haplotype pairs that explain the (unphased) genotype counts; in case of ambiguity, return "ambiguous call" instead of a forced result. citeturn26search0turn25search1
5) **Function + activity score**: Use CPIC allele_functionality_reference + diplotype_phenotype tables. For CYP2D6, an activity value is assigned per allele, multiplied for duplications, then summed (activity score). citeturn21view2turn28search1turn28search0turn28search2
6) **Phenotype & drug guidance**: Gene-specific CPIC recommendation logic/tables; for many genes, CPIC "gene_cds" tables exist as structured decision support. citeturn27search6turn28search7turn23view1turn23view2turn31search0

**Pseudocode (diplotype + phenotype, data-driven)**
```pseudo
# data-driven: do not hardcode, instead load CPIC tables

function call_gene_pgx(sample_vcf, gene):
    defs  = load_cpic_allele_definitions(gene)        # allele_definition_table.xlsx
    func  = load_cpic_allele_function(gene)           # allele_functionality_reference.xlsx
    d2p   = load_cpic_diplotype_to_phenotype(gene)    # Diplotype_Phenotype_Table.xlsx

    # 1) extract relevant variants
    g = subset_variants(sample_vcf, defs.all_required_sites)

    # 2) determine possible haplotypes (compatible with observed genotypes)
    possible_haplotypes = []
    for star in defs.star_alleles:
        if star_is_compatible(star, g):   # accounts for missingness + conflicts
            possible_haplotypes.append(star)

    # 3) enumerate diplotypes (unphased)
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
    phenotype = d2p.lookup(diplotype)  # preferred: CPIC table, not self-computed
    return {status:"ok", diplotype:diplotype, phenotype:phenotype}
```

**Missing SNPs: safe default rules**
"Not tested" is not "wildtype". CPIC points out in multiple guidelines that genotype-only tests do not capture rare/newly discovered variants; consequently: if defining sites are missing, "Normal" is often not derivable and must be returned as "Unclear/Indeterminate". citeturn21view2turn23view0turn28search0 In practice: Output a **coverage metric** for each gene ("X% of defining sites for common core alleles observed"), and only output a phenotype if the CPIC tables allow the call as determined. citeturn27search2turn28search2

**CNVs / structural variants (especially CYP2D6): hard limits of consumer arrays**
CYP2D6 phenotype depends heavily on deletions and duplications ("xN"); CPIC describes deletion (*5) and duplications explicitly and that activity values are multiplied for multi-copy alleles. citeturn21view2turn28search1 Many consumer arrays do not reliably detect these CNVs; therefore, your tool must carry a **"CNV unknown" status** for CYP2D6 by default and must not claim Ultrarapid calls when no copy number information is available. Tools like Stargazer can in principle also work with SNP array data, but CNV resolution remains limited depending on input; Aldy/Astrolabe are primarily designed for sequence data and structural models. citeturn15search4turn15search3turn15search2

**CPIC recommendations: highly relevant examples (reference data-driven, not paraphrased)**
- CYP2C19–Clopidogrel: CPIC 2022 classifies *1/*17 as Rapid, *17/*17 as Ultrarapid, and *2/*17 as Intermediate; treatment recommendations are derived from this. citeturn21view0turn20view0
- DPYD–Fluoropyrimidines: CPIC uses the DPYD Activity Score; Table 2 provides dose-reduced starting doses for partial DPD deficiency (e.g., 50% reduction for certain AS constellations) and emphasizes titration/monitoring. citeturn23view1
- TPMT/NUDT15–Thiopurines: Table 1 maps diplotypes to phenotypes, recommendations reduce starting doses depending on combination status and allow dose adjustment based on myelosuppression. citeturn23view0
- SLCO1B1–Simvastatin: rs4149056 genotype is translated into phenotype classes; recommendations avoid 80 mg and suggest lower doses or alternative statins, especially for C alleles. citeturn23view2turn17search3
- CYP2D6–Codeine/Tramadol: CPIC recommends not using codeine/tramadol in Ultrarapid metabolizers (toxicity risk) and alternatives in Poor metabolizers; the 2021 opioid guideline formulates this with activity score thresholds. citeturn31search0turn31search4

**Data sources/implementation artifacts (CPIC "Full Tables", machine-readable)**
```text
# CPIC: CYP2C19 (Definition/Function/Diplotype→Phenotype)
https://files.cpicpgx.org/data/report/current/allele_definition/CYP2C19_allele_definition_table.xlsx
https://files.cpicpgx.org/data/report/current/allele_function_reference/CYP2C19_allele_functionality_reference.xlsx
https://files.cpicpgx.org/data/report/current/diplotype_phenotype/CYP2C19_Diplotype_Phenotype_Table.xlsx
https://cpicpgx.org/gene/cyp2c19/

# CPIC: CYP2D6 (Definition/Function/Diplotype→Phenotype)
https://files.cpicpgx.org/data/report/current/allele_definition/CYP2D6_allele_definition_table.xlsx
https://files.cpicpgx.org/data/report/current/allele_function_reference/CYP2D6_allele_functionality_reference.xlsx
https://files.cpicpgx.org/data/report/current/diplotype_phenotype/CYP2D6_Diplotype_Phenotype_Table.xlsx
https://files.cpicpgx.org/data/report/current/gene_phenotype/CYP2D6_phenotypes.xlsx

# PharmCAT (Reference Implementation)
https://pharmcat.clinpgx.org/using/
https://github.com/PharmGKB/PharmCAT-tutorial
https://pmc.ncbi.nlm.nih.gov/articles/PMC10121724/
```
citeturn27search0turn27search1turn27search2turn27search4turn28search11turn28search1turn28search0turn28search2turn1search1turn1search5turn1search13

**Common pitfalls**
(1) Interpreting "absent in data" as *1 → false Normal/Rapid/Ultrarapid calls. citeturn23view0turn21view2
(2) Forcing TPMT*3A without phasing → incorrect diplotype (cis/trans). citeturn26search0turn26search4
(3) Reporting CYP2D6 as complete without CNV status → high-risk misassessments. citeturn21view2turn15search9
(4) Ignoring star definition updates → drift between report and current CPIC/PharmVar nomenclature. citeturn15search17turn28search24

**Priority (patient safety)**
High to critical. Incorrect PGx phenotypes can lead to real medication mismanagement; the FDA explicitly warns against non-validated PGx claims and also software interpretations in the DTC context. citeturn30view0turn9search0

## Correctly Interpreting ClinVar Pathogenicity

**Scientifically correct approach (not "contains pathogenic")**
ClinVar is an archive of submitted interpretations (SCVs), aggregated into variant/variant-condition records (VCV/RCV). citeturn32search12turn24view2 Since 2024, ClinVar separates clinical classification types (germline, somatic clinical impact, oncogenicity) into separate fields; "clinical_significance" must therefore be contextualized. citeturn24view2

The correct logic has three stages:

1) **Check alleles**: Only if the user carries the ALT allele(s) classified in ClinVar is the ClinVar classification genotype-relevant at all (otherwise "homozygous reference"/"non-carrier"). citeturn13search9turn21view4
2) **Interpret classification type + term**: ClinVar terms include, among others, Pathogenic/Likely pathogenic/Benign/Likely benign/VUS as well as "risk factor", "drug response", "association/protective/other". citeturn24view1turn24view2
3) **Weight evidence quality (review status/stars)**: ClinVar defines stars/review status: 4 = practice guideline, 3 = expert panel, 2 = multiple submitters/no conflicts, 1 = criteria provided (single submitter) or criteria provided (conflicting), 0 = no criteria/no classification, etc. citeturn24view0

**Reliability thresholds (safety-oriented, ClinVar-compliant)**
For a consumer tool, a conservative policy is necessary:

- **High confidence**: 3-4 stars (Expert Panel/Practice Guideline). citeturn24view0turn32search5
- **Medium**: 2 stars (multiple submitters, criteria, no conflict). citeturn24view0
- **Low/Informational**: 1 star (single submitter with criteria) -- display, but clearly label as non-consensus. citeturn24view0
- **Conflict case**: 1 star "criteria provided, conflicting classifications" -- never present as unambiguous pathogenicity; instead display the conflict structure (who says what, which review status, which date). citeturn24view0turn24view2
- **0 stars / no criteria**: do not present as a clinical assertion; at most as a raw hint. citeturn24view0

**Zygosity + inheritance: correct, programmable integration**
Genotype (0/1/2 ALT copies) is interpretable for monogenic diseases only in the context of **inheritance (AD/AR/X-linked/mitochondrial)**. Fundamentals of inheritance modes are standardized (autosomal dominant/recessive, X-linked, etc.). citeturn32search15 Programmable sources for "Mode of Inheritance" (MOI):

- ClinVar can contain MOI at the submission level; when a submitter specifies MOI, it is displayed on variant pages. citeturn32search1turn24view2
- ClinVar "properties" and filter terms contain MOI categories (moi autosomal dominant/recessive, etc.). citeturn32search8
- MedGen supports "mode of inheritance" as a search/property field and processes, among others, Orphanet/ORDO including MOI. citeturn32search0turn32search4
- ClinGen Gene-Disease Validity Knowledge Base maintains "Mode of Inheritance" per gene-disease assertion and is publicly accessible. citeturn2search19turn2search7

**Correct ClinVar interpretation logic (pseudocode, low-risk)**
```pseudo
# Inputs:
#   variant_call: {user_count_alt, genotype_quality, allele_normalized_ok}
#   clinvar_records: list of RCV-like entries {condition_id, clinsig_term, review_status, last_eval_date, moi?}
# Output:
#   structured interpretation objects (not "diagnosis", but variant hint)

function interpret_clinvar(variant_call, clinvar_records):
    if variant_call.allele_normalized_ok != true:
        return {status:"cannot_interpret_without_allele_match"}

    if variant_call.user_count_alt == 0:
        return {status:"non_carrier"}  # no ALT/risk allele copy

    # filter for germline classification (not somatic clinical impact / oncogenicity)
    records = filter_germline_records(clinvar_records)

    # weighted selection by review status
    # 4>3>2>1>0; conflict status separate
    best = pick_by_highest_review_status(records)

    # If best is "conflicting classifications": never decide binarily
    if best.review_status_contains("conflicting"):
        return {status:"conflicting", details: summarize_conflict(records)}

    # correctly handle clinsig terminology
    if best.clinsig in {"Pathogenic","Likely pathogenic"}:
        return {status:"P_or_LP", zygosity:variant_call.user_count_alt, moi:best.moi, caveats:...}
    if best.clinsig in {"Benign","Likely benign"}:
        return {status:"B_or_LB", ...}
    if best.clinsig in {"Uncertain significance"}:
        return {status:"VUS", ...}
    if best.clinsig in {"risk factor","association","protective","drug response","other"}:
        return {status:"non_mendelian_or_pgx_term", ...}
```

**ACMG/AMP framework: what it is and how it relates to ClinVar**
ACMG/AMP (Richards et al., 2015) defines the standard terminology (pathogenic/likely pathogenic/VUS/likely benign/benign) and evidence-based criteria including population frequency rules (BA1/BS1/PM2, etc.). citeturn2search2turn8search1turn8search5 ClinVar terms are aligned with this terminology, but ClinVar is not itself an ACMG classifier; it aggregates submitter assertions and weights aggregation based on review status, among other factors. citeturn24view2turn24view0 The ClinVar star system is therefore **evidence/process quality**, not a formal mapping to ACMG criteria fulfillment. citeturn24view0turn2search2

**Data sources/implementation artifacts**
```text
# ClinVar Docs: Review Status / Stars, Classification Types, Terms
https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
https://www.ncbi.nlm.nih.gov/clinvar/docs/properties/
https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/

# ACMG/AMP 2015 (Primary Source)
https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/

# MedGen MOI
https://www.ncbi.nlm.nih.gov/medgen/docs/search/
https://www.ncbi.nlm.nih.gov/medgen/docs/data/

# ClinGen Gene-Disease Validity KB (MOI visible)
https://search.clinicalgenome.org/
```
citeturn24view0turn24view2turn24view1turn13search9turn2search2turn32search0turn32search4turn2search19

**Common pitfalls**
(1) Displaying "Pathogenic" despite 0 copies of ALT; this is exactly your rsID-only bug. citeturn13search9turn24view1
(2) Treating "Conflicting interpretations" as "Pathogenic"; ClinVar explicitly distinguishes conflicts in the review status. citeturn24view0turn24view1
(3) Presenting somatic/oncogenicity classifications as germline "Pathogenic"; ClinVar separates these classification types. citeturn24view2
(4) Ignoring MOI/zygosity → AR diseases falsely reported as affected in heterozygous carriers. citeturn32search15turn32search0

**Priority (patient safety)**
Critical. Mendelian "Pathogenic" labels without allele and review status controls are the fastest path to highly harmful false alarms.

## GWAS Risk Interpretation and PRS

**Scientifically correct approach (single-variant and multi-variant)**
GWAS top hits are associations with typically small effects; correct user interpretation requires: (1) effect/risk allele, (2) effect size (OR or beta) referenced to that allele, (3) genotype (0/1/2 copies), (4) population frequencies where applicable, and (5) baseline prevalence for absolute risks. citeturn7search0turn7search20turn37search1

**Which alleles in the GWAS Catalog are "risk" alleles (and why this is tricky)**
In the curated GWAS Catalog top hits, "STRONGEST SNP-RISK ALLELE" denotes the variant plus risk/effect allele; "?" if unknown. citeturn7search0turn7search2 Additionally, the effect column "OR or BETA" is context-dependent; for studies curated before January 2021, OR<1 was sometimes inverted and the reported allele flipped accordingly, so that OR>1 is stored. citeturn16search0turn16search12 Consequently: The "risk allele" in top hits is not automatically a consistently plus-strand-harmonized effect allele in the VCF sense.

For robust, scalable interpretation, harmonized summary statistics (GWAS-SSF) are preferred, as they explicitly capture **effect_allele/other_allele, effect size, SE, p**, etc. and can be harmonized through pipelines. citeturn7search8turn3search5turn7search5

**Genotype → risk allele copies (0/1/2): correct logic**
After strand/build normalization, the GWAS side is analogous to ClinVar:

- 0 copies of effect allele: reference for the (study-)defined baseline comparison.
- 1 copy: heterozygous.
- 2 copies: homozygous.

This counting is only valid if the effect allele and user alleles are on the same reference orientation. citeturn18search7turn7search15turn7search5

**Correctly translating OR into personal (relative/absolute) risk measures**
OR is an odds ratio, not directly "risk". For an additive log-odds assumption per effect allele, the typical formula is: **Odds_multiplier = OR^k** (k = 0/1/2 effect alleles). citeturn7search20turn3search17 For an absolute risk approximation, you need a baseline risk assumption p0 (e.g., prevalence/lifetime risk in the target population): odds0 = p0/(1-p0); odds = odds0 * OR^k; p = odds/(1+odds). The relationship OR↔RR depends strongly on p0; at high baseline rates, OR and RR diverge substantially. citeturn3search2turn3search10turn3search17 Without p0, "absolute risk" cannot be seriously calculated.

**Beta (quantitative traits) treated differently from OR (binary)**
Beta is an additive effect per effect allele on the trait scale (or a transformed scale), so the naive expected shift is **ΔTrait ≈ beta * k**; additionally, one needs unit/scaling and population reference (mean/SD) to enable patient understanding. citeturn3search5turn7search20

**Reporting thresholds: what is scientifically defensible (consumer context)**
GWAS literature emphasizes the need for standardized reporting components (including alleles/strand/effect sizes) and that minimum information must be present. citeturn7search20turn3search5 For consumer reporting, a defensible minimum filter is: only associations with genome-wide significance (classically p<5x10^-8) and clear effect allele definition; everything else should be downgraded to "exploratory". citeturn7search20turn16search3 Effect sizes are often small; isolated single-SNP statements are usually weakly predictive. citeturn3search12turn3search18

**Polygenic Risk Scores (PRS): standard method and required metadata**
Standard PRS is a weighted sum score: **PRS = Σ (β_i * G_i)**, where G_i is the number of effect alleles. citeturn3search3turn3search18 For interpretable percentiles, a reference distribution (mean/SD) in the matching ancestry cohort is needed; without a reference, the score is a raw index without a clinical benchmark. citeturn3search3turn7search7

**Ancestry transferability: central limitation**
GWAS/PRS transferability is population-dependent; many GWAS originate from European cohorts, which can lead to poorer performance in non-European groups. citeturn7search20turn3search18turn3search12

**Data sources/implementation artifacts**
```text
# GWAS Catalog: Top Hits Field Definitions (Risk Allele, OR/BETA)
https://www.ebi.ac.uk/gwas/docs/fileheaders

# GWAS Summary Statistics Standards / Harmonization
https://pmc.ncbi.nlm.nih.gov/articles/PMC11526975/      # GWAS-SSF in GWAS Catalog
https://www.sciencedirect.com/science/article/pii/S2666979X21000045  # Workshop/Standards
https://github.com/EBISPOT/gwas-sumstats-harmoniser

# PRS Best Practice
https://pmc.ncbi.nlm.nih.gov/articles/PMC7612115/
```
citeturn7search0turn7search8turn3search5turn7search15turn3search3

**Common pitfalls**
(1) Not counting the risk allele against the user allele (your current state) → all hits are falsely reported as "relevant". citeturn7search0turn18search7
(2) Interpreting OR from top hits without considering historical inversion/allele swap. citeturn16search0turn16search12
(3) Presenting OR as "absolute risk" without baseline p0 and without explaining the RR/OR difference. citeturn3search2turn3search10
(4) Applying GWAS effects without ancestry context. citeturn7search20turn3search18

**Priority (patient safety)**
Medium to high. Direct medical mismanagement is less common than with ClinVar/PGx, but misinterpretation can lead to risky behavior and poor decisions; PRS/single-SNP are limited in clinical utility. citeturn3search18turn11search19

## Using Population Frequencies Meaningfully

**Scientifically correct approach (ACMG/ClinGen-compliant)**
Population allele frequencies are central to making pathogenicity plausible or refuting it. ACMG/AMP contains explicit frequency criteria; BA1 is "stand-alone benign" at high frequency, BS1 is "benign strong" at frequency higher than expected for the disease, while PM2 uses "absent/very rare in controls" as supporting pathogenicity. citeturn2search2turn8search1turn8search5 An updated BA1 recommendation specifies that BA1 can be applied at AF>0.05 (5%) in an appropriate reference dataset. citeturn8search1

**Which frequency number is the right one? (global vs ancestry-matched, popmax)**
For filtering and plausibility checks, it is recommended to use **popmax** (maximum across continental populations), because a variant that is common in one population generally cannot be considered highly penetrant Mendelian disease-causing. citeturn8search0 gnomAD itself uses/explains population categories and provides AF per population. citeturn8search2turn8search6

**How frequencies belong in your pipeline (concrete)**
1) Only classify ClinVar P/LP as "highly relevant" when popmax is below disease-specific thresholds; at very high frequency, automatically mark as "penetrance low / classification questionable / re-evaluate". citeturn8search0turn8search5turn2search2
2) For GWAS: frequencies are necessary to correctly model baseline/genotype distributions (especially when approximating absolute risks). citeturn7search20turn3search2
3) For PGx: frequency/ancestry context is relevant because relevant alleles vary strongly by population; CPIC provides frequency tables for this purpose. citeturn27search3turn28search24

**Local ancestry inference: cautious, but possible**
PCA/admixture-based methods are standard for genetic ancestry description, but there are known risks of misuse and misinterpretation; in particular, such tools must not be misunderstood as historical/ethnic statements. citeturn8search7turn8search3 For a privacy-first tool, the technically clean approach is: purely statistical "gnomAD superpopulation nearest" selection for frequency display, without identity-adjacent labels, and always with an "uncertain" option.

**Data sources/implementation artifacts**
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
citeturn2search2turn8search1turn8search0turn8search2turn8search6

**Common pitfalls**
(1) Only displaying frequencies but not integrating them into the interpretation → ClinVar false positives remain unfiltered. citeturn8search0turn2search2
(2) Using global AF instead of popmax/ancestry-matched → incorrect assertions for population-specific variants. citeturn8search0turn8search6
(3) Implementing BA1 bluntly as ">5% = benign" without disease context; ClinGen/ACMG emphasize disease-specific thresholds. citeturn8search5turn8search1

**Priority (patient safety)**
High. Frequency logic is a primary defense against "Pathogenic" overinterpretation of common variants and against overcalling in DTC raw data. citeturn8search1turn33search2

## Responsible Results Presentation and Regulatory Compliance

**Scientifically correct approach (risk- and compliance-driven)**
A non-clinical interpretation tool must strictly separate (a) "information presentation" from (b) "medical recommendation". This is not only ethical but regulatory relevant: the FDA publicly warned that many PGx test claims (including tests and software interpretations marketed directly to consumers) have not been reviewed by the FDA and may be insufficiently supported scientifically; therapy changes based on such claims can cause patient harm. citeturn30view0turn9search0

**FDA regulatory framework / labeling logic (relevant for DTC-adjacent reports)**
21 CFR Part 809 describes labeling requirements for in vitro diagnostics; 21 CFR 809.10 requires, among other things, information on intended use, limitations, and performance characteristics. citeturn9search1turn9search9 An FDA training in 2024 on 809.10(b) emphasizes that labeling should provide consistent core information (intended use, limitations/warnings, performance) and that these elements can also be reflected in test report templates. citeturn30view1turn9search5

**FDA and pharmacogenetics: safe reference points**
The FDA maintains a "Table of Pharmacogenetic Associations" (informational, not automatically a DTC free pass) and reiterates therein that genotyping does not replace clinical vigilance and patient management. citeturn19search2 This is useful as a template for safety language in a tool: PGx is context-dependent and not absolute.

**Disclaimers: how established services draw the line (evidence)**
- Nebula explicitly states "informational/educational only", "not intended for diagnostic purpose", "no medical advice". citeturn10search1turn10search5
- SelfDecode describes that the product is not intended for diagnosis/treatment and that no medical decisions should be derived from it. citeturn10search6
Such disclaimers are necessary but not sufficient: without correct allele and quality logic, the output remains dangerous.

**Why "clinical confirmation" must be part of the UI**
A study on clinical confirmation of DTC raw data variants found that a large proportion of variants reported in DTC raw data were false positives upon clinical confirmation and that some variants marked as "increased risk" were clinically classified as benign; the authors emphasize the necessity of clinical confirmation testing. citeturn33search2turn9search7 This is a direct argument for UI flags: "Raw data are not clinically validated; confirmation in a qualified laboratory is required."

**Communication to non-experts: evidence-based report formats**
Patient-friendly reports benefit from a prominent results summary in understandable language and clear interpretation sections; the academic literature explicitly recommends structured, clear presentation to reduce misinterpretation. citeturn19search0turn19search17turn19search3 ClinGen also provides communication/consent frameworks (e.g., CADRe) that structure disclosure strategies. citeturn19search12

**Flagging high-impact findings separately (and why)**
For highly consequential, medically "actionable" genes/conditions, clinical genomics has the ACMG Secondary Findings (SF) policy with curated lists (e.g., SF v3.2, v3.3) as a framework for responsible return in clinical sequencing. citeturn33search14turn33search1 For a consumer tool, this implies a safety policy: **if** (and only if) a variant appears highly relevant after (i) allele match, (ii) high ClinVar review quality, (iii) frequency plausibility, and (iv) plausible inheritance/genotype constellation, it must be output as "High impact, confirm clinically", not as a diagnosis. The DTC false-positive data support this rigor. citeturn33search2turn8search1turn24view0

**Data sources/implementation artifacts**
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
citeturn30view0turn9search0turn9search1turn9search9turn19search2turn33search2turn19search0

**Common pitfalls**
(1) Treating disclaimers as a substitute for technical correctness. citeturn30view0turn33search2
(2) Presenting PGx as "automatic dosing instructions" without context; the FDA warns against non-validated claims. citeturn30view0turn19search2
(3) Displaying high-impact variants like normal "info" hits instead of as "confirm clinically". citeturn33search2turn33search14

**Priority (patient safety)**
Critical. This area determines whether a technically correct tool is clinically misused.

## Reference Implementations and Best Practices

**Scientifically correct approach: derive from reference tools, do not reinvent**
For PGx, PharmCAT is the reference implementation: PharmCAT accepts VCF/"outside call", identifies PGx genotypes, infers star alleles, and generates reports with guideline recommendations (CPIC/DPWG); it consists of a preprocessor and a Named Allele Matcher. citeturn1search1turn1search13turn1search5 From this, the direct architectural derivation follows: (a) VCF normalization, (b) named allele matching data-driven from tables, (c) report generation strictly separated from matching.

For generic variant annotation (ClinVar, gnomAD, custom tables), established annotators are relevant:

- OpenCRAVAT is modular, can be run locally, and is pipeline-capable. citeturn12search5turn12search1
- Ensembl VEP documents very concretely how VCF entries are internally normalized/trimmed and that allele comparison can optionally be disabled ("don't compare alleles"), which serves as a warning signal: "colocated rsID" is not automatically the same allele state. citeturn12search7turn12search11
- GA4GH VRS defines normalization as a canonical representation for cross-system comparability. citeturn12search0turn12search8

**GA4GH/VCF normalization: what is mandatory and when**
Variation Normalization (VRS) aims at canonical forms to make "equivalent" variants unambiguous. citeturn12search0turn12search4 For VCFs, the classic steps are: (i) multi-allelic split, (ii) REF check against FASTA, (iii) left-alignment of indels in repetitive regions, (iv) trimming of shared bases; bcftools norm is the de facto standard utility for such operations. citeturn12search2turn12search14turn12search7

For pure SNP array genotypes, left-alignment/trimming is rarely relevant (because they are mostly SNVs), but as soon as you generate VCF from array data or ingest actual VCFs, normalization must be part of the standard pipeline. citeturn1search1turn12search2turn12search0

**ClinVar data: VCF vs Full Release**
The ClinVar FTP primer makes clear that ClinVar VCF covers only certain variant types (simple alleles, <10kb, precise); "complete" ClinVar coverage may require the Full Release (XML/TXT). citeturn13search9turn13search17 For a local consumer tool, VCF is often sufficient, but the limitations must be visible in the UI ("coverage limitations").

**Best-practice pipeline (compact, implementable)**
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
  run_named_allele_matching_data_driven (PharmCAT-like)
  output_diplotype + phenotype + guideline pointers
  carry_forward "missing sites" + "CNV unknown" flags
```

**Data sources/implementation artifacts**
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
citeturn1search1turn1search5turn1search13turn12search5turn12search7turn12search0turn12search2

**Common pitfalls**
(1) Treating annotator output (rsID co-location) as identical alleles; VEP explicitly documents that allele comparison is a separate decision, not implicitly safe. citeturn12search11turn12search7
(2) Implementing PGx as "rule-based hardcoding" instead of data-driven tables (CPIC xlsx) → update drift. citeturn27search0turn28search24
(3) Treating ClinVar VCF as complete → coverage gaps. citeturn13search9turn13search17
(4) Skipping normalization/REF check → silent mismatches for indels/repeats. citeturn12search2turn12search0

**Priority (patient safety)**
High. Reference architectures (PharmCAT, VEP/VRS normalization) provide the safest blueprints for avoiding the main classes of silent failure (allele/strand/representation mismatch). citeturn1search13turn12search0turn12search11
