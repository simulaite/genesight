#!/usr/bin/env python3
"""Build curated seed SQLite databases for testing the GeneSight CLI tool.

Creates two databases:
  - genesight.db  — Main database (variants, clinvar, gwas, frequencies, pharmacogenomics)
  - snpedia.db    — Optional SNPedia database (snpedia)

All data is curated from public genome databases using realistic values.
No real personal DNA data is included.

Usage:
    python3 build_seed_db.py [--output-dir DIR]
"""

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

# ============================================================
# Variant definitions — GRCh38 coordinates
# ============================================================

VARIANTS = [
    # rsid, chromosome, position, ref_allele, alt_allele
    # --- ClinVar monogenic / carrier ---
    ("rs28897696", "17", 43093449, "A", "G"),       # BRCA1
    ("rs80357906", "17", 43045682, "C", "T"),       # BRCA1 (5382insC region)
    ("rs80356862", "13", 32337326, "G", "A"),       # BRCA2
    ("rs75930367", "7", 117559590, "G", "A"),       # CFTR (G542X)
    ("rs121908769", "7", 117587811, "G", "A"),      # CFTR (G551D)
    ("rs334", "11", 5227002, "T", "A"),             # HBB (sickle cell)
    ("rs63750447", "21", 25891796, "C", "T"),       # APP (Alzheimer)
    ("rs121909001", "15", 72346580, "C", "T"),      # HEXA (Tay-Sachs)
    ("rs1800553", "13", 20189473, "G", "A"),        # GJB2 (hearing loss)
    ("rs11571833", "13", 32398489, "A", "T"),       # BRCA2 K3326X
    # --- Pharmacogenomics ---
    ("rs3892097", "22", 42128945, "G", "A"),        # CYP2D6*4
    ("rs4244285", "10", 94781859, "G", "A"),        # CYP2C19*2
    ("rs4986893", "10", 94780653, "G", "A"),        # CYP2C19*3
    ("rs1799853", "10", 94942290, "C", "T"),        # CYP2C9*2
    ("rs1057910", "10", 94981296, "A", "C"),        # CYP2C9*3
    ("rs9923231", "16", 31096368, "C", "T"),        # VKORC1
    ("rs1045642", "7", 87509329, "A", "G"),         # ABCB1
    ("rs12248560", "10", 94761900, "C", "T"),       # CYP2C19*17
    ("rs4149056", "12", 21178615, "T", "C"),         # SLCO1B1*5
    # --- GWAS polygenic ---
    ("rs9939609", "16", 53786615, "T", "A"),        # FTO
    ("rs7903146", "10", 112998590, "C", "T"),       # TCF7L2
    ("rs10811661", "9", 22134094, "T", "C"),        # CDKN2A/B (T2D)
    ("rs4402960", "3", 185511687, "G", "T"),        # IGF2BP2
    ("rs1333049", "9", 22125503, "C", "G"),         # CDKN2A/B (CAD)
    ("rs10757274", "9", 22098619, "A", "G"),        # 9p21
    ("rs2187668", "6", 32605884, "T", "C"),         # HLA-DQ2.5
    ("rs7412", "19", 44908822, "C", "T"),           # APOE e2
    ("rs429358", "19", 44908684, "T", "C"),         # APOE e4
    ("rs1800497", "11", 113400106, "G", "A"),       # DRD2/ANKK1
    ("rs6265", "11", 27658369, "C", "T"),           # BDNF
    ("rs53576", "3", 8762685, "G", "A"),            # OXTR
    # --- SNPedia-only ---
    ("rs1815739", "11", 66560624, "C", "T"),        # ACTN3
    ("rs4988235", "2", 136608646, "G", "A"),        # MCM6/LCT
    ("rs1805007", "16", 89919709, "C", "T"),        # MC1R
    ("rs12913832", "15", 28120472, "A", "G"),       # HERC2/OCA2
    ("rs1426654", "15", 48426484, "A", "G"),        # SLC24A5
    ("rs4680", "22", 19963748, "G", "A"),           # COMT
    # --- New PGx expansion (CYP3A5, DPYD, TPMT, NUDT15) ---
    ("rs776746", "7", 99652770, "A", "G"),            # CYP3A5*3
    ("rs3918290", "1", 97915614, "C", "T"),           # DPYD*2A
    ("rs55886062", "1", 97981343, "A", "C"),          # DPYD*13
    ("rs67376798", "1", 97450058, "T", "A"),          # DPYD c.2846A>T
    ("rs75017182", "1", 98348885, "C", "T"),          # DPYD HapB3
    ("rs1800462", "6", 18130918, "G", "C"),           # TPMT*2
    ("rs1800460", "6", 18139422, "C", "T"),           # TPMT*3B/*3A
    ("rs1142345", "6", 18143606, "T", "C"),           # TPMT*3C/*3A
    ("rs116855232", "13", 45003967, "C", "T"),        # NUDT15*3
]

# ============================================================
# ClinVar data
# ============================================================

CLINVAR_ENTRIES = [
    # rsid, clinical_significance, review_status (stars), conditions (JSON), gene_symbol, last_updated, classification_type
    (
        "rs28897696", "Pathogenic", 4,
        json.dumps(["Hereditary breast and ovarian cancer syndrome", "Breast neoplasm"]),
        "BRCA1", "2024-06-15", "germline",
    ),
    (
        "rs80357906", "Pathogenic", 4,
        json.dumps(["Hereditary breast and ovarian cancer syndrome", "Ovarian neoplasm"]),
        "BRCA1", "2024-06-15", "germline",
    ),
    (
        "rs80356862", "Pathogenic", 3,
        json.dumps(["Hereditary breast cancer", "Familial cancer of breast"]),
        "BRCA2", "2024-05-20", "germline",
    ),
    (
        "rs75930367", "Pathogenic", 4,
        json.dumps(["Cystic fibrosis"]),
        "CFTR", "2024-07-01", "germline",
    ),
    (
        "rs121908769", "Pathogenic", 4,
        json.dumps(["Cystic fibrosis"]),
        "CFTR", "2024-07-01", "germline",
    ),
    (
        "rs334", "Pathogenic", 4,
        json.dumps(["Sickle cell disease", "Sickle cell trait"]),
        "HBB", "2024-04-10", "germline",
    ),
    (
        "rs63750447", "Pathogenic", 3,
        json.dumps(["Alzheimer disease, familial, type 1"]),
        "APP", "2024-03-20", "germline",
    ),
    (
        "rs121909001", "Pathogenic", 3,
        json.dumps(["Tay-Sachs disease"]),
        "HEXA", "2024-02-15", "germline",
    ),
    (
        "rs1800553", "Pathogenic", 3,
        json.dumps(["Nonsyndromic hearing loss, autosomal recessive"]),
        "GJB2", "2024-01-30", "germline",
    ),
    (
        "rs11571833", "Pathogenic", 2,
        json.dumps(["Hereditary breast cancer", "Fanconi anemia"]),
        "BRCA2", "2024-05-20", "germline",
    ),
]

# ============================================================
# Pharmacogenomics data
# ============================================================

PHARMACOGENOMICS_ENTRIES = [
    # rsid, drug, phenotype_category, evidence_level, clinical_recommendation, gene_symbol
    (
        "rs3892097", "codeine", "Poor Metabolizer", "1A",
        "Avoid codeine. Use alternative analgesics such as morphine or a non-opioid. "
        "CYP2D6 poor metabolizers cannot convert codeine to its active metabolite morphine.",
        "CYP2D6",
    ),
    (
        "rs3892097", "tramadol", "Poor Metabolizer", "1A",
        "Avoid tramadol. Use alternative analgesics. CYP2D6 poor metabolizers have "
        "significantly reduced analgesic effect.",
        "CYP2D6",
    ),
    (
        "rs4244285", "clopidogrel", "Poor Metabolizer", "1A",
        "Use alternative antiplatelet therapy (e.g., prasugrel, ticagrelor). "
        "CYP2C19 poor metabolizers have reduced clopidogrel activation.",
        "CYP2C19",
    ),
    (
        "rs4986893", "clopidogrel", "Poor Metabolizer", "1A",
        "Use alternative antiplatelet therapy (e.g., prasugrel, ticagrelor). "
        "CYP2C19*3 is a no-function allele.",
        "CYP2C19",
    ),
    (
        "rs1799853", "warfarin", "Intermediate Metabolizer", "1A",
        "Reduce initial warfarin dose by 20-40%. CYP2C9*2 carriers metabolize "
        "warfarin more slowly, increasing bleeding risk.",
        "CYP2C9",
    ),
    (
        "rs1057910", "warfarin", "Poor Metabolizer", "1A",
        "Reduce initial warfarin dose by 40-70%. CYP2C9*3 causes substantially "
        "reduced warfarin metabolism and significantly increased bleeding risk.",
        "CYP2C9",
    ),
    (
        "rs9923231", "warfarin", "Increased Sensitivity", "1A",
        "VKORC1 -1639G>A. Carriers of the A allele require lower warfarin doses. "
        "Use pharmacogenomic dosing algorithms (e.g., IWPC or Gage).",
        "VKORC1",
    ),
    (
        "rs1045642", "multiple drugs", "Altered Drug Transport", "2A",
        "ABCB1 3435C>T polymorphism affects P-glycoprotein expression and may alter "
        "bioavailability of digoxin, cyclosporine, and other P-gp substrates.",
        "ABCB1",
    ),
    (
        "rs12248560", "clopidogrel", "Ultrarapid Metabolizer", "1A",
        "CYP2C19*17 increases enzyme activity. Standard clopidogrel dosing is effective. "
        "May have increased bleeding risk with standard doses of some medications.",
        "CYP2C19",
    ),
]

# ============================================================
# GWAS data
# ============================================================

GWAS_ENTRIES = [
    # rsid, trait, p_value, odds_ratio, beta, risk_allele, risk_allele_frequency, pubmed_id, mapped_gene
    (
        "rs9939609", "Body mass index", 1e-120, 1.31, None,
        "A", 0.42, "17554300", "FTO",
    ),
    (
        "rs7903146", "Type 2 diabetes", 1e-75, 1.37, None,
        "T", 0.30, "17463246", "TCF7L2",
    ),
    (
        "rs10811661", "Type 2 diabetes", 1e-30, 1.20, None,
        "T", 0.83, "17463249", "CDKN2A/CDKN2B",
    ),
    (
        "rs4402960", "Type 2 diabetes", 1e-15, 1.14, None,
        "T", 0.30, "17463249", "IGF2BP2",
    ),
    (
        "rs1333049", "Coronary artery disease", 1e-20, 1.29, None,
        "C", 0.47, "17634449", "CDKN2A/CDKN2B",
    ),
    (
        "rs10757274", "Coronary artery disease", 1e-25, 1.28, None,
        "G", 0.49, "17478679", "CDKN2A/CDKN2B",
    ),
    (
        "rs2187668", "Celiac disease", 1e-200, 7.0, None,
        "T", 0.12, "20190752", "HLA-DQA1",
    ),
    (
        "rs7412", "LDL cholesterol levels", 1e-100, None, -0.45,
        "T", 0.08, "20686565", "APOE",
    ),
    (
        "rs429358", "Alzheimer disease", 1e-50, 3.68, None,
        "C", 0.15, "24162737", "APOE",
    ),
    (
        "rs1800497", "Addiction susceptibility", 5e-8, 1.18, None,
        "A", 0.33, "21242757", "ANKK1/DRD2",
    ),
    (
        "rs6265", "Cognitive function", 1e-10, None, -0.05,
        "T", 0.20, "22472876", "BDNF",
    ),
    (
        "rs53576", "Social behavior", 5e-8, None, None,
        "A", 0.40, "22069391", "OXTR",
    ),
]

# ============================================================
# Allele frequencies (gnomAD-like values)
# ============================================================

FREQUENCY_ENTRIES = [
    # rsid, af_total, af_afr, af_amr, af_eas, af_eur, af_sas, source
    # --- ClinVar variants (rare pathogenic) ---
    ("rs28897696",  0.0003,  0.0001, 0.0002, 0.0001, 0.0005, 0.0001, "gnomad"),
    ("rs80357906",  0.0004,  0.0001, 0.0002, 0.0000, 0.0010, 0.0001, "gnomad"),
    ("rs80356862",  0.0002,  0.0001, 0.0001, 0.0000, 0.0004, 0.0001, "gnomad"),
    ("rs75930367",  0.0008,  0.0001, 0.0004, 0.0001, 0.0020, 0.0003, "gnomad"),
    ("rs121908769", 0.0003,  0.0000, 0.0001, 0.0000, 0.0008, 0.0001, "gnomad"),
    ("rs334",       0.0150,  0.0700, 0.0050, 0.0001, 0.0005, 0.0050, "gnomad"),
    ("rs63750447",  0.0001,  0.0000, 0.0001, 0.0000, 0.0002, 0.0001, "gnomad"),
    ("rs121909001", 0.0005,  0.0001, 0.0003, 0.0000, 0.0012, 0.0001, "gnomad"),
    ("rs1800553",   0.0003,  0.0002, 0.0002, 0.0001, 0.0006, 0.0001, "gnomad"),
    ("rs11571833",  0.0060,  0.0020, 0.0040, 0.0010, 0.0110, 0.0020, "gnomad"),
    # --- Pharmacogenomics variants ---
    ("rs3892097",   0.2000,  0.0700, 0.1200, 0.0100, 0.2800, 0.1500, "gnomad"),
    ("rs4244285",   0.1500,  0.1800, 0.1200, 0.3000, 0.1500, 0.3400, "gnomad"),
    ("rs4986893",   0.0100,  0.0010, 0.0030, 0.0700, 0.0004, 0.0100, "gnomad"),
    ("rs1799853",   0.0600,  0.0100, 0.0500, 0.0010, 0.1100, 0.0400, "gnomad"),
    ("rs1057910",   0.0400,  0.0100, 0.0400, 0.0400, 0.0700, 0.0800, "gnomad"),
    ("rs9923231",   0.3200,  0.1000, 0.4000, 0.9200, 0.3900, 0.1500, "gnomad"),
    ("rs1045642",   0.4800,  0.1800, 0.4300, 0.6000, 0.5300, 0.5500, "gnomad"),
    ("rs12248560",  0.2100,  0.2400, 0.1400, 0.0100, 0.2100, 0.1500, "gnomad"),
    ("rs4149056",   0.0800,  0.0200, 0.0600, 0.1000, 0.0800, 0.0500, "gnomad"),
    # --- GWAS variants (common) ---
    ("rs9939609",   0.4200,  0.4900, 0.3000, 0.1400, 0.4100, 0.3100, "gnomad"),
    ("rs7903146",   0.3000,  0.2800, 0.2400, 0.0400, 0.3000, 0.3100, "gnomad"),
    ("rs10811661",  0.8300,  0.9300, 0.8700, 0.5800, 0.8100, 0.8800, "gnomad"),
    ("rs4402960",   0.3000,  0.3600, 0.2100, 0.2700, 0.3000, 0.3200, "gnomad"),
    ("rs1333049",   0.4700,  0.3600, 0.4500, 0.5200, 0.4800, 0.5100, "gnomad"),
    ("rs10757274",  0.4900,  0.3400, 0.4800, 0.5500, 0.5000, 0.4900, "gnomad"),
    ("rs2187668",   0.1200,  0.0300, 0.0600, 0.0100, 0.1500, 0.0500, "gnomad"),
    ("rs7412",      0.0800,  0.1100, 0.0500, 0.0900, 0.0700, 0.0500, "gnomad"),
    ("rs429358",    0.1500,  0.2700, 0.1000, 0.0900, 0.1500, 0.0800, "gnomad"),
    ("rs1800497",   0.3300,  0.3200, 0.3900, 0.3800, 0.3200, 0.3600, "gnomad"),
    ("rs6265",      0.2000,  0.0200, 0.1000, 0.4500, 0.2000, 0.1600, "gnomad"),
    ("rs53576",     0.4000,  0.7000, 0.4200, 0.4000, 0.3600, 0.3800, "gnomad"),
    # --- SNPedia-only variants ---
    ("rs1815739",   0.4200,  0.1300, 0.4300, 0.5300, 0.4500, 0.3900, "gnomad"),
    ("rs4988235",   0.2600,  0.1400, 0.3000, 0.0100, 0.4900, 0.2200, "gnomad"),
    ("rs1805007",   0.0600,  0.0000, 0.0200, 0.0000, 0.1100, 0.0100, "gnomad"),
    ("rs12913832",  0.2500,  0.0200, 0.0900, 0.0100, 0.7200, 0.0500, "gnomad"),
    ("rs1426654",   0.4500,  0.0400, 0.5500, 0.0100, 0.9900, 0.6500, "gnomad"),
    ("rs4680",      0.4800,  0.3300, 0.4700, 0.2800, 0.5000, 0.4300, "gnomad"),
    # --- New PGx expansion variants ---
    ("rs776746",    0.9300,  0.2700, 0.8200, 0.7100, 0.9500, 0.8100, "gnomad"),
    ("rs3918290",   0.0100,  0.0020, 0.0050, 0.0010, 0.0180, 0.0030, "gnomad"),
    ("rs55886062",  0.0030,  0.0010, 0.0020, 0.0010, 0.0060, 0.0010, "gnomad"),
    ("rs67376798",  0.0120,  0.0030, 0.0080, 0.0020, 0.0200, 0.0050, "gnomad"),
    ("rs75017182",  0.0160,  0.0040, 0.0100, 0.0030, 0.0260, 0.0070, "gnomad"),
    ("rs1800462",   0.0030,  0.0020, 0.0030, 0.0010, 0.0050, 0.0020, "gnomad"),
    ("rs1800460",   0.0400,  0.0200, 0.0300, 0.0050, 0.0450, 0.0200, "gnomad"),
    ("rs1142345",   0.0450,  0.0550, 0.0350, 0.0220, 0.0480, 0.0250, "gnomad"),
    ("rs116855232", 0.0400,  0.0080, 0.0800, 0.1200, 0.0040, 0.0200, "gnomad"),
]

# ============================================================
# SNPedia data
# ============================================================

SNPEDIA_ENTRIES = [
    # rsid, magnitude, repute, summary, genotype_descriptions (dict)
    (
        "rs1815739", 3.0, "good",
        "Enhanced sprint/power muscle performance",
        {
            "CC": "Functional ACTN3. Likely enhanced sprint/power muscle fiber performance.",
            "CT": "One copy of functional ACTN3. Intermediate sprint/power performance.",
            "TT": "No functional ACTN3 (alpha-actinin-3 deficiency). Endurance-oriented muscle fibers.",
        },
    ),
    (
        "rs4988235", 3.0, "good",
        "Lactose tolerant (European variant)",
        {
            "GG": "Likely lactose intolerant (ancestral). Cannot digest lactose well in adulthood.",
            "GA": "Lactose tolerant. One copy of the European lactase persistence allele.",
            "AA": "Lactose tolerant. Homozygous for European lactase persistence.",
        },
    ),
    (
        "rs1800497", 2.5, "bad",
        "Reduced dopamine receptor density",
        {
            "GG": "Normal DRD2 dopamine receptor density.",
            "GA": "Somewhat reduced DRD2 receptor density. Mildly increased addiction susceptibility.",
            "AA": "Reduced DRD2 receptor density (Taq1A A1/A1). Increased addiction susceptibility.",
        },
    ),
    (
        "rs429358", 4.0, "bad",
        "Increased Alzheimer risk (APOE e4)",
        {
            "TT": "No APOE-e4 alleles. Baseline Alzheimer risk.",
            "TC": "One APOE-e4 allele. Approximately 3x increased Alzheimer risk.",
            "CC": "Two APOE-e4 alleles (homozygous). Approximately 12x increased Alzheimer risk.",
        },
    ),
    (
        "rs7412", 3.0, "good",
        "Protective against Alzheimer's (APOE e2)",
        {
            "CC": "No APOE-e2 alleles. Baseline lipid metabolism.",
            "CT": "One APOE-e2 allele. Mildly protective against Alzheimer disease and lower LDL.",
            "TT": "Two APOE-e2 alleles. Strong protection against Alzheimer disease. Lower LDL cholesterol.",
        },
    ),
    (
        "rs1805007", 2.5, "good",
        "Red hair, fair skin (MC1R R151C variant)",
        {
            "CC": "Normal MC1R function. Typical pigmentation for ethnic background.",
            "CT": "One MC1R variant. Possibly lighter skin, increased freckles, some red hair tendency.",
            "TT": "Homozygous MC1R R151C. Red/auburn hair, very fair skin, increased sun sensitivity.",
        },
    ),
    (
        "rs12913832", 3.0, "good",
        "Blue eye color (HERC2/OCA2 regulatory variant)",
        {
            "AA": "Likely brown eyes. Ancestral allele associated with darker eye color.",
            "AG": "Likely green or hazel eyes. Heterozygous for the blue-eye variant.",
            "GG": "Likely blue eyes. Homozygous for the HERC2 variant that reduces OCA2 expression.",
        },
    ),
    (
        "rs1426654", 2.0, "good",
        "European skin pigmentation variant (SLC24A5 Thr111Ala)",
        {
            "AA": "Ancestral allele. Darker skin pigmentation typical of African/East Asian ancestry.",
            "AG": "Heterozygous. Intermediate skin pigmentation.",
            "GG": "European-associated lighter skin pigmentation. Nearly fixed in European populations.",
        },
    ),
    (
        "rs4680", 2.0, None,
        "Warrior vs Worrier variant (COMT Val158Met)",
        {
            "GG": "Val/Val (Warrior). Higher COMT activity, lower prefrontal dopamine. "
                  "Better stress resilience, slightly lower working memory.",
            "GA": "Val/Met. Intermediate COMT activity. Balanced dopamine levels.",
            "AA": "Met/Met (Worrier). Lower COMT activity, higher prefrontal dopamine. "
                  "Better working memory and attention, but higher stress anxiety.",
        },
    ),
    (
        "rs53576", 2.0, None,
        "Empathy and social behavior variant (OXTR)",
        {
            "GG": "Higher empathy and social sensitivity. Better at reading emotional cues. "
                  "More responsive to social support.",
            "GA": "Intermediate oxytocin receptor expression. Moderate social sensitivity.",
            "AA": "Lower empathy scores in some studies. Less responsive to social cues. "
                  "May be more resilient to social stress.",
        },
    ),
]

# ============================================================
# PGx Allele Definitions (CPIC-style star allele mapping)
# ============================================================

PGX_ALLELE_DEFINITIONS = [
    # gene, allele_name, rsid, alt_allele, function, activity_score
    # --- CYP2C19 ---
    ("CYP2C19", "*2", "rs4244285", "A", "No Function", 0.0),
    ("CYP2C19", "*3", "rs4986893", "A", "No Function", 0.0),
    ("CYP2C19", "*17", "rs12248560", "T", "Increased Function", 1.5),
    # --- CYP2D6 ---
    ("CYP2D6", "*4", "rs3892097", "A", "No Function", 0.0),
    # --- CYP2C9 ---
    ("CYP2C9", "*2", "rs1799853", "T", "Decreased Function", 0.5),
    ("CYP2C9", "*3", "rs1057910", "C", "No Function", 0.0),
    # --- SLCO1B1 ---
    ("SLCO1B1", "*5", "rs4149056", "C", "Decreased Function", 0.5),
    # --- CYP3A5 ---
    ("CYP3A5", "*3", "rs776746", "G", "No Function", 0.0),
    # --- DPYD ---
    ("DPYD", "*2A", "rs3918290", "T", "No Function", 0.0),
    ("DPYD", "*13", "rs55886062", "C", "No Function", 0.0),
    ("DPYD", "c.2846A>T", "rs67376798", "A", "Decreased Function", 0.5),
    ("DPYD", "HapB3", "rs75017182", "T", "Decreased Function", 0.5),
    # --- TPMT ---
    ("TPMT", "*2", "rs1800462", "C", "No Function", 0.0),
    ("TPMT", "*3A", "rs1800460", "T", "No Function", 0.0),
    ("TPMT", "*3A", "rs1142345", "C", "No Function", 0.0),
    ("TPMT", "*3B", "rs1800460", "T", "No Function", 0.0),
    ("TPMT", "*3C", "rs1142345", "C", "No Function", 0.0),
    # --- NUDT15 ---
    ("NUDT15", "*3", "rs116855232", "T", "No Function", 0.0),
    # --- VKORC1 ---
    ("VKORC1", "-1639A", "rs9923231", "T", "Low Warfarin Dose", 0.5),
]

# ============================================================
# PGx Diplotype-to-Phenotype mapping
# ============================================================

PGX_DIPLOTYPE_PHENOTYPES = [
    # gene, diplotype, phenotype, activity_score
    # --- CYP2C19 ---
    ("CYP2C19", "*1/*1", "Normal Metabolizer", 2.0),
    ("CYP2C19", "*1/*2", "Intermediate Metabolizer", 1.0),
    ("CYP2C19", "*2/*2", "Poor Metabolizer", 0.0),
    ("CYP2C19", "*1/*17", "Rapid Metabolizer", 2.5),
    ("CYP2C19", "*17/*17", "Ultrarapid Metabolizer", 3.0),
    # --- CYP2D6 ---
    ("CYP2D6", "*1/*1", "Normal Metabolizer", 2.0),
    ("CYP2D6", "*1/*4", "Intermediate Metabolizer", 1.0),
    ("CYP2D6", "*4/*4", "Poor Metabolizer", 0.0),
    # --- CYP2C9 ---
    ("CYP2C9", "*1/*1", "Normal Metabolizer", 2.0),
    ("CYP2C9", "*1/*2", "Intermediate Metabolizer", 1.5),
    ("CYP2C9", "*1/*3", "Intermediate Metabolizer", 1.0),
    ("CYP2C9", "*2/*3", "Poor Metabolizer", 0.5),
    ("CYP2C9", "*3/*3", "Poor Metabolizer", 0.0),
    # --- SLCO1B1 ---
    ("SLCO1B1", "*1/*1", "Normal Function", 2.0),
    ("SLCO1B1", "*1/*5", "Intermediate Function", 1.5),
    ("SLCO1B1", "*5/*5", "Poor Function", 1.0),
    # --- CYP3A5 ---
    ("CYP3A5", "*1/*1", "Extensive Metabolizer", 2.0),
    ("CYP3A5", "*1/*3", "Intermediate Metabolizer", 1.0),
    ("CYP3A5", "*3/*3", "Poor Metabolizer", 0.0),
    # --- DPYD ---
    ("DPYD", "*1/*1", "Normal DPD Activity", 2.0),
    ("DPYD", "*1/*2A", "Intermediate DPD Activity", 1.0),
    ("DPYD", "*2A/*2A", "Poor DPD Activity (DPD Deficient)", 0.0),
    ("DPYD", "*1/*13", "Intermediate DPD Activity", 1.0),
    ("DPYD", "*1/c.2846A>T", "Intermediate DPD Activity", 1.5),
    ("DPYD", "*1/HapB3", "Intermediate DPD Activity", 1.5),
    # --- TPMT ---
    ("TPMT", "*1/*1", "Normal Metabolizer", 2.0),
    ("TPMT", "*1/*2", "Intermediate Metabolizer", 1.0),
    ("TPMT", "*1/*3A", "Intermediate Metabolizer", 1.0),
    ("TPMT", "*1/*3B", "Intermediate Metabolizer", 1.0),
    ("TPMT", "*1/*3C", "Intermediate Metabolizer", 1.0),
    ("TPMT", "*3A/*3A", "Poor Metabolizer", 0.0),
    ("TPMT", "*3A/*3C", "Poor Metabolizer", 0.0),
    # --- NUDT15 ---
    ("NUDT15", "*1/*1", "Normal Metabolizer", 2.0),
    ("NUDT15", "*1/*3", "Intermediate Metabolizer", 1.0),
    ("NUDT15", "*3/*3", "Poor Metabolizer", 0.0),
    # --- VKORC1 ---
    ("VKORC1", "GG", "Normal Warfarin Sensitivity", 2.0),
    ("VKORC1", "GA", "Increased Warfarin Sensitivity", 1.5),
    ("VKORC1", "AA", "Highly Increased Warfarin Sensitivity", 1.0),
]

# ============================================================
# PGx Drug Recommendations
# ============================================================

PGX_DRUG_RECOMMENDATIONS = [
    # gene, drug, phenotype, recommendation, evidence_level
    # --- CYP2C19 + clopidogrel ---
    ("CYP2C19", "clopidogrel", "Normal Metabolizer",
     "Use clopidogrel at standard dose.", "1A"),
    ("CYP2C19", "clopidogrel", "Intermediate Metabolizer",
     "Consider alternative antiplatelet (prasugrel, ticagrelor) due to reduced clopidogrel activation.", "1A"),
    ("CYP2C19", "clopidogrel", "Poor Metabolizer",
     "Avoid clopidogrel. Use prasugrel or ticagrelor. CYP2C19 poor metabolizers have greatly reduced activation.", "1A"),
    # --- CYP2D6 + codeine ---
    ("CYP2D6", "codeine", "Normal Metabolizer",
     "Use codeine at standard dose.", "1A"),
    ("CYP2D6", "codeine", "Intermediate Metabolizer",
     "Use codeine with caution at reduced dose, or consider alternative analgesic.", "1A"),
    ("CYP2D6", "codeine", "Poor Metabolizer",
     "Avoid codeine. Use alternative analgesics such as morphine or a non-opioid.", "1A"),
    # --- CYP2C9 + warfarin ---
    ("CYP2C9", "warfarin", "Normal Metabolizer",
     "Use standard warfarin dosing algorithm.", "1A"),
    ("CYP2C9", "warfarin", "Intermediate Metabolizer",
     "Reduce initial warfarin dose by 20-40%. Monitor INR closely.", "1A"),
    ("CYP2C9", "warfarin", "Poor Metabolizer",
     "Reduce initial warfarin dose by 40-70%. Significantly increased bleeding risk.", "1A"),
    # --- CYP3A5 + tacrolimus ---
    ("CYP3A5", "tacrolimus", "Extensive Metabolizer",
     "Increase starting dose by 1.5-2x. CYP3A5 expressers metabolize tacrolimus rapidly.", "1A"),
    ("CYP3A5", "tacrolimus", "Intermediate Metabolizer",
     "Increase starting dose by 1.25-1.5x. Monitor trough levels closely.", "1A"),
    ("CYP3A5", "tacrolimus", "Poor Metabolizer",
     "Use standard recommended starting dose. CYP3A5 non-expressers have normal tacrolimus metabolism.", "1A"),
    # --- DPYD + fluorouracil ---
    ("DPYD", "fluorouracil", "Normal DPD Activity",
     "Use standard fluorouracil dosing.", "1A"),
    ("DPYD", "fluorouracil", "Intermediate DPD Activity",
     "Reduce fluorouracil dose by 50%. Increased risk of severe/fatal toxicity.", "1A"),
    ("DPYD", "fluorouracil", "Poor DPD Activity (DPD Deficient)",
     "Avoid fluorouracil and fluorouracil-based regimens. Use alternative chemotherapy.", "1A"),
    # --- DPYD + capecitabine ---
    ("DPYD", "capecitabine", "Normal DPD Activity",
     "Use standard capecitabine dosing.", "1A"),
    ("DPYD", "capecitabine", "Intermediate DPD Activity",
     "Reduce capecitabine dose by 50%. Increased risk of severe/fatal toxicity.", "1A"),
    ("DPYD", "capecitabine", "Poor DPD Activity (DPD Deficient)",
     "Avoid capecitabine. Use alternative chemotherapy.", "1A"),
    # --- TPMT + azathioprine ---
    ("TPMT", "azathioprine", "Normal Metabolizer",
     "Use standard azathioprine dose.", "1A"),
    ("TPMT", "azathioprine", "Intermediate Metabolizer",
     "Reduce azathioprine dose by 30-70%. Start low and titrate based on tolerance.", "1A"),
    ("TPMT", "azathioprine", "Poor Metabolizer",
     "Avoid azathioprine or use 10% of standard dose with frequent monitoring.", "1A"),
    # --- TPMT + mercaptopurine ---
    ("TPMT", "mercaptopurine", "Normal Metabolizer",
     "Use standard mercaptopurine dose.", "1A"),
    ("TPMT", "mercaptopurine", "Intermediate Metabolizer",
     "Reduce mercaptopurine dose by 30-70%. Start low and titrate based on tolerance.", "1A"),
    ("TPMT", "mercaptopurine", "Poor Metabolizer",
     "Avoid mercaptopurine or use 10% of standard dose with frequent monitoring.", "1A"),
    # --- NUDT15 + azathioprine ---
    ("NUDT15", "azathioprine", "Normal Metabolizer",
     "Use standard azathioprine dose.", "1A"),
    ("NUDT15", "azathioprine", "Intermediate Metabolizer",
     "Reduce azathioprine dose by 25-50%. Monitor for myelosuppression.", "1A"),
    ("NUDT15", "azathioprine", "Poor Metabolizer",
     "Avoid azathioprine or use drastically reduced dose (10%). High risk of myelosuppression.", "1A"),
    # --- NUDT15 + mercaptopurine ---
    ("NUDT15", "mercaptopurine", "Normal Metabolizer",
     "Use standard mercaptopurine dose.", "1A"),
    ("NUDT15", "mercaptopurine", "Intermediate Metabolizer",
     "Reduce mercaptopurine dose by 25-50%. Monitor for myelosuppression.", "1A"),
    ("NUDT15", "mercaptopurine", "Poor Metabolizer",
     "Avoid mercaptopurine or use drastically reduced dose (10%). High risk of myelosuppression.", "1A"),
    # --- VKORC1 + warfarin ---
    ("VKORC1", "warfarin", "Normal Warfarin Sensitivity",
     "Use standard warfarin dosing algorithm.", "1A"),
    ("VKORC1", "warfarin", "Increased Warfarin Sensitivity",
     "Reduce initial warfarin dose. Use pharmacogenomic dosing algorithm (IWPC or Gage).", "1A"),
    ("VKORC1", "warfarin", "Highly Increased Warfarin Sensitivity",
     "Significantly reduce initial warfarin dose. Homozygous VKORC1 -1639A requires lowest dose tier.", "1A"),
]


# ============================================================
# Database creation functions
# ============================================================

def create_main_db(db_path: str) -> dict:
    """Create and populate the main genesight.db database.

    Returns a dict of table names to row counts for reporting.
    """
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    cur = conn.cursor()

    # --- Create schema ---
    cur.execute("""
        CREATE TABLE IF NOT EXISTS schema_version (
            version INTEGER PRIMARY KEY,
            applied_at TEXT NOT NULL DEFAULT (datetime('now'))
        )
    """)
    cur.execute("INSERT INTO schema_version (version) VALUES (1)")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS db_metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        )
    """)
    cur.execute("INSERT INTO db_metadata (key, value) VALUES ('assembly', 'GRCh38')")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            rsid TEXT PRIMARY KEY,
            chromosome TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref_allele TEXT,
            alt_allele TEXT
        )
    """)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_variants_chr_pos ON variants(chromosome, position)"
    )

    cur.execute("""
        CREATE TABLE IF NOT EXISTS clinvar (
            rsid TEXT REFERENCES variants(rsid),
            clinical_significance TEXT,
            review_status INTEGER,
            conditions TEXT,
            gene_symbol TEXT,
            last_updated DATE,
            classification_type TEXT DEFAULT 'germline'
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_clinvar_rsid ON clinvar(rsid)")

    cur.execute("""
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
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON gwas(rsid)")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS frequencies (
            rsid TEXT REFERENCES variants(rsid),
            af_total REAL,
            af_afr REAL,
            af_amr REAL,
            af_eas REAL,
            af_eur REAL,
            af_sas REAL,
            source TEXT
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_freq_rsid ON frequencies(rsid)")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS pharmacogenomics (
            rsid TEXT REFERENCES variants(rsid),
            drug TEXT NOT NULL,
            phenotype_category TEXT,
            evidence_level TEXT,
            clinical_recommendation TEXT,
            gene_symbol TEXT
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_pharma_rsid ON pharmacogenomics(rsid)")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS pgx_allele_definitions (
            gene TEXT NOT NULL,
            allele_name TEXT NOT NULL,
            rsid TEXT REFERENCES variants(rsid),
            alt_allele TEXT NOT NULL,
            function TEXT NOT NULL,
            activity_score REAL NOT NULL
        )
    """)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_pgx_allele_gene ON pgx_allele_definitions(gene)"
    )
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_pgx_allele_rsid ON pgx_allele_definitions(rsid)"
    )

    cur.execute("""
        CREATE TABLE IF NOT EXISTS pgx_diplotype_phenotypes (
            gene TEXT NOT NULL,
            diplotype TEXT NOT NULL,
            phenotype TEXT NOT NULL,
            activity_score REAL NOT NULL
        )
    """)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_pgx_diplo_gene ON pgx_diplotype_phenotypes(gene)"
    )

    cur.execute("""
        CREATE TABLE IF NOT EXISTS pgx_drug_recommendations (
            gene TEXT NOT NULL,
            drug TEXT NOT NULL,
            phenotype TEXT NOT NULL,
            recommendation TEXT NOT NULL,
            evidence_level TEXT NOT NULL
        )
    """)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_pgx_drug_gene ON pgx_drug_recommendations(gene)"
    )

    # --- Insert data ---
    cur.executemany(
        "INSERT INTO variants (rsid, chromosome, position, ref_allele, alt_allele) "
        "VALUES (?, ?, ?, ?, ?)",
        VARIANTS,
    )

    cur.executemany(
        "INSERT INTO clinvar (rsid, clinical_significance, review_status, conditions, "
        "gene_symbol, last_updated, classification_type) VALUES (?, ?, ?, ?, ?, ?, ?)",
        CLINVAR_ENTRIES,
    )

    cur.executemany(
        "INSERT INTO gwas (rsid, trait, p_value, odds_ratio, beta, risk_allele, "
        "risk_allele_frequency, pubmed_id, mapped_gene) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        GWAS_ENTRIES,
    )

    cur.executemany(
        "INSERT INTO frequencies (rsid, af_total, af_afr, af_amr, af_eas, af_eur, "
        "af_sas, source) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        FREQUENCY_ENTRIES,
    )

    cur.executemany(
        "INSERT INTO pharmacogenomics (rsid, drug, phenotype_category, evidence_level, "
        "clinical_recommendation, gene_symbol) VALUES (?, ?, ?, ?, ?, ?)",
        PHARMACOGENOMICS_ENTRIES,
    )

    cur.executemany(
        "INSERT INTO pgx_allele_definitions (gene, allele_name, rsid, alt_allele, "
        "function, activity_score) VALUES (?, ?, ?, ?, ?, ?)",
        PGX_ALLELE_DEFINITIONS,
    )

    cur.executemany(
        "INSERT INTO pgx_diplotype_phenotypes (gene, diplotype, phenotype, "
        "activity_score) VALUES (?, ?, ?, ?)",
        PGX_DIPLOTYPE_PHENOTYPES,
    )

    cur.executemany(
        "INSERT INTO pgx_drug_recommendations (gene, drug, phenotype, "
        "recommendation, evidence_level) VALUES (?, ?, ?, ?, ?)",
        PGX_DRUG_RECOMMENDATIONS,
    )

    conn.commit()

    # Collect statistics
    stats = {}
    for table in ["schema_version", "db_metadata", "variants", "clinvar", "gwas",
                   "frequencies", "pharmacogenomics", "pgx_allele_definitions",
                   "pgx_diplotype_phenotypes", "pgx_drug_recommendations"]:
        count = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]  # noqa: S608
        stats[table] = count

    conn.close()
    return stats


def create_snpedia_db(db_path: str) -> dict:
    """Create and populate the optional snpedia.db database.

    Returns a dict of table names to row counts for reporting.
    """
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS snpedia (
            rsid TEXT PRIMARY KEY,
            magnitude REAL,
            repute TEXT,
            summary TEXT,
            genotype_descriptions TEXT
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_snpedia_rsid ON snpedia(rsid)")

    rows = [
        (entry[0], entry[1], entry[2], entry[3], json.dumps(entry[4]))
        for entry in SNPEDIA_ENTRIES
    ]

    cur.executemany(
        "INSERT INTO snpedia (rsid, magnitude, repute, summary, genotype_descriptions) "
        "VALUES (?, ?, ?, ?, ?)",
        rows,
    )

    conn.commit()

    count = cur.execute("SELECT COUNT(*) FROM snpedia").fetchone()[0]
    stats = {"snpedia": count}

    conn.close()
    return stats


# ============================================================
# Main
# ============================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build curated seed SQLite databases for GeneSight testing.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory where genesight.db and snpedia.db will be created (default: current directory)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    main_db_path = output_dir / "genesight.db"
    snpedia_db_path = output_dir / "snpedia.db"

    print(f"Building seed databases in: {output_dir}")
    print()

    # Build main database
    print("Creating genesight.db ...")
    main_stats = create_main_db(str(main_db_path))
    print(f"  genesight.db created: {main_db_path}")
    for table, count in main_stats.items():
        print(f"    {table}: {count} rows")

    print()

    # Build SNPedia database
    print("Creating snpedia.db ...")
    snpedia_stats = create_snpedia_db(str(snpedia_db_path))
    print(f"  snpedia.db created: {snpedia_db_path}")
    for table, count in snpedia_stats.items():
        print(f"    {table}: {count} rows")

    print()
    total_rows = sum(main_stats.values()) + sum(snpedia_stats.values())
    print(f"Done. Total: {total_rows} rows across {len(main_stats) + len(snpedia_stats)} tables.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
