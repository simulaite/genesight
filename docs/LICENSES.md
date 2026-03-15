# GeneSight — License Overview

## Project License

**GeneSight** is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

### Why GPL-3.0?

1. **Copyleft:** Ensures that forks and derivative works also remain open source
2. **Compatibility:** GPL-3.0 is compatible with CC-BY-NC-SA 3.0 (SNPedia) in a non-commercial context
3. **Community:** Encourages contributions back to the project
4. **Privacy signal:** Underscores that the code is inspectable — important when dealing with DNA data

---

## Database Licenses in Detail

### Public Domain (no restrictions)

| Database | Rationale |
|----------|-----------|
| **ClinVar** | US Government Work — not subject to copyright per 17 U.S.C. § 105 |
| **dbSNP** | US Government Work — same rationale |

These data may be used, distributed, and commercially exploited without any restriction.
Attribution is not legally required but is considered good scientific practice.

### CC-BY-NC-SA 3.0 US (SNPedia)

**Permitted:**
- Share — copy and redistribute in any format
- Adapt — remix, transform, and build upon the material

**Conditions:**
- **BY (Attribution):** SNPedia must be credited as the source
- **NC (Non-commercial):** No commercial use without a separate license
- **SA (Share-Alike):** Derivative works must be released under the same license

**What this means for GeneSight:**
- ✅ The open-source project may freely use SNPedia data
- ✅ Users may use the tool for personal DNA analysis
- ✅ Academic research is permitted
- ❌ A commercial fork would need to remove SNPedia data
- → **Architectural decision:** SNPedia data is treated as a separate, optional download rather than bundled in the repository

### CC-BY-SA 4.0 (PharmGKB)

**Permitted:**
- Share and adapt, including for commercial purposes

**Conditions:**
- **BY (Attribution):** PharmGKB must be credited as the source
- **SA (Share-Alike):** Derivative works must be released under the same or a compatible license

**Note:** PharmGKB has additional terms of use for commercial usage.
For academic and non-commercial open-source use: freely available.

### Open Access (GWAS Catalog, gnomAD)

| Database | License | Details |
|----------|---------|---------|
| **GWAS Catalog** | EMBL-EBI Terms of Use | Free for all purposes; attribution requested |
| **gnomAD** | ODC Open Database License | Free for all purposes including commercial use |

---

## License Compatibility Matrix

```
GPL-3.0 (GeneSight Code)
├── ✅ Public Domain (ClinVar, dbSNP) — no conflict
├── ✅ CC-BY-SA 4.0 (PharmGKB) — compatible with GPL-3.0
├── ✅ ODC-ODbL (gnomAD) — compatible
├── ✅ Open Access (GWAS Catalog) — compatible
└── ⚠️ CC-BY-NC-SA 3.0 (SNPedia) — compatible ONLY if:
    - The overall project is used non-commercially, OR
    - SNPedia data is treated as a separate, optional download
```

### Architectural Solution for CC-BY-NC-SA

```
genesight.db (Main database)
├── clinvar     → Public Domain ✅
├── dbsnp       → Public Domain ✅
├── gwas        → Open Access ✅
├── gnomad      → ODC-ODbL ✅
└── pharmgkb    → CC-BY-SA 4.0 ✅

snpedia.db (Separate, optional database)
└── snpedia     → CC-BY-NC-SA 3.0 ⚠️
```

The CLI tool works without `snpedia.db` — the user can optionally
download it via `genesight fetch --include-snpedia`.

---

## Attribution in Code

Every generated report MUST include the following attribution block:

```markdown
---
## Data Sources

This report was generated with GeneSight (GPL-3.0).
The following data sources were used:

- **ClinVar** — NCBI, National Library of Medicine (Public Domain)
- **SNPedia** — snpedia.com (CC-BY-NC-SA 3.0) [if used]
- **GWAS Catalog** — NHGRI-EBI (Open Access)
- **gnomAD** — Broad Institute (ODC-ODbL)
- **PharmGKB** — pharmgkb.org (CC-BY-SA 4.0) [if used]

This report is NOT diagnostic. Consult a physician or
genetic counselor for medical decisions.
---
```

---

## Third-Party Dependencies (Rust Crates)

All Rust crates must be GPL-3.0-compatible. Permitted licenses:

- MIT ✅
- Apache-2.0 ✅
- BSD-2-Clause / BSD-3-Clause ✅
- ISC ✅
- MPL-2.0 ✅ (file-level copyleft)
- GPL-3.0 ✅
- Unlicense ✅

**Not permitted:**
- GPL-2.0-only (without "or later") — incompatible with GPL-3.0
- AGPL-3.0 — would require server operators to disclose source
- Proprietary licenses
