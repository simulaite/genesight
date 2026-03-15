//! Parser for VCF (Variant Call Format) files.
//!
//! Standard bioinformatics format. We extract rsIDs and genotype calls
//! from the GT field in the sample column. Multi-allelic records are split
//! into separate biallelic variants, and alleles are trimmed to their
//! minimal representation via the [`normalizer`](crate::normalizer) module.

use crate::models::{Genotype, SourceFormat, Variant};
use crate::normalizer::{normalize_vcf_record, NormalizationStatus};

use super::ParseError;

/// Parse VCF content into variants.
///
/// Multi-allelic records (e.g., `ALT=G,T`) are split into separate `Variant`
/// entries, each with a correctly remapped genotype. Allele trimming is applied
/// to produce minimal representations. Symbolic alleles (e.g., `<DEL>`) produce
/// a variant with `Genotype::NoCall`.
pub fn parse(content: &str) -> Result<Vec<Variant>, ParseError> {
    let mut variants = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 10 {
            return Err(ParseError::InvalidLine {
                line: line_num + 1,
                reason: format!("expected at least 10 columns, found {}", cols.len()),
            });
        }

        let chromosome = cols[0];
        let position: u64 = cols[1].parse().map_err(|_| ParseError::InvalidLine {
            line: line_num + 1,
            reason: format!("invalid position: {}", cols[1]),
        })?;
        let id = cols[2];
        let ref_allele = cols[3];
        let alt_field = cols[4];

        // rsID: "." means no ID assigned
        let rsid = if id == "." { None } else { Some(id) };

        // Split ALT field on comma for multi-allelic sites
        let alt_alleles: Vec<&str> = alt_field.split(',').collect();

        // Parse FORMAT and SAMPLE columns to extract GT
        let format_fields: Vec<&str> = cols[8].split(':').collect();
        let sample_fields: Vec<&str> = cols[9].split(':').collect();

        let gt_str = format_fields
            .iter()
            .position(|&f| f == "GT")
            .and_then(|idx| sample_fields.get(idx).copied())
            .unwrap_or("./.");

        // Normalize: split multi-allelic and trim alleles
        let normalized = normalize_vcf_record(position, ref_allele, &alt_alleles, gt_str);

        for record in normalized {
            let genotype = resolve_genotype(
                &record.gt,
                &record.ref_allele,
                &record.alt_allele,
                &record.status,
            );

            variants.push(Variant {
                rsid: rsid.map(|s| s.to_string()),
                chromosome: chromosome.to_string(),
                position: record.position,
                genotype,
                source_format: SourceFormat::Vcf,
            });
        }
    }

    Ok(variants)
}

/// Resolve a normalized GT string into a `Genotype` value.
///
/// This replaces the old `parse_gt` function. It handles:
/// - Diploid GT: `"0/1"`, `"1/1"`, `"0|1"`, etc.
/// - Hemizygous GT: `"1"` (single allele, no separator)
/// - Missing calls: `"."`, `"./."`, etc.
/// - Multi-base alleles (indels): uses `Genotype::Indel`
/// - Normalization failures: always `Genotype::NoCall`
fn resolve_genotype(
    gt: &str,
    ref_allele: &str,
    alt_allele: &str,
    status: &NormalizationStatus,
) -> Genotype {
    // Symbolic or otherwise failed normalization -> NoCall
    if matches!(status, NormalizationStatus::NormalizationFailed(_)) {
        return Genotype::NoCall;
    }

    let resolve_allele = |idx: &str| -> Option<&str> {
        match idx {
            "0" => Some(ref_allele),
            "1" => Some(alt_allele),
            "." => None,
            _ => None,
        }
    };

    // Split GT into allele indices
    let allele_indices: Vec<&str> = if gt.contains('/') {
        gt.split('/').collect()
    } else if gt.contains('|') {
        gt.split('|').collect()
    } else {
        // Hemizygous: single allele, no separator
        match resolve_allele(gt) {
            Some(allele) => return allele_to_genotype_hom(allele),
            None => return Genotype::NoCall,
        }
    };

    match allele_indices.len() {
        1 => match resolve_allele(allele_indices[0]) {
            Some(allele) => allele_to_genotype_hom(allele),
            None => Genotype::NoCall,
        },
        2 => {
            let a = resolve_allele(allele_indices[0]);
            let b = resolve_allele(allele_indices[1]);
            match (a, b) {
                (Some(a_str), Some(b_str)) => alleles_to_genotype(a_str, b_str),
                (Some(a_str), None) => allele_to_genotype_hom(a_str),
                (None, Some(b_str)) => allele_to_genotype_hom(b_str),
                (None, None) => Genotype::NoCall,
            }
        }
        _ => Genotype::NoCall,
    }
}

/// Convert a single allele string to a homozygous genotype.
fn allele_to_genotype_hom(allele: &str) -> Genotype {
    if allele.len() == 1 {
        match allele.chars().next() {
            Some(c) => Genotype::Homozygous(c),
            None => Genotype::NoCall,
        }
    } else {
        // Multi-base allele (indel)
        Genotype::Indel(format!("{allele}{allele}"))
    }
}

/// Convert two allele strings to a genotype.
fn alleles_to_genotype(a: &str, b: &str) -> Genotype {
    match (a.len(), b.len()) {
        (1, 1) => {
            let ca = a.chars().next().unwrap_or('?');
            let cb = b.chars().next().unwrap_or('?');
            if ca == cb {
                Genotype::Homozygous(ca)
            } else {
                Genotype::Heterozygous(ca, cb)
            }
        }
        _ => {
            // At least one allele is multi-base (indel)
            Genotype::Indel(format!("{a}/{b}"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Test 9: Existing VCF test fixture still parses correctly
    // -----------------------------------------------------------------------
    #[test]
    fn parse_basic_vcf() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t82154\trs4477212\tG\tA\t.\tPASS\t.\tGT\t0/1
1\t752566\trs3094315\tA\tG\t.\tPASS\t.\tGT\t1/1
1\t800000\t.\tC\tT\t.\tPASS\t.\tGT\t0/0
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 3);
        assert_eq!(variants[0].rsid, Some("rs4477212".to_string()));
        assert_eq!(variants[0].genotype, Genotype::Heterozygous('G', 'A'));
        assert_eq!(variants[1].genotype, Genotype::Homozygous('G'));
        // Novel variant (no rsID) still parsed
        assert_eq!(variants[2].rsid, None);
        assert_eq!(variants[2].genotype, Genotype::Homozygous('C'));
    }

    // -----------------------------------------------------------------------
    // Multi-allelic splitting through the full parser
    // -----------------------------------------------------------------------
    #[test]
    fn parse_multiallelic_splits_into_two_variants() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs123\tA\tG,T\t50\tPASS\t.\tGT\t0/1
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 2);

        // First split: REF=A, ALT=G, het A/G
        assert_eq!(variants[0].rsid, Some("rs123".to_string()));
        assert_eq!(variants[0].genotype, Genotype::Heterozygous('A', 'G'));
        assert_eq!(variants[0].position, 100);

        // Second split: REF=A, ALT=T, GT remapped to 0/. -> hom REF
        assert_eq!(variants[1].rsid, Some("rs123".to_string()));
        assert_eq!(variants[1].genotype, Genotype::Homozygous('A'));
        assert_eq!(variants[1].position, 100);
    }

    #[test]
    fn parse_multiallelic_1_2_both_het() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs456\tA\tG,T\t50\tPASS\t.\tGT\t1/2
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 2);

        // Split 0: GT 1/. -> hom G (one allele is G, other is missing -> treat as hom)
        assert_eq!(variants[0].genotype, Genotype::Homozygous('G'));

        // Split 1: GT ./1 -> hom T
        assert_eq!(variants[1].genotype, Genotype::Homozygous('T'));
    }

    // -----------------------------------------------------------------------
    // Symbolic allele -> NoCall
    // -----------------------------------------------------------------------
    #[test]
    fn parse_symbolic_allele_nocall() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs789\tA\t<DEL>\t50\tPASS\t.\tGT\t0/1
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype, Genotype::NoCall);
    }

    // -----------------------------------------------------------------------
    // Hemizygous GT "1" -> Homozygous(alt)
    // -----------------------------------------------------------------------
    #[test]
    fn parse_hemizygous_gt() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
X\t100\trs999\tA\tG\t50\tPASS\t.\tGT\t1
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype, Genotype::Homozygous('G'));
        assert_eq!(variants[0].chromosome, "X");
    }

    // -----------------------------------------------------------------------
    // Trimming through the full parser
    // -----------------------------------------------------------------------
    #[test]
    fn parse_trimmed_alleles() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs100\tATG\tACG\t50\tPASS\t.\tGT\t0/1
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 1);
        // After trimming: REF=T, ALT=C at pos 101
        assert_eq!(variants[0].position, 101);
        assert_eq!(variants[0].genotype, Genotype::Heterozygous('T', 'C'));
    }

    // -----------------------------------------------------------------------
    // Phased genotype
    // -----------------------------------------------------------------------
    #[test]
    fn parse_phased_gt() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs200\tA\tG\t50\tPASS\t.\tGT\t1|0
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype, Genotype::Heterozygous('G', 'A'));
    }

    // -----------------------------------------------------------------------
    // No GT field -> NoCall
    // -----------------------------------------------------------------------
    #[test]
    fn parse_missing_gt_field() {
        let content = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\trs300\tA\tG\t50\tPASS\t.\tDP\t30
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 1);
        // GT not found in FORMAT -> ./. default -> NoCall
        assert_eq!(variants[0].genotype, Genotype::NoCall);
    }
}
