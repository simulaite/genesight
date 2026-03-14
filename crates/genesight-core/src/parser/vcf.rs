//! Parser for VCF (Variant Call Format) files.
//!
//! Standard bioinformatics format. We extract rsIDs and genotype calls
//! from the GT field in the sample column.

use crate::models::{Genotype, SourceFormat, Variant};

use super::ParseError;

/// Parse VCF content into variants.
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

        let chromosome = cols[0].to_string();
        let position: u64 = cols[1].parse().map_err(|_| ParseError::InvalidLine {
            line: line_num + 1,
            reason: format!("invalid position: {}", cols[1]),
        })?;
        let id = cols[2];
        let ref_allele = cols[3];
        let alt_allele = cols[4];

        // Skip entries without a proper rsID
        let rsid = if id == "." {
            None
        } else {
            Some(id.to_string())
        };

        // Parse FORMAT and SAMPLE columns to extract GT
        let format_fields: Vec<&str> = cols[8].split(':').collect();
        let sample_fields: Vec<&str> = cols[9].split(':').collect();

        let gt_index = format_fields.iter().position(|&f| f == "GT");
        let genotype = match gt_index {
            Some(idx) if idx < sample_fields.len() => {
                parse_gt(sample_fields[idx], ref_allele, alt_allele)
            }
            _ => Genotype::NoCall,
        };

        variants.push(Variant {
            rsid,
            chromosome,
            position,
            genotype,
            source_format: SourceFormat::Vcf,
        });
    }

    Ok(variants)
}

fn parse_gt(gt: &str, ref_allele: &str, alt_allele: &str) -> Genotype {
    let alleles: Vec<&str> = if gt.contains('/') {
        gt.split('/').collect()
    } else if gt.contains('|') {
        gt.split('|').collect()
    } else {
        return Genotype::NoCall;
    };

    let resolve = |idx: &str| -> Option<char> {
        match idx {
            "0" => ref_allele.chars().next(),
            "1" => alt_allele.chars().next(),
            "." => None,
            _ => None,
        }
    };

    match alleles.len() {
        1 => match resolve(alleles[0]) {
            Some(a) => Genotype::Homozygous(a),
            None => Genotype::NoCall,
        },
        2 => match (resolve(alleles[0]), resolve(alleles[1])) {
            (Some(a), Some(b)) if a == b => Genotype::Homozygous(a),
            (Some(a), Some(b)) => Genotype::Heterozygous(a, b),
            (Some(a), None) => Genotype::Homozygous(a),
            _ => Genotype::NoCall,
        },
        _ => Genotype::NoCall,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
