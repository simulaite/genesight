//! Parser for AncestryDNA raw data format.
//!
//! Format: tab-separated with columns `rsid`, `chromosome`, `position`, `allele1`, `allele2`.
//! Comment lines start with `#`.

use crate::models::{Genotype, SourceFormat, Variant};

use super::ParseError;

/// Parse AncestryDNA raw data content into variants.
pub fn parse(content: &str) -> Result<Vec<Variant>, ParseError> {
    let mut variants = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Skip header line
        if trimmed.starts_with("rsid") {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() != 5 {
            return Err(ParseError::InvalidLine {
                line: line_num + 1,
                reason: format!("expected 5 columns, found {}", cols.len()),
            });
        }

        let rsid = cols[0].to_string();
        let chromosome = cols[1].to_string();
        let position: u64 = cols[2].parse().map_err(|_| ParseError::InvalidLine {
            line: line_num + 1,
            reason: format!("invalid position: {}", cols[2]),
        })?;

        let genotype = parse_alleles(cols[3], cols[4]);

        variants.push(Variant {
            rsid: Some(rsid),
            chromosome,
            position,
            genotype,
            source_format: SourceFormat::AncestryDNA,
        });
    }

    Ok(variants)
}

fn parse_alleles(a1: &str, a2: &str) -> Genotype {
    match (a1, a2) {
        ("0", "0") | ("-", "-") => Genotype::NoCall,
        (a, b) if a.len() == 1 && b.len() == 1 => {
            let c1 = a.chars().next().unwrap();
            let c2 = b.chars().next().unwrap();
            if c1 == '0' || c2 == '0' || c1 == '-' || c2 == '-' {
                if c1 != '0' && c1 != '-' {
                    Genotype::Homozygous(c1)
                } else if c2 != '0' && c2 != '-' {
                    Genotype::Homozygous(c2)
                } else {
                    Genotype::NoCall
                }
            } else if c1 == c2 {
                Genotype::Homozygous(c1)
            } else {
                Genotype::Heterozygous(c1, c2)
            }
        }
        _ => Genotype::NoCall,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_ancestry() {
        let content = "\
#AncestryDNA raw data
rsid\tchromosome\tposition\tallele1\tallele2
rs4477212\t1\t82154\tA\tA
rs3094315\t1\t752566\tA\tG
rs9999999\t1\t800000\t0\t0
";
        let variants = parse(content).unwrap();
        assert_eq!(variants.len(), 3);
        assert_eq!(variants[0].genotype, Genotype::Homozygous('A'));
        assert_eq!(variants[0].source_format, SourceFormat::AncestryDNA);
        assert_eq!(variants[1].genotype, Genotype::Heterozygous('A', 'G'));
        assert_eq!(variants[2].genotype, Genotype::NoCall);
    }
}
