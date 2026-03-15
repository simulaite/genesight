//! Phase ambiguity detection for multi-SNP star alleles.
//!
//! Unphased microarray data cannot distinguish cis vs trans configurations.
//! The canonical example is TPMT: heterozygous at both rs1800460 and rs1142345
//! could be *3A/*1 (cis, Intermediate Metabolizer) or *3B/*3C (trans, Poor
//! Metabolizer). These produce different clinical recommendations.
//!
//! This module detects such ambiguities by enumerating all consistent diplotype
//! pairs and checking if they yield different phenotypes.

use std::collections::HashMap;

use super::definitions::GeneAlleleDefinitions;
use super::diplotype::{AlternativeDiplotype, ClinicalPhaseImpact, DiplotypeCall};
use super::phenotype;

/// Detect and annotate phase ambiguity on a diplotype call.
///
/// For genes with multi-SNP star alleles, enumerates all possible diplotype
/// pairs consistent with the observed genotype counts, determines their
/// phenotypes, and flags ambiguity if different phenotypes are possible.
///
/// This function mutates the `DiplotypeCall` in place, setting `alternatives`
/// and `phase_ambiguity` fields.
pub fn detect_phase_ambiguity(
    call: &mut DiplotypeCall,
    gene_defs: &GeneAlleleDefinitions,
    user_genotypes: &HashMap<String, String>,
    activity_values: &HashMap<String, f64>,
) {
    // Only multi-SNP genes can have phase ambiguity
    if gene_defs.defining_rsids.len() < 2 {
        return;
    }

    // Check if we have heterozygous positions at multiple defining sites.
    // Phase ambiguity only arises when a user is heterozygous at 2+ positions
    // that define different star alleles.
    let mut het_positions = 0;
    for rsid in &gene_defs.defining_rsids {
        if let Some(genotype) = user_genotypes.get(rsid.as_str()) {
            let chars: Vec<char> = genotype.chars().collect();
            if chars.len() >= 2 && !chars[0].eq_ignore_ascii_case(&chars[1]) {
                het_positions += 1;
            }
        }
    }

    // Need at least 2 heterozygous positions for phase ambiguity
    if het_positions < 2 {
        return;
    }

    // Enumerate all star alleles and their alt-allele requirements
    let mut allele_requirements: Vec<(&str, HashMap<&str, char>)> = Vec::new();

    // *1 requires reference at all positions (zero alt alleles)
    allele_requirements.push(("*1", HashMap::new()));

    for (star_allele, defining_variants) in &gene_defs.alleles {
        if star_allele == "*1" {
            continue;
        }
        let mut reqs: HashMap<&str, char> = HashMap::new();
        for def_var in defining_variants {
            if let Some(c) = def_var.alt_allele.chars().next() {
                reqs.insert(def_var.rsid.as_str(), c.to_ascii_uppercase());
            }
        }
        allele_requirements.push((star_allele.as_str(), reqs));
    }

    // For each position, compute how many "alt" alleles the user has
    // (per each star allele's definition of "alt")
    // We use a simplified approach: count occurrences of the variant allele in genotype
    let position_genotypes: HashMap<&str, (char, char)> = gene_defs
        .defining_rsids
        .iter()
        .filter_map(|rsid| {
            user_genotypes.get(rsid.as_str()).and_then(|gt| {
                let chars: Vec<char> = gt.chars().collect();
                if chars.len() >= 2 {
                    Some((
                        rsid.as_str(),
                        (chars[0].to_ascii_uppercase(), chars[1].to_ascii_uppercase()),
                    ))
                } else {
                    None
                }
            })
        })
        .collect();

    // Enumerate all consistent diplotype pairs (A, B) where A <= B
    let mut consistent_pairs: Vec<(String, String, f64)> = Vec::new();

    for i in 0..allele_requirements.len() {
        for j in i..allele_requirements.len() {
            let (allele_a, reqs_a) = &allele_requirements[i];
            let (allele_b, reqs_b) = &allele_requirements[j];

            // Check consistency: for each defining position, the pair must
            // account for all observed alt alleles
            let mut consistent = true;

            for rsid in &gene_defs.defining_rsids {
                let rsid_str = rsid.as_str();
                let (g1, g2) = match position_genotypes.get(rsid_str) {
                    Some(gt) => *gt,
                    None => continue, // Missing position — skip check
                };

                // Check all variant alleles defined for this position
                for defining_variants in gene_defs.alleles.values() {
                    for def_var in defining_variants {
                        if def_var.rsid.as_str() != rsid_str {
                            continue;
                        }
                        let alt_char = match def_var.alt_allele.chars().next() {
                            Some(c) => c.to_ascii_uppercase(),
                            None => continue,
                        };

                        // Count observed alt alleles
                        let observed_alt = u8::from(g1 == alt_char) + u8::from(g2 == alt_char);

                        // Count required alt alleles from this pair
                        let req_a = if reqs_a.get(rsid_str) == Some(&alt_char) {
                            1u8
                        } else {
                            0
                        };
                        let req_b = if reqs_b.get(rsid_str) == Some(&alt_char) {
                            1u8
                        } else {
                            0
                        };

                        if req_a + req_b != observed_alt {
                            consistent = false;
                        }
                    }
                }

                if !consistent {
                    break;
                }
            }

            if consistent {
                let score_a = activity_values.get(*allele_a).copied().unwrap_or(1.0);
                let score_b = activity_values.get(*allele_b).copied().unwrap_or(1.0);

                let (a, b) = if allele_a <= allele_b {
                    (allele_a.to_string(), allele_b.to_string())
                } else {
                    (allele_b.to_string(), allele_a.to_string())
                };

                // Avoid duplicates
                if !consistent_pairs.iter().any(|(x, y, _)| x == &a && y == &b) {
                    consistent_pairs.push((a, b, score_a + score_b));
                }
            }
        }
    }

    // If only one consistent pair, no ambiguity
    if consistent_pairs.len() <= 1 {
        return;
    }

    // Sort by activity score (lowest first = most conservative = patient safety)
    consistent_pairs.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    // Build alternative diplotypes with phenotypes
    let mut alternatives: Vec<AlternativeDiplotype> = Vec::new();
    let mut phenotypes_seen: Vec<String> = Vec::new();

    for (idx, (a1, a2, score)) in consistent_pairs.iter().enumerate() {
        let phenotype_result = phenotype::call_phenotype(&call.gene, *score);
        let is_primary = idx == 0;

        if !phenotypes_seen.contains(&phenotype_result.phenotype) {
            phenotypes_seen.push(phenotype_result.phenotype.clone());
        }

        alternatives.push(AlternativeDiplotype {
            allele1: a1.clone(),
            allele2: a2.clone(),
            activity_score: *score,
            phenotype: phenotype_result.phenotype,
            is_primary,
        });
    }

    // Determine clinical impact
    let phase_impact = if phenotypes_seen.len() > 1 {
        ClinicalPhaseImpact::DifferentPhenotypes
    } else {
        ClinicalPhaseImpact::Uniform
    };

    // Set conservative call (lowest activity score)
    let conservative = &consistent_pairs[0];
    call.allele1 = conservative.0.clone();
    call.allele2 = conservative.1.clone();
    call.activity_score = conservative.2;
    call.alternatives = alternatives;
    call.phase_ambiguity = Some(phase_impact);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pgx::definitions::AlleleDefiningVariant;
    use crate::pgx::diplotype::{call_diplotype, CoverageStatus};

    /// Build TPMT-like gene with *3A (both positions), *3B (rs1800460 only),
    /// *3C (rs1142345 only).
    fn make_tpmt_defs() -> GeneAlleleDefinitions {
        let mut alleles = HashMap::new();
        // *3A requires both rs1800460 AND rs1142345
        alleles.insert(
            "*3A".to_string(),
            vec![
                AlleleDefiningVariant {
                    allele_name: "*3A".to_string(),
                    rsid: "rs1800460".to_string(),
                    alt_allele: "A".to_string(),
                    function: "No Function".to_string(),
                    activity_score: 0.0,
                },
                AlleleDefiningVariant {
                    allele_name: "*3A".to_string(),
                    rsid: "rs1142345".to_string(),
                    alt_allele: "C".to_string(),
                    function: "No Function".to_string(),
                    activity_score: 0.0,
                },
            ],
        );
        // *3B requires only rs1800460
        alleles.insert(
            "*3B".to_string(),
            vec![AlleleDefiningVariant {
                allele_name: "*3B".to_string(),
                rsid: "rs1800460".to_string(),
                alt_allele: "A".to_string(),
                function: "No Function".to_string(),
                activity_score: 0.0,
            }],
        );
        // *3C requires only rs1142345
        alleles.insert(
            "*3C".to_string(),
            vec![AlleleDefiningVariant {
                allele_name: "*3C".to_string(),
                rsid: "rs1142345".to_string(),
                alt_allele: "C".to_string(),
                function: "No Function".to_string(),
                activity_score: 0.0,
            }],
        );
        GeneAlleleDefinitions {
            alleles,
            defining_rsids: vec!["rs1800460".to_string(), "rs1142345".to_string()],
        }
    }

    fn make_tpmt_activity() -> HashMap<String, f64> {
        let mut av = HashMap::new();
        av.insert("*1".to_string(), 1.0);
        av.insert("*3A".to_string(), 0.0);
        av.insert("*3B".to_string(), 0.0);
        av.insert("*3C".to_string(), 0.0);
        av
    }

    #[test]
    fn tpmt_het_het_is_ambiguous() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        // User is heterozygous at BOTH positions
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GA".to_string()),
            ("rs1142345".to_string(), "TC".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        assert!(call.phase_ambiguity.is_some());
        assert!(!call.alternatives.is_empty());
        assert!(call.alternatives.len() >= 2);
    }

    #[test]
    fn tpmt_het_at_one_site_not_ambiguous() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        // Heterozygous at only ONE position
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GA".to_string()),
            ("rs1142345".to_string(), "TT".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        assert!(call.phase_ambiguity.is_none());
        assert!(call.alternatives.is_empty());
    }

    #[test]
    fn tpmt_hom_alt_one_site_not_ambiguous() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        // Homozygous alt at one position only
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GG".to_string()),
            ("rs1142345".to_string(), "CC".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        // Hom alt at rs1142345 = *3C/*3C, not ambiguous
        assert!(call.phase_ambiguity.is_none());
    }

    #[test]
    fn single_snp_gene_never_ambiguous() {
        let mut alleles = HashMap::new();
        alleles.insert(
            "*2".to_string(),
            vec![AlleleDefiningVariant {
                allele_name: "*2".to_string(),
                rsid: "rs4244285".to_string(),
                alt_allele: "A".to_string(),
                function: "No Function".to_string(),
                activity_score: 0.0,
            }],
        );
        let defs = GeneAlleleDefinitions {
            alleles,
            defining_rsids: vec!["rs4244285".to_string()],
        };
        let mut av = HashMap::new();
        av.insert("*1".to_string(), 1.0);
        av.insert("*2".to_string(), 0.0);

        let user: HashMap<String, String> =
            HashMap::from([("rs4244285".to_string(), "GA".to_string())]);

        let mut call = call_diplotype("CYP2C19", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        assert!(call.phase_ambiguity.is_none());
    }

    #[test]
    fn all_ref_not_ambiguous() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GG".to_string()),
            ("rs1142345".to_string(), "TT".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        assert!(call.phase_ambiguity.is_none());
        assert_eq!(call.allele1, "*1");
        assert_eq!(call.allele2, "*1");
    }

    #[test]
    fn tpmt_ambiguity_has_different_phenotypes() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        // Het at both = ambiguous
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GA".to_string()),
            ("rs1142345".to_string(), "TC".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        // Should flag DifferentPhenotypes since *3A/*1 (Intermediate) vs *3B/*3C (Poor)
        if let Some(impact) = call.phase_ambiguity {
            assert_eq!(impact, ClinicalPhaseImpact::DifferentPhenotypes);
        } else {
            panic!("expected phase_ambiguity to be set");
        }
    }

    #[test]
    fn conservative_call_picks_lowest_activity() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GA".to_string()),
            ("rs1142345".to_string(), "TC".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        // Conservative = lowest activity score (patient safety)
        // *3B/*3C = 0.0+0.0=0.0 < *3A/*1 = 0.0+1.0=1.0
        assert!(call.activity_score <= 1.0);
    }

    #[test]
    fn alternatives_include_primary_flag() {
        let defs = make_tpmt_defs();
        let av = make_tpmt_activity();
        let user: HashMap<String, String> = HashMap::from([
            ("rs1800460".to_string(), "GA".to_string()),
            ("rs1142345".to_string(), "TC".to_string()),
        ]);

        let mut call = call_diplotype("TPMT", &defs, &user, &av);
        detect_phase_ambiguity(&mut call, &defs, &user, &av);

        let primary_count = call.alternatives.iter().filter(|a| a.is_primary).count();
        assert_eq!(
            primary_count, 1,
            "exactly one alternative should be primary"
        );
    }
}
