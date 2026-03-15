//! Diplotype calling algorithm for pharmacogenomics.
//!
//! Determines which star alleles a user carries for a given gene based on
//! their genotype at defining SNP positions. Since consumer array data is
//! unphased, we use the most parsimonious diplotype assignment.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::allele::{match_single_allele, AlleleMatch};

use super::definitions::GeneAlleleDefinitions;

/// Coverage status of PGx defining positions in the user's data.
///
/// When a defining SNP for a star allele is absent from the user's array data,
/// the pipeline must NOT silently default to *1 (wildtype). This enum tracks
/// how much of the gene's defining positions were covered.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum CoverageStatus {
    /// All defining positions for this gene were observed in the user's data.
    Complete,
    /// Some defining positions are missing. The diplotype call is tentative.
    Partial {
        /// rsIDs not present in the user's data.
        missing: Vec<String>,
        /// Fraction of positions observed (0.0 to 1.0).
        coverage_pct: f64,
    },
    /// Too few defining positions were observed to make any meaningful call.
    /// The phenotype should be reported as "Indeterminate".
    Insufficient {
        /// rsIDs not present in the user's data.
        missing: Vec<String>,
    },
}

/// Minimum fraction of defining positions required for a partial call.
/// Below this threshold, coverage is `Insufficient`.
const MINIMUM_COVERAGE: f64 = 0.5;

/// Result of calling a diplotype for a gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiplotypeCall {
    /// Gene symbol (e.g., "CYP2C19")
    pub gene: String,
    /// First allele (e.g., "*1")
    pub allele1: String,
    /// Second allele (e.g., "*17")
    pub allele2: String,
    /// Sum of activity values for both alleles
    pub activity_score: f64,
    /// Coverage status of defining positions.
    pub coverage: CoverageStatus,
    /// Alternative diplotype interpretations due to phasing ambiguity.
    /// Empty if the call is unambiguous.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub alternatives: Vec<AlternativeDiplotype>,
    /// Whether phasing ambiguity affects clinical interpretation.
    pub phase_ambiguity: Option<ClinicalPhaseImpact>,
}

/// An alternative diplotype interpretation from phase ambiguity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlternativeDiplotype {
    /// First allele
    pub allele1: String,
    /// Second allele
    pub allele2: String,
    /// Activity score for this interpretation
    pub activity_score: f64,
    /// Phenotype for this interpretation
    pub phenotype: String,
    /// Whether this is the primary (conservative) call
    pub is_primary: bool,
}

/// Whether phasing ambiguity has clinical impact.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ClinicalPhaseImpact {
    /// All alternative interpretations yield the same phenotype.
    Uniform,
    /// Alternative interpretations yield different phenotypes — confirmation needed.
    DifferentPhenotypes,
}

impl DiplotypeCall {
    /// Format the diplotype as "allele1/allele2" (e.g., "*1/*17").
    pub fn diplotype_string(&self) -> String {
        format!("{}/{}", self.allele1, self.allele2)
    }

    /// Whether coverage is complete (all defining positions observed).
    pub fn is_complete(&self) -> bool {
        matches!(self.coverage, CoverageStatus::Complete)
    }

    /// Whether the call has insufficient coverage.
    pub fn is_insufficient(&self) -> bool {
        matches!(self.coverage, CoverageStatus::Insufficient { .. })
    }
}

/// Call the diplotype for a gene given user genotype data.
///
/// Algorithm (simplified PharmCAT approach for array data):
/// 1. For each star allele definition, check if the user has the defining variant(s)
/// 2. A star allele is "called" if the user has alt alleles at ALL its defining positions
/// 3. Since data is unphased, construct possible diplotype pairs
/// 4. Pick the most parsimonious diplotype (fewest rare alleles)
/// 5. If core defining SNPs are missing, flag with limitation
///
/// # Arguments
///
/// * `gene` - Gene symbol (e.g., "CYP2C19")
/// * `gene_defs` - Star allele definitions for this gene
/// * `user_genotypes` - Maps rsID -> observed genotype string (e.g., "AG", "CC")
/// * `activity_values` - Maps star allele name -> activity value (e.g., "*2" -> 0.0)
pub fn call_diplotype(
    gene: &str,
    gene_defs: &GeneAlleleDefinitions,
    user_genotypes: &HashMap<String, String>,
    activity_values: &HashMap<String, f64>,
) -> DiplotypeCall {
    // Compute coverage status
    let total_positions = gene_defs.defining_rsids.len();
    let missing_positions: Vec<String> = gene_defs
        .defining_rsids
        .iter()
        .filter(|rsid| !user_genotypes.contains_key(rsid.as_str()))
        .cloned()
        .collect();

    let coverage = if missing_positions.is_empty() {
        CoverageStatus::Complete
    } else if total_positions == 0 {
        CoverageStatus::Insufficient {
            missing: missing_positions.clone(),
        }
    } else {
        let observed = total_positions - missing_positions.len();
        let coverage_pct = observed as f64 / total_positions as f64;
        if coverage_pct >= MINIMUM_COVERAGE {
            CoverageStatus::Partial {
                missing: missing_positions.clone(),
                coverage_pct,
            }
        } else {
            CoverageStatus::Insufficient {
                missing: missing_positions.clone(),
            }
        }
    };

    // Determine which star alleles the user could carry
    let mut candidate_alleles: Vec<(String, f64)> = Vec::new();

    for (star_allele, defining_variants) in &gene_defs.alleles {
        if star_allele == "*1" {
            continue; // *1 is the default/reference allele
        }

        // Check if user has the defining variant(s) for this star allele
        let mut all_match = true;
        let mut any_present = false;
        let mut total_alt_copies: usize = 0;

        for def_var in defining_variants {
            match user_genotypes.get(&def_var.rsid) {
                Some(genotype) => {
                    any_present = true;
                    let alt_char = match def_var.alt_allele.chars().next() {
                        Some(c) => c.to_ascii_uppercase(),
                        None => {
                            all_match = false;
                            continue;
                        }
                    };
                    let count = genotype
                        .chars()
                        .filter(|&c| {
                            let m = match_single_allele(alt_char, c);
                            matches!(m, AlleleMatch::DirectMatch | AlleleMatch::ComplementMatch)
                        })
                        .count();
                    if count == 0 {
                        // User is homozygous reference at this position
                        all_match = false;
                    }
                    total_alt_copies += count;
                }
                None => {
                    // Position missing from user data — can't confirm this allele
                    all_match = false;
                }
            }
        }

        if any_present && all_match {
            let activity = activity_values
                .get(star_allele.as_str())
                .copied()
                .unwrap_or(1.0);

            // If user has 2 copies of alt at all defining positions,
            // this allele appears on both chromosomes
            if total_alt_copies >= 2 * defining_variants.len() {
                candidate_alleles.push((star_allele.clone(), activity));
                candidate_alleles.push((star_allele.clone(), activity));
            } else {
                candidate_alleles.push((star_allele.clone(), activity));
            }
        }
    }

    // Build the diplotype
    // Most parsimonious: the remaining chromosome gets *1
    let (allele1, allele2, activity_score) = match candidate_alleles.len() {
        0 => {
            let a1_act = activity_values.get("*1").copied().unwrap_or(1.0);
            ("*1".to_string(), "*1".to_string(), a1_act * 2.0)
        }
        1 => {
            let a1_act = activity_values.get("*1").copied().unwrap_or(1.0);
            let a2_act = candidate_alleles[0].1;
            (
                candidate_alleles[0].0.clone(),
                "*1".to_string(),
                a1_act + a2_act,
            )
        }
        _ => {
            // Sort by activity (rarest/most impactful first)
            candidate_alleles
                .sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            candidate_alleles.truncate(2);
            let a1 = &candidate_alleles[0];
            let a2 = &candidate_alleles[1];
            (a1.0.clone(), a2.0.clone(), a1.1 + a2.1)
        }
    };

    // Normalize order: alphabetically smaller allele first
    let (allele1, allele2) = if allele1 <= allele2 {
        (allele1, allele2)
    } else {
        (allele2, allele1)
    };

    DiplotypeCall {
        gene: gene.to_string(),
        allele1,
        allele2,
        activity_score,
        coverage,
        alternatives: Vec::new(),
        phase_ambiguity: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pgx::definitions::AlleleDefiningVariant;

    fn make_gene_defs() -> GeneAlleleDefinitions {
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
        alleles.insert(
            "*17".to_string(),
            vec![AlleleDefiningVariant {
                allele_name: "*17".to_string(),
                rsid: "rs12248560".to_string(),
                alt_allele: "T".to_string(),
                function: "Increased Function".to_string(),
                activity_score: 1.5,
            }],
        );
        GeneAlleleDefinitions {
            alleles,
            defining_rsids: vec!["rs4244285".to_string(), "rs12248560".to_string()],
        }
    }

    fn make_activity_values() -> HashMap<String, f64> {
        let mut av = HashMap::new();
        av.insert("*1".to_string(), 1.0);
        av.insert("*2".to_string(), 0.0);
        av.insert("*17".to_string(), 1.5);
        av
    }

    #[test]
    fn wildtype_when_no_variants() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        let user: HashMap<String, String> = HashMap::from([
            ("rs4244285".to_string(), "GG".to_string()),
            ("rs12248560".to_string(), "CC".to_string()),
        ]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert_eq!(call.allele1, "*1");
        assert_eq!(call.allele2, "*1");
        assert!((call.activity_score - 2.0).abs() < f64::EPSILON);
        assert!(call.is_complete());
        assert_eq!(call.coverage, CoverageStatus::Complete);
    }

    #[test]
    fn heterozygous_star2() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        let user: HashMap<String, String> = HashMap::from([
            ("rs4244285".to_string(), "GA".to_string()),
            ("rs12248560".to_string(), "CC".to_string()),
        ]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert!(call.allele1 == "*1" || call.allele2 == "*1");
        assert!(call.allele1 == "*2" || call.allele2 == "*2");
        assert!((call.activity_score - 1.0).abs() < f64::EPSILON); // *1(1.0) + *2(0.0)
    }

    #[test]
    fn rapid_metabolizer_star17() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        let user: HashMap<String, String> = HashMap::from([
            ("rs4244285".to_string(), "GG".to_string()),
            ("rs12248560".to_string(), "CT".to_string()),
        ]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert!(call.allele1 == "*1" || call.allele2 == "*1");
        assert!(call.allele1 == "*17" || call.allele2 == "*17");
        assert!((call.activity_score - 2.5).abs() < f64::EPSILON); // *1(1.0) + *17(1.5)
    }

    #[test]
    fn missing_positions_flagged_partial() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        // Only one of two positions present -> 50% coverage = Partial
        let user: HashMap<String, String> =
            HashMap::from([("rs4244285".to_string(), "GG".to_string())]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert!(!call.is_complete());
        match &call.coverage {
            CoverageStatus::Partial {
                missing,
                coverage_pct,
            } => {
                assert!(missing.contains(&"rs12248560".to_string()));
                assert!((*coverage_pct - 0.5).abs() < f64::EPSILON);
            }
            other => panic!("expected Partial, got {other:?}"),
        }
    }

    #[test]
    fn all_positions_missing_is_insufficient() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        // No positions present at all
        let user: HashMap<String, String> = HashMap::new();

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert!(call.is_insufficient());
        match &call.coverage {
            CoverageStatus::Insufficient { missing } => {
                assert_eq!(missing.len(), 2);
            }
            other => panic!("expected Insufficient, got {other:?}"),
        }
    }

    #[test]
    fn ref_genotype_not_called_ultrarapid() {
        // Critical test: user with rs12248560 CC (ref allele) should NOT be *17
        let defs = make_gene_defs();
        let av = make_activity_values();
        let user: HashMap<String, String> = HashMap::from([
            ("rs4244285".to_string(), "GG".to_string()),
            ("rs12248560".to_string(), "CC".to_string()), // CC = reference
        ]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert_eq!(call.allele1, "*1");
        assert_eq!(call.allele2, "*1");
        // Should NOT be *17 — the *17 allele is T, not C
    }

    #[test]
    fn homozygous_alt_calls_two_copies() {
        let defs = make_gene_defs();
        let av = make_activity_values();
        let user: HashMap<String, String> = HashMap::from([
            ("rs4244285".to_string(), "AA".to_string()), // homozygous *2
            ("rs12248560".to_string(), "CC".to_string()),
        ]);

        let call = call_diplotype("CYP2C19", &defs, &user, &av);
        assert_eq!(call.allele1, "*2");
        assert_eq!(call.allele2, "*2");
        assert!((call.activity_score - 0.0).abs() < f64::EPSILON); // *2(0.0) + *2(0.0)
    }

    #[test]
    fn diplotype_string_format() {
        let call = DiplotypeCall {
            gene: "CYP2C19".to_string(),
            allele1: "*1".to_string(),
            allele2: "*17".to_string(),
            activity_score: 2.5,
            coverage: CoverageStatus::Complete,
            alternatives: vec![],
            phase_ambiguity: None,
        };
        assert_eq!(call.diplotype_string(), "*1/*17");
    }
}
