//! Gene-specific phenotype classification from activity scores.
//!
//! Each gene has different phenotype thresholds and nomenclature based on
//! CPIC guidelines. This module translates numeric activity scores into
//! standard phenotype names and attaches gene-specific limitation notes.

/// Result of phenotype classification for a gene.
#[derive(Debug, Clone)]
pub struct PhenotypeResult {
    /// Classified phenotype name (e.g., "Poor Metabolizer")
    pub phenotype: String,
    /// Limitations or caveats about this classification
    pub limitations: Vec<String>,
}

/// Classify a phenotype from an activity score, with gene-specific logic.
///
/// Combines activity-score-based classification with gene-specific limitation
/// notes about phasing, single-variant coverage, and other caveats.
///
/// # Arguments
///
/// * `gene` - Gene symbol (e.g., "CYP2D6", "DPYD")
/// * `activity_score` - Combined diplotype activity score (sum of both alleles)
///
/// # Returns
///
/// A `PhenotypeResult` with the phenotype name and any applicable limitations.
pub fn call_phenotype(gene: &str, activity_score: f64) -> PhenotypeResult {
    let phenotype = classify_by_activity_score(gene, activity_score);
    let mut limitations = common_limitations();

    match gene {
        "TPMT" => {
            limitations.push(
                "TPMT *3A requires both rs1800460 and rs1142345 variants in cis (on the \
                 same chromosome). Without phasing data, *3A cannot be distinguished from \
                 compound heterozygous *3B/*3C. This may affect phenotype classification."
                    .to_string(),
            );
        }
        "VKORC1" => {
            limitations.push(
                "VKORC1 warfarin sensitivity is based on a single promoter variant \
                 (rs9923231, -1639G>A). Other factors including CYP2C9 genotype, age, \
                 weight, diet (vitamin K intake), and concomitant medications substantially \
                 influence warfarin dose requirements."
                    .to_string(),
            );
        }
        _ => {}
    }

    PhenotypeResult {
        phenotype: phenotype.to_string(),
        limitations,
    }
}

/// Classify phenotype from activity score with gene-specific thresholds.
///
/// Each gene uses different phenotype names and score boundaries per CPIC
/// guidelines. Falls back to generic metabolizer categories for unknown genes.
fn classify_by_activity_score(gene: &str, score: f64) -> &'static str {
    match gene {
        "CYP2C19" => {
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 1.5 {
                "Intermediate Metabolizer"
            } else if score <= 2.0 {
                "Normal Metabolizer"
            } else if score < 3.0 {
                "Rapid Metabolizer"
            } else {
                "Ultrarapid Metabolizer"
            }
        }
        "CYP2D6" => {
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 1.5 {
                "Intermediate Metabolizer"
            } else if score <= 2.0 {
                "Normal Metabolizer"
            } else {
                "Ultrarapid Metabolizer"
            }
        }
        "CYP2C9" => {
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 2.0 {
                "Intermediate Metabolizer"
            } else {
                "Normal Metabolizer"
            }
        }
        "SLCO1B1" => {
            if score <= 1.0 {
                "Poor Function"
            } else if score < 2.0 {
                "Intermediate Function"
            } else {
                "Normal Function"
            }
        }
        "CYP3A5" => {
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 2.0 {
                "Intermediate Metabolizer"
            } else {
                "Extensive Metabolizer"
            }
        }
        "DPYD" => {
            if score <= 0.0 {
                "Poor DPD Activity (DPD Deficient)"
            } else if score < 2.0 {
                "Intermediate DPD Activity"
            } else {
                "Normal DPD Activity"
            }
        }
        "TPMT" | "NUDT15" => {
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 2.0 {
                "Intermediate Metabolizer"
            } else {
                "Normal Metabolizer"
            }
        }
        "VKORC1" => {
            if score >= 2.0 {
                "Normal Warfarin Sensitivity"
            } else if score >= 1.5 {
                "Increased Warfarin Sensitivity"
            } else {
                "Highly Increased Warfarin Sensitivity"
            }
        }
        _ => {
            // Generic fallback for unknown genes
            if score <= 0.0 {
                "Poor Metabolizer"
            } else if score < 2.0 {
                "Intermediate Metabolizer"
            } else {
                "Normal Metabolizer"
            }
        }
    }
}

/// Common limitations that apply to all PGx calls from microarray data.
fn common_limitations() -> Vec<String> {
    vec![
        "Microarray genotyping does not detect gene deletions or duplications (copy number variants).".to_string(),
        "Diplotype inference from unphased data is approximate; clinical-grade testing recommended for actionable results.".to_string(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- CYP2C19 ----

    #[test]
    fn cyp2c19_poor_metabolizer_at_zero() {
        let result = call_phenotype("CYP2C19", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn cyp2c19_intermediate_at_one() {
        let result = call_phenotype("CYP2C19", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn cyp2c19_normal_at_two() {
        let result = call_phenotype("CYP2C19", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    #[test]
    fn cyp2c19_rapid_at_two_point_five() {
        let result = call_phenotype("CYP2C19", 2.5);
        assert_eq!(result.phenotype, "Rapid Metabolizer");
    }

    #[test]
    fn cyp2c19_ultrarapid_at_three() {
        let result = call_phenotype("CYP2C19", 3.0);
        assert_eq!(result.phenotype, "Ultrarapid Metabolizer");
    }

    // ---- CYP2D6 ----

    #[test]
    fn cyp2d6_poor_metabolizer_at_zero() {
        let result = call_phenotype("CYP2D6", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn cyp2d6_intermediate_at_one() {
        let result = call_phenotype("CYP2D6", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn cyp2d6_normal_at_two() {
        let result = call_phenotype("CYP2D6", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    #[test]
    fn cyp2d6_ultrarapid_above_two() {
        let result = call_phenotype("CYP2D6", 2.5);
        assert_eq!(result.phenotype, "Ultrarapid Metabolizer");
    }

    // ---- CYP2C9 ----

    #[test]
    fn cyp2c9_poor_metabolizer_at_zero() {
        let result = call_phenotype("CYP2C9", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn cyp2c9_intermediate_at_one() {
        let result = call_phenotype("CYP2C9", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn cyp2c9_normal_at_two() {
        let result = call_phenotype("CYP2C9", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    // ---- SLCO1B1 ----

    #[test]
    fn slco1b1_poor_function_at_one() {
        let result = call_phenotype("SLCO1B1", 1.0);
        assert_eq!(result.phenotype, "Poor Function");
    }

    #[test]
    fn slco1b1_intermediate_at_one_point_five() {
        let result = call_phenotype("SLCO1B1", 1.5);
        assert_eq!(result.phenotype, "Intermediate Function");
    }

    #[test]
    fn slco1b1_normal_at_two() {
        let result = call_phenotype("SLCO1B1", 2.0);
        assert_eq!(result.phenotype, "Normal Function");
    }

    // ---- CYP3A5 ----

    #[test]
    fn cyp3a5_poor_metabolizer_at_zero() {
        let result = call_phenotype("CYP3A5", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn cyp3a5_intermediate_at_one() {
        let result = call_phenotype("CYP3A5", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn cyp3a5_extensive_at_two() {
        let result = call_phenotype("CYP3A5", 2.0);
        assert_eq!(result.phenotype, "Extensive Metabolizer");
    }

    // ---- DPYD ----

    #[test]
    fn dpyd_poor_activity_at_zero() {
        let result = call_phenotype("DPYD", 0.0);
        assert_eq!(result.phenotype, "Poor DPD Activity (DPD Deficient)");
    }

    #[test]
    fn dpyd_intermediate_at_one() {
        let result = call_phenotype("DPYD", 1.0);
        assert_eq!(result.phenotype, "Intermediate DPD Activity");
    }

    #[test]
    fn dpyd_normal_at_two() {
        let result = call_phenotype("DPYD", 2.0);
        assert_eq!(result.phenotype, "Normal DPD Activity");
    }

    // ---- TPMT ----

    #[test]
    fn tpmt_poor_metabolizer_at_zero() {
        let result = call_phenotype("TPMT", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn tpmt_intermediate_at_one() {
        let result = call_phenotype("TPMT", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn tpmt_normal_at_two() {
        let result = call_phenotype("TPMT", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    #[test]
    fn tpmt_has_phasing_limitation() {
        let result = call_phenotype("TPMT", 1.0);
        assert!(result.limitations.iter().any(|l| l.contains("*3A")));
        assert!(result.limitations.iter().any(|l| l.contains("phasing")));
    }

    // ---- NUDT15 ----

    #[test]
    fn nudt15_poor_metabolizer_at_zero() {
        let result = call_phenotype("NUDT15", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");
    }

    #[test]
    fn nudt15_intermediate_at_one() {
        let result = call_phenotype("NUDT15", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn nudt15_normal_at_two() {
        let result = call_phenotype("NUDT15", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    // ---- VKORC1 ----

    #[test]
    fn vkorc1_normal_sensitivity_at_two() {
        let result = call_phenotype("VKORC1", 2.0);
        assert_eq!(result.phenotype, "Normal Warfarin Sensitivity");
    }

    #[test]
    fn vkorc1_increased_sensitivity_at_one_point_five() {
        let result = call_phenotype("VKORC1", 1.5);
        assert_eq!(result.phenotype, "Increased Warfarin Sensitivity");
    }

    #[test]
    fn vkorc1_highly_increased_at_one() {
        let result = call_phenotype("VKORC1", 1.0);
        assert_eq!(result.phenotype, "Highly Increased Warfarin Sensitivity");
    }

    #[test]
    fn vkorc1_has_single_variant_limitation() {
        let result = call_phenotype("VKORC1", 1.5);
        assert!(result
            .limitations
            .iter()
            .any(|l| l.contains("single promoter variant")));
    }

    // ---- Generic / unknown gene ----

    #[test]
    fn unknown_gene_uses_generic_classification() {
        let result = call_phenotype("UNKNOWN_GENE", 0.0);
        assert_eq!(result.phenotype, "Poor Metabolizer");

        let result = call_phenotype("UNKNOWN_GENE", 1.0);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");

        let result = call_phenotype("UNKNOWN_GENE", 2.0);
        assert_eq!(result.phenotype, "Normal Metabolizer");
    }

    // ---- Common limitations ----

    #[test]
    fn all_results_have_common_limitations() {
        let result = call_phenotype("CYP2D6", 1.0);
        assert!(result.limitations.iter().any(|l| l.contains("copy number")));
        assert!(result.limitations.iter().any(|l| l.contains("unphased")));
    }

    // ---- Boundary tests ----

    #[test]
    fn boundary_cyp3a5_just_below_two() {
        let result = call_phenotype("CYP3A5", 1.999);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn boundary_dpyd_just_below_two() {
        let result = call_phenotype("DPYD", 1.999);
        assert_eq!(result.phenotype, "Intermediate DPD Activity");
    }

    #[test]
    fn boundary_tpmt_just_below_two() {
        let result = call_phenotype("TPMT", 1.999);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn boundary_nudt15_just_below_two() {
        let result = call_phenotype("NUDT15", 1.999);
        assert_eq!(result.phenotype, "Intermediate Metabolizer");
    }

    #[test]
    fn boundary_vkorc1_just_below_one_point_five() {
        let result = call_phenotype("VKORC1", 1.499);
        assert_eq!(result.phenotype, "Highly Increased Warfarin Sensitivity");
    }

    #[test]
    fn boundary_vkorc1_at_zero() {
        let result = call_phenotype("VKORC1", 0.0);
        assert_eq!(result.phenotype, "Highly Increased Warfarin Sensitivity");
    }
}
