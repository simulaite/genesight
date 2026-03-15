//! Allele matching and strand resolution for genotype comparisons.
//!
//! DNA genotyping arrays report alleles on an arbitrary strand. When comparing
//! a user's genotype against a database reference allele, we must account for
//! the possibility that the alleles are reported on opposite strands. For most
//! SNPs this is straightforward (complement the bases), but **palindromic SNPs**
//! (A/T and C/G) are inherently ambiguous because both strands yield the same
//! pair of bases.
//!
//! This module provides:
//! - [`match_alleles`] — strand-aware matching that flags palindromic SNPs
//! - [`match_alleles_with_frequency`] — uses allele frequency to resolve
//!   palindromic ambiguity when possible
//! - [`AlleleMatch`] — the result of an allele comparison

pub mod strand;

use strand::{complement, is_palindromic};

/// Result of comparing a user's genotype alleles against a database reference allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleMatch {
    /// The user carries the reference allele on the forward strand.
    DirectMatch,
    /// The user carries the complement of the reference allele, indicating
    /// the genotype was reported on the opposite strand.
    ComplementMatch,
    /// The SNP is palindromic (A/T or C/G) and strand cannot be determined
    /// from sequence alone. Frequency-based resolution may help.
    StrandAmbiguous,
    /// The user does not carry the reference allele on either strand.
    Mismatch,
}

/// Compare user genotype alleles against a database reference allele.
///
/// For non-palindromic SNPs, checks both direct and complement matching.
/// For palindromic SNPs (A/T, C/G), returns [`AlleleMatch::StrandAmbiguous`]
/// because the strand cannot be determined from sequence alone.
///
/// # Arguments
///
/// * `user_alleles` - The user's genotype alleles (e.g., `('A', 'G')`)
/// * `ref_allele` - The database reference/risk allele (e.g., `'A'`)
/// * `alt_allele` - The alternate allele at this position (e.g., `'G'`)
///
/// # Examples
///
/// ```
/// use genesight_core::allele::{match_alleles, AlleleMatch};
///
/// // Non-palindromic: A/G — clear strand resolution
/// assert_eq!(match_alleles(('A', 'G'), 'A', 'G'), AlleleMatch::DirectMatch);
///
/// // Palindromic: A/T — ambiguous without frequency data
/// assert_eq!(match_alleles(('A', 'T'), 'A', 'T'), AlleleMatch::StrandAmbiguous);
/// ```
pub fn match_alleles(
    user_alleles: (char, char),
    ref_allele: char,
    alt_allele: char,
) -> AlleleMatch {
    let (a1, a2) = (
        user_alleles.0.to_ascii_uppercase(),
        user_alleles.1.to_ascii_uppercase(),
    );
    let ref_upper = ref_allele.to_ascii_uppercase();
    let alt_upper = alt_allele.to_ascii_uppercase();

    // If the ref/alt pair is palindromic, strand is ambiguous
    if is_palindromic(ref_upper, alt_upper) {
        return AlleleMatch::StrandAmbiguous;
    }

    // Check direct match (same strand)
    if a1 == ref_upper || a2 == ref_upper {
        return AlleleMatch::DirectMatch;
    }

    // Check complement match (opposite strand)
    if let Some(ref_comp) = complement(ref_upper) {
        if a1 == ref_comp || a2 == ref_comp {
            return AlleleMatch::ComplementMatch;
        }
    }

    AlleleMatch::Mismatch
}

/// Threshold for allele frequency difference to resolve palindromic ambiguity.
///
/// If the absolute difference between the user-side frequency and the
/// database frequency exceeds this threshold, we consider the strand
/// resolved. A value of 0.10 (10%) provides reasonable confidence while
/// avoiding false resolution for common variants near 50% frequency.
const AF_RESOLUTION_THRESHOLD: f64 = 0.10;

/// Compare user alleles against a reference allele, using allele frequency
/// to resolve palindromic SNP ambiguity when possible.
///
/// For non-palindromic SNPs, this behaves identically to [`match_alleles`].
/// For palindromic SNPs, it attempts to resolve strand using the difference
/// between the user-observed allele frequency and the database allele
/// frequency:
///
/// - If the frequencies are close (within [`AF_RESOLUTION_THRESHOLD`]),
///   the alleles are likely on the same strand: returns `DirectMatch`.
/// - If the frequencies are far apart, strand is likely flipped: returns
///   `ComplementMatch`.
/// - If frequency data is insufficient, returns `StrandAmbiguous`.
///
/// # Arguments
///
/// * `user_alleles` - The user's genotype alleles
/// * `ref_allele` - Database reference/risk allele
/// * `alt_allele` - Alternate allele
/// * `user_af` - Observed frequency of the reference allele in the user's
///   population (if known). In single-sample mode this is typically `None`.
/// * `db_af` - Allele frequency from the database (gnomAD/dbSNP)
pub fn match_alleles_with_frequency(
    user_alleles: (char, char),
    ref_allele: char,
    alt_allele: char,
    user_af: Option<f64>,
    db_af: Option<f64>,
) -> AlleleMatch {
    let basic = match_alleles(user_alleles, ref_allele, alt_allele);

    // Only attempt frequency resolution for palindromic SNPs
    if basic != AlleleMatch::StrandAmbiguous {
        return basic;
    }

    // Need both frequencies to resolve
    let (u_af, d_af) = match (user_af, db_af) {
        (Some(u), Some(d)) => (u, d),
        // With only db_af, we can still heuristically resolve: if the db AF
        // is far from 0.5, the palindromic allele assignment is less ambiguous.
        // However, without a user AF to compare, we stay conservative.
        (None, Some(d)) => {
            // If AF is very far from 0.5, the chance of strand flip being
            // undetectable is low. Use a stricter threshold for single-sample.
            let dist_from_half = (d - 0.5).abs();
            if dist_from_half > 0.40 {
                // AF < 0.10 or > 0.90: very likely same strand
                return AlleleMatch::DirectMatch;
            }
            return AlleleMatch::StrandAmbiguous;
        }
        _ => return AlleleMatch::StrandAmbiguous,
    };

    let diff = (u_af - d_af).abs();

    if diff <= AF_RESOLUTION_THRESHOLD {
        // Frequencies agree: same strand
        AlleleMatch::DirectMatch
    } else if (u_af - (1.0 - d_af)).abs() <= AF_RESOLUTION_THRESHOLD {
        // User AF matches complement of db AF: flipped strand
        AlleleMatch::ComplementMatch
    } else {
        // Cannot confidently resolve
        AlleleMatch::StrandAmbiguous
    }
}

/// Result of counting risk allele copies in a user's genotype.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RiskAlleleCopies {
    /// Determined copy count (0, 1, or 2) with the match type used.
    Determined { copies: u8, match_type: AlleleMatch },
    /// The risk allele is palindromic with the other allele, so strand
    /// cannot be resolved from sequence alone.
    Palindromic { copies: u8 },
    /// Risk allele was not provided in the GWAS data.
    Indeterminate,
}

/// Match a single risk allele character against an observed allele.
///
/// Returns `DirectMatch` if the bases match (case-insensitive),
/// `ComplementMatch` if the complement matches, or `Mismatch`.
pub fn match_single_allele(risk: char, observed: char) -> AlleleMatch {
    let r = risk.to_ascii_uppercase();
    let o = observed.to_ascii_uppercase();

    if r == o {
        AlleleMatch::DirectMatch
    } else if complement(r) == Some(o) {
        AlleleMatch::ComplementMatch
    } else {
        AlleleMatch::Mismatch
    }
}

/// Count how many copies of a risk allele a genotype carries,
/// with strand-awareness and palindromic detection.
///
/// # Arguments
///
/// * `risk_allele` - The risk allele string from GWAS (e.g., "A", "G").
///   If `None`, returns `Indeterminate`.
/// * `allele1` - First allele of the user's genotype.
/// * `allele2` - Second allele of the user's genotype.
pub fn count_risk_allele_copies(
    risk_allele: Option<&str>,
    allele1: char,
    allele2: char,
) -> RiskAlleleCopies {
    let risk_str = match risk_allele {
        Some(r) if !r.is_empty() => r.trim(),
        _ => return RiskAlleleCopies::Indeterminate,
    };

    let risk_char = match risk_str.chars().next() {
        Some(c) if c.is_ascii_alphabetic() => c.to_ascii_uppercase(),
        _ => return RiskAlleleCopies::Indeterminate,
    };

    let a1 = allele1.to_ascii_uppercase();
    let a2 = allele2.to_ascii_uppercase();

    // Check if this is a palindromic situation
    let is_palindrome = if a1 == a2 {
        false // Homozygous — cannot determine palindrome from genotype alone
    } else {
        is_palindromic(a1, a2)
    };

    let match1 = match_single_allele(risk_char, a1);
    let match2 = match_single_allele(risk_char, a2);

    let direct_copies =
        u8::from(match1 == AlleleMatch::DirectMatch) + u8::from(match2 == AlleleMatch::DirectMatch);
    let complement_copies = u8::from(match1 == AlleleMatch::ComplementMatch)
        + u8::from(match2 == AlleleMatch::ComplementMatch);

    if is_palindrome {
        let copies = if direct_copies > 0 {
            direct_copies
        } else {
            complement_copies
        };
        RiskAlleleCopies::Palindromic { copies }
    } else if direct_copies > 0 {
        RiskAlleleCopies::Determined {
            copies: direct_copies,
            match_type: AlleleMatch::DirectMatch,
        }
    } else if complement_copies > 0 {
        RiskAlleleCopies::Determined {
            copies: complement_copies,
            match_type: AlleleMatch::ComplementMatch,
        }
    } else {
        RiskAlleleCopies::Determined {
            copies: 0,
            match_type: AlleleMatch::Mismatch,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- match_alleles tests ---

    #[test]
    fn non_palindromic_direct_match() {
        // A/G SNP, user has A — direct match
        assert_eq!(
            match_alleles(('A', 'G'), 'A', 'G'),
            AlleleMatch::DirectMatch
        );
    }

    #[test]
    fn non_palindromic_complement_match() {
        // A/G SNP but user reports on opposite strand (T/C), ref allele A
        // complement of A is T, which the user has
        assert_eq!(
            match_alleles(('T', 'C'), 'A', 'G'),
            AlleleMatch::ComplementMatch
        );
    }

    #[test]
    fn non_palindromic_mismatch() {
        // A/G SNP, user has C/T — neither matches A or complement T...
        // wait, complement of A is T and user has T. Let me pick better values.
        // ref=A, alt=G, user=(C,C): C is not A, complement(A)=T, C!=T => mismatch
        assert_eq!(match_alleles(('C', 'C'), 'A', 'G'), AlleleMatch::Mismatch);
    }

    #[test]
    fn palindromic_at_returns_ambiguous() {
        // A/T SNP — palindromic, cannot resolve strand
        assert_eq!(
            match_alleles(('A', 'T'), 'A', 'T'),
            AlleleMatch::StrandAmbiguous
        );
    }

    #[test]
    fn palindromic_ta_returns_ambiguous() {
        assert_eq!(
            match_alleles(('T', 'A'), 'T', 'A'),
            AlleleMatch::StrandAmbiguous
        );
    }

    #[test]
    fn palindromic_cg_returns_ambiguous() {
        // C/G SNP — palindromic
        assert_eq!(
            match_alleles(('C', 'G'), 'C', 'G'),
            AlleleMatch::StrandAmbiguous
        );
    }

    #[test]
    fn palindromic_gc_returns_ambiguous() {
        assert_eq!(
            match_alleles(('G', 'C'), 'G', 'C'),
            AlleleMatch::StrandAmbiguous
        );
    }

    #[test]
    fn non_palindromic_ag_not_ambiguous() {
        // A/G is NOT palindromic (complement of A is T, not G)
        assert_ne!(
            match_alleles(('A', 'G'), 'A', 'G'),
            AlleleMatch::StrandAmbiguous
        );
    }

    // --- match_alleles_with_frequency tests ---

    #[test]
    fn freq_resolves_palindromic_same_strand() {
        // A/T palindromic, user_af=0.30, db_af=0.28 => close => DirectMatch
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', Some(0.30), Some(0.28)),
            AlleleMatch::DirectMatch,
        );
    }

    #[test]
    fn freq_resolves_palindromic_flipped_strand() {
        // A/T palindromic, user_af=0.70, db_af=0.28
        // user_af (0.70) close to 1.0 - db_af (0.72) => ComplementMatch
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', Some(0.70), Some(0.28)),
            AlleleMatch::ComplementMatch,
        );
    }

    #[test]
    fn freq_ambiguous_when_neither_match_nor_complement() {
        // A/T palindromic, user_af=0.35, db_af=0.50 — diff=0.15 (too far for
        // direct), 1.0-db_af=0.50 so complement diff=0.15 (also too far).
        // Neither strand hypothesis is supported.
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', Some(0.35), Some(0.50)),
            AlleleMatch::StrandAmbiguous,
        );
    }

    #[test]
    fn freq_no_user_af_extreme_db_af_resolves() {
        // A/T palindromic, no user_af, db_af=0.05 (very far from 0.5)
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', None, Some(0.05)),
            AlleleMatch::DirectMatch,
        );
    }

    #[test]
    fn freq_no_user_af_moderate_db_af_stays_ambiguous() {
        // A/T palindromic, no user_af, db_af=0.30 (not extreme enough)
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', None, Some(0.30)),
            AlleleMatch::StrandAmbiguous,
        );
    }

    #[test]
    fn freq_no_data_stays_ambiguous() {
        // A/T palindromic, no frequency data at all
        assert_eq!(
            match_alleles_with_frequency(('A', 'T'), 'A', 'T', None, None),
            AlleleMatch::StrandAmbiguous,
        );
    }

    #[test]
    fn freq_non_palindromic_ignores_frequency() {
        // A/G non-palindromic — frequency irrelevant, returns DirectMatch
        assert_eq!(
            match_alleles_with_frequency(('A', 'G'), 'A', 'G', Some(0.50), Some(0.50)),
            AlleleMatch::DirectMatch,
        );
    }

    #[test]
    fn cg_palindromic_with_frequency_resolves() {
        // C/G palindromic, user_af=0.15, db_af=0.12 => close => DirectMatch
        assert_eq!(
            match_alleles_with_frequency(('C', 'G'), 'C', 'G', Some(0.15), Some(0.12)),
            AlleleMatch::DirectMatch,
        );
    }

    // --- count_risk_allele_copies tests ---

    #[test]
    fn risk_allele_direct_match_one_copy() {
        // User A/G, risk allele A -> 1 copy
        let result = count_risk_allele_copies(Some("A"), 'A', 'G');
        assert_eq!(
            result,
            RiskAlleleCopies::Determined {
                copies: 1,
                match_type: AlleleMatch::DirectMatch
            }
        );
    }

    #[test]
    fn risk_allele_direct_match_two_copies() {
        // User A/A, risk allele A -> 2 copies
        let result = count_risk_allele_copies(Some("A"), 'A', 'A');
        assert_eq!(
            result,
            RiskAlleleCopies::Determined {
                copies: 2,
                match_type: AlleleMatch::DirectMatch
            }
        );
    }

    #[test]
    fn risk_allele_complement_match() {
        // User T/G, risk allele A -> complement match (T is complement of A) -> 1 copy
        let result = count_risk_allele_copies(Some("A"), 'T', 'G');
        assert_eq!(
            result,
            RiskAlleleCopies::Determined {
                copies: 1,
                match_type: AlleleMatch::ComplementMatch
            }
        );
    }

    #[test]
    fn risk_allele_zero_copies() {
        // User G/G, risk allele A -> 0 copies
        let result = count_risk_allele_copies(Some("A"), 'G', 'G');
        assert_eq!(
            result,
            RiskAlleleCopies::Determined {
                copies: 0,
                match_type: AlleleMatch::Mismatch
            }
        );
    }

    #[test]
    fn risk_allele_palindromic_het() {
        // User A/T, risk allele A -> palindromic
        let result = count_risk_allele_copies(Some("A"), 'A', 'T');
        assert!(matches!(
            result,
            RiskAlleleCopies::Palindromic { copies: 1 }
        ));
    }

    #[test]
    fn risk_allele_none_is_indeterminate() {
        let result = count_risk_allele_copies(None, 'A', 'G');
        assert_eq!(result, RiskAlleleCopies::Indeterminate);
    }

    #[test]
    fn risk_allele_empty_is_indeterminate() {
        let result = count_risk_allele_copies(Some(""), 'A', 'G');
        assert_eq!(result, RiskAlleleCopies::Indeterminate);
    }
}
