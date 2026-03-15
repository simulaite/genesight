//! Low-level strand utilities for DNA allele handling.
//!
//! Provides complement mapping and palindromic SNP detection. These are the
//! building blocks used by the allele matching logic in the parent module.

/// Return the Watson-Crick complement of a DNA base.
///
/// Maps A<->T, C<->G (case-insensitive). Returns `None` for
/// non-standard bases (N, indels, etc.).
///
/// # Examples
///
/// ```
/// use genesight_core::allele::strand::complement;
/// assert_eq!(complement('A'), Some('T'));
/// assert_eq!(complement('C'), Some('G'));
/// assert_eq!(complement('N'), None);
/// ```
pub fn complement(base: char) -> Option<char> {
    match base.to_ascii_uppercase() {
        'A' => Some('T'),
        'T' => Some('A'),
        'C' => Some('G'),
        'G' => Some('C'),
        _ => None,
    }
}

/// Check whether a pair of alleles forms a palindromic SNP.
///
/// Palindromic SNPs have alleles that are complements of each other (A/T or
/// C/G). These are problematic because the genotyping array cannot
/// distinguish the forward strand from the reverse strand based on sequence
/// alone -- both strands read the same pair of bases.
///
/// # Examples
///
/// ```
/// use genesight_core::allele::strand::is_palindromic;
/// assert!(is_palindromic('A', 'T'));
/// assert!(is_palindromic('C', 'G'));
/// assert!(!is_palindromic('A', 'G'));
/// ```
pub fn is_palindromic(allele1: char, allele2: char) -> bool {
    let a = allele1.to_ascii_uppercase();
    let b = allele2.to_ascii_uppercase();
    matches!((a, b), ('A', 'T') | ('T', 'A') | ('C', 'G') | ('G', 'C'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_standard_bases() {
        assert_eq!(complement('A'), Some('T'));
        assert_eq!(complement('T'), Some('A'));
        assert_eq!(complement('C'), Some('G'));
        assert_eq!(complement('G'), Some('C'));
    }

    #[test]
    fn complement_case_insensitive() {
        assert_eq!(complement('a'), Some('T'));
        assert_eq!(complement('t'), Some('A'));
        assert_eq!(complement('c'), Some('G'));
        assert_eq!(complement('g'), Some('C'));
    }

    #[test]
    fn complement_non_standard_returns_none() {
        assert_eq!(complement('N'), None);
        assert_eq!(complement('X'), None);
        assert_eq!(complement('-'), None);
    }

    #[test]
    fn palindromic_at_pairs() {
        assert!(is_palindromic('A', 'T'));
        assert!(is_palindromic('T', 'A'));
    }

    #[test]
    fn palindromic_cg_pairs() {
        assert!(is_palindromic('C', 'G'));
        assert!(is_palindromic('G', 'C'));
    }

    #[test]
    fn non_palindromic_pairs() {
        assert!(!is_palindromic('A', 'G'));
        assert!(!is_palindromic('A', 'C'));
        assert!(!is_palindromic('T', 'G'));
        assert!(!is_palindromic('T', 'C'));
        assert!(!is_palindromic('G', 'A'));
        assert!(!is_palindromic('C', 'A'));
    }

    #[test]
    fn palindromic_case_insensitive() {
        assert!(is_palindromic('a', 'T'));
        assert!(is_palindromic('c', 'g'));
    }
}
