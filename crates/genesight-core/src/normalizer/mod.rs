//! VCF variant normalization: multi-allelic splitting and allele trimming.
//!
//! When a VCF record has multiple ALT alleles (e.g., `ALT=G,T`), it must be
//! split into separate biallelic records, each with a correctly remapped GT
//! field. After splitting, alleles are trimmed to their minimal representation
//! (right-trim shared suffix, then left-trim shared prefix while incrementing
//! position).
//!
//! This module is used internally by the VCF parser and is not intended for
//! direct use by downstream consumers.

use serde::{Deserialize, Serialize};

/// Describes what normalization was applied to produce a record.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NormalizationStatus {
    /// No normalization was needed; the record was already canonical.
    Original,
    /// Alleles were trimmed (shared prefix/suffix removed).
    Trimmed {
        /// Number of leading bases removed (position was incremented by this amount).
        leading: u8,
        /// Number of trailing bases removed.
        trailing: u8,
    },
    /// Record was produced by splitting a multi-allelic site.
    MultiAllelicSplit {
        /// Zero-based index of the ALT allele this record represents.
        alt_index: usize,
        /// Total number of ALT alleles in the original record.
        total_alts: usize,
    },
    /// Normalization could not be applied (e.g., symbolic alleles like `<DEL>`).
    NormalizationFailed(String),
}

/// A single normalized VCF record ready for conversion to a `Variant`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedRecord {
    /// Reference allele after normalization.
    pub ref_allele: String,
    /// Alternate allele after normalization.
    pub alt_allele: String,
    /// Genomic position after normalization (may differ from original if trimmed).
    pub position: u64,
    /// Genotype string remapped for this specific split (e.g., "0/1").
    pub gt: String,
    /// What normalization was applied.
    pub status: NormalizationStatus,
}

/// Normalize a VCF record, splitting multi-allelic sites and trimming alleles.
///
/// For a record with N ALT alleles, returns N `NormalizedRecord`s (one per ALT).
/// A simple biallelic record returns exactly one record.
///
/// # Arguments
///
/// * `pos` - 1-based genomic position from the VCF POS column
/// * `ref_allele` - Reference allele string from the VCF REF column
/// * `alt_alleles` - ALT alleles (the ALT column split on `,`)
/// * `gt` - Genotype string from the GT field (e.g., "0/1", "1/2", "1|0")
///
/// # Returns
///
/// A vector of normalized records, one per ALT allele.
pub fn normalize_vcf_record(
    pos: u64,
    ref_allele: &str,
    alt_alleles: &[&str],
    gt: &str,
) -> Vec<NormalizedRecord> {
    let total_alts = alt_alleles.len();

    // Single ALT allele: no splitting needed
    if total_alts == 1 {
        let alt = alt_alleles[0];

        // Check for symbolic alleles
        if is_symbolic_allele(alt) {
            return vec![NormalizedRecord {
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                position: pos,
                gt: gt.to_string(),
                status: NormalizationStatus::NormalizationFailed(format!("symbolic allele: {alt}")),
            }];
        }

        // Try trimming
        let (trimmed_ref, trimmed_alt, trimmed_pos, leading, trailing) =
            trim_alleles(ref_allele, alt, pos);

        let status = if leading > 0 || trailing > 0 {
            NormalizationStatus::Trimmed {
                leading: leading as u8,
                trailing: trailing as u8,
            }
        } else {
            NormalizationStatus::Original
        };

        return vec![NormalizedRecord {
            ref_allele: trimmed_ref,
            alt_allele: trimmed_alt,
            position: trimmed_pos,
            gt: gt.to_string(),
            status,
        }];
    }

    // Multi-allelic: split into one record per ALT
    let mut records = Vec::with_capacity(total_alts);

    for (alt_idx, alt) in alt_alleles.iter().enumerate() {
        // Check for symbolic alleles
        if is_symbolic_allele(alt) {
            records.push(NormalizedRecord {
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                position: pos,
                gt: resolve_gt_for_split(gt, alt_idx, total_alts),
                status: NormalizationStatus::NormalizationFailed(format!("symbolic allele: {alt}")),
            });
            continue;
        }

        let remapped_gt = resolve_gt_for_split(gt, alt_idx, total_alts);

        // Trim alleles for this split
        let (trimmed_ref, trimmed_alt, trimmed_pos, _leading, _trailing) =
            trim_alleles(ref_allele, alt, pos);

        records.push(NormalizedRecord {
            ref_allele: trimmed_ref,
            alt_allele: trimmed_alt,
            position: trimmed_pos,
            gt: remapped_gt,
            status: NormalizationStatus::MultiAllelicSplit {
                alt_index: alt_idx,
                total_alts,
            },
        });
    }

    records
}

/// Check whether an allele is a symbolic allele (e.g., `<DEL>`, `<INS>`, `<DUP>`).
fn is_symbolic_allele(allele: &str) -> bool {
    allele.starts_with('<') && allele.ends_with('>')
}

/// Trim shared prefix and suffix from REF and ALT alleles.
///
/// VCF records often include flanking context bases. This function removes:
/// 1. First: shared suffix (right-trim), preserving at least 1 base on each side
/// 2. Then: shared prefix (left-trim), incrementing `pos` accordingly
///
/// Returns `(trimmed_ref, trimmed_alt, new_pos, leading_trimmed, trailing_trimmed)`.
pub fn trim_alleles(
    ref_allele: &str,
    alt_allele: &str,
    pos: u64,
) -> (String, String, u64, usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    // Both alleles must have at least 1 base
    if ref_bytes.is_empty() || alt_bytes.is_empty() {
        return (ref_allele.to_string(), alt_allele.to_string(), pos, 0, 0);
    }

    // Step 1: Right-trim shared suffix (but keep at least 1 base on each side)
    let max_suffix = ref_bytes.len().min(alt_bytes.len()) - 1; // -1 to keep at least 1 base
    let mut suffix_len = 0;
    for i in 0..max_suffix {
        if ref_bytes[ref_bytes.len() - 1 - i] == alt_bytes[alt_bytes.len() - 1 - i] {
            suffix_len += 1;
        } else {
            break;
        }
    }

    let ref_after_suffix = &ref_bytes[..ref_bytes.len() - suffix_len];
    let alt_after_suffix = &alt_bytes[..alt_bytes.len() - suffix_len];

    // Step 2: Left-trim shared prefix (but keep at least 1 base on each side)
    let max_prefix = ref_after_suffix.len().min(alt_after_suffix.len()) - 1;
    let mut prefix_len = 0;
    for i in 0..max_prefix {
        if ref_after_suffix[i] == alt_after_suffix[i] {
            prefix_len += 1;
        } else {
            break;
        }
    }

    let trimmed_ref = &ref_after_suffix[prefix_len..];
    let trimmed_alt = &alt_after_suffix[prefix_len..];

    // Safety: VCF alleles are ASCII, so this conversion is safe
    let trimmed_ref_str = std::str::from_utf8(trimmed_ref).unwrap_or(ref_allele);
    let trimmed_alt_str = std::str::from_utf8(trimmed_alt).unwrap_or(alt_allele);

    (
        trimmed_ref_str.to_string(),
        trimmed_alt_str.to_string(),
        pos + prefix_len as u64,
        prefix_len,
        suffix_len,
    )
}

/// Remap a genotype string for a specific ALT allele after multi-allelic splitting.
///
/// In a multi-allelic VCF record, the GT field uses indices where 0 = REF,
/// 1 = first ALT, 2 = second ALT, etc. When splitting into biallelic records,
/// each split record needs its GT adjusted so the target ALT becomes index 1
/// and all other ALTs become `.` (missing).
///
/// # Arguments
///
/// * `gt` - Original GT string (e.g., "0/1", "1/2", "1|0")
/// * `target_alt_index` - Zero-based index of the ALT allele for this split
/// * `total_alts` - Total number of ALT alleles in the original record
///
/// # Examples
///
/// For `ALT=G,T` (total_alts=2):
/// - GT `"0/1"`, target=0 -> `"0/1"` (this split carries the first ALT)
/// - GT `"0/1"`, target=1 -> `"0/."` (second ALT not in this GT)
/// - GT `"1/2"`, target=0 -> `"1/."` (first ALT present, remap 1->1)
/// - GT `"1/2"`, target=1 -> `"./1"` (second ALT present, remap 2->1)
/// - GT `"1/1"`, target=0 -> `"1/1"` (homozygous for first ALT)
/// - GT `"1/1"`, target=1 -> `"./."` (not homozygous for second ALT)
pub fn resolve_gt_for_split(gt: &str, target_alt_index: usize, _total_alts: usize) -> String {
    // The 1-based ALT index in the original GT that we are targeting
    let target_gt_index = target_alt_index + 1;

    // Determine the separator
    let (separator, alleles_str): (char, Vec<&str>) = if gt.contains('|') {
        ('|', gt.split('|').collect())
    } else if gt.contains('/') {
        ('/', gt.split('/').collect())
    } else {
        // Hemizygous (single allele, no separator)
        return remap_single_allele(gt, target_gt_index);
    };

    let remapped: Vec<String> = alleles_str
        .iter()
        .map(|&a| remap_allele_index(a, target_gt_index))
        .collect();

    remapped.join(&separator.to_string())
}

/// Remap a single allele index for a split record.
///
/// - `"0"` (REF) stays `"0"`
/// - `"."` (missing) stays `"."`
/// - index matching `target_gt_index` becomes `"1"` (the ALT in this split)
/// - any other index becomes `"."` (not relevant to this split)
fn remap_allele_index(allele: &str, target_gt_index: usize) -> String {
    match allele {
        "0" => "0".to_string(),
        "." => ".".to_string(),
        other => match other.parse::<usize>() {
            Ok(idx) if idx == target_gt_index => "1".to_string(),
            Ok(_) => ".".to_string(),
            Err(_) => ".".to_string(),
        },
    }
}

/// Handle hemizygous (single allele, no separator) GT for split records.
fn remap_single_allele(allele: &str, target_gt_index: usize) -> String {
    match allele {
        "0" => "0".to_string(),
        "." => ".".to_string(),
        other => match other.parse::<usize>() {
            Ok(idx) if idx == target_gt_index => "1".to_string(),
            Ok(_) => ".".to_string(),
            Err(_) => ".".to_string(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Test 1: Simple biallelic SNP -> unchanged, NormalizationStatus::Original
    // -----------------------------------------------------------------------
    #[test]
    fn biallelic_snp_unchanged() {
        let records = normalize_vcf_record(82154, "G", &["A"], "0/1");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ref_allele, "G");
        assert_eq!(records[0].alt_allele, "A");
        assert_eq!(records[0].position, 82154);
        assert_eq!(records[0].gt, "0/1");
        assert_eq!(records[0].status, NormalizationStatus::Original);
    }

    // -----------------------------------------------------------------------
    // Test 2: Multi-allelic ALT=G,T with GT 0/1 -> 2 records
    // -----------------------------------------------------------------------
    #[test]
    fn multiallelic_gt_0_1() {
        let records = normalize_vcf_record(100, "A", &["G", "T"], "0/1");
        assert_eq!(records.len(), 2);

        // Split 0: ALT=G, GT should be 0/1 (carries first ALT)
        assert_eq!(records[0].alt_allele, "G");
        assert_eq!(records[0].gt, "0/1");
        assert_eq!(
            records[0].status,
            NormalizationStatus::MultiAllelicSplit {
                alt_index: 0,
                total_alts: 2,
            }
        );

        // Split 1: ALT=T, GT should be 0/. (second ALT not in GT)
        assert_eq!(records[1].alt_allele, "T");
        assert_eq!(records[1].gt, "0/.");
        assert_eq!(
            records[1].status,
            NormalizationStatus::MultiAllelicSplit {
                alt_index: 1,
                total_alts: 2,
            }
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: Multi-allelic with GT 1/2 -> 2 records, each het for their ALT
    // -----------------------------------------------------------------------
    #[test]
    fn multiallelic_gt_1_2() {
        let records = normalize_vcf_record(100, "A", &["G", "T"], "1/2");
        assert_eq!(records.len(), 2);

        // Split 0: ALT=G (index 1 in original), GT should be 1/.
        assert_eq!(records[0].alt_allele, "G");
        assert_eq!(records[0].gt, "1/.");

        // Split 1: ALT=T (index 2 in original), GT should be ./1
        assert_eq!(records[1].alt_allele, "T");
        assert_eq!(records[1].gt, "./1");
    }

    // -----------------------------------------------------------------------
    // Test 4: Multi-allelic with GT 1/1 -> split 0 hom, split 1 nocall
    // -----------------------------------------------------------------------
    #[test]
    fn multiallelic_gt_1_1() {
        let records = normalize_vcf_record(100, "A", &["G", "T"], "1/1");
        assert_eq!(records.len(), 2);

        // Split 0: hom for first ALT
        assert_eq!(records[0].alt_allele, "G");
        assert_eq!(records[0].gt, "1/1");

        // Split 1: no call (not hom for second ALT)
        assert_eq!(records[1].alt_allele, "T");
        assert_eq!(records[1].gt, "./.");
    }

    // -----------------------------------------------------------------------
    // Test 5: Trimming REF=ATG ALT=ACG at pos 100 -> REF=T ALT=C at pos 101
    // -----------------------------------------------------------------------
    #[test]
    fn trim_shared_prefix_suffix() {
        let records = normalize_vcf_record(100, "ATG", &["ACG"], "0/1");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ref_allele, "T");
        assert_eq!(records[0].alt_allele, "C");
        assert_eq!(records[0].position, 101);
        assert!(matches!(
            records[0].status,
            NormalizationStatus::Trimmed {
                leading: 1,
                trailing: 1
            }
        ));
    }

    // -----------------------------------------------------------------------
    // Test 6: Deletion REF=AT ALT=A -> stays REF=AT ALT=A (already minimal
    //         with the anchor base)
    // -----------------------------------------------------------------------
    #[test]
    fn deletion_already_canonical() {
        let records = normalize_vcf_record(100, "AT", &["A"], "0/1");
        assert_eq!(records.len(), 1);
        // After trimming: right-trim nothing (no shared suffix beyond anchor),
        // left-trim nothing (A vs A is the anchor base, but we keep at least 1)
        // Actually: REF=AT, ALT=A. Shared suffix: none (T != A). Shared prefix: A.
        // But we must keep 1 base, so prefix trim = 0 for ALT (len 1 - 1 = 0 max prefix).
        // The result is unchanged.
        assert_eq!(records[0].ref_allele, "AT");
        assert_eq!(records[0].alt_allele, "A");
        assert_eq!(records[0].position, 100);
        assert_eq!(records[0].status, NormalizationStatus::Original);
    }

    // -----------------------------------------------------------------------
    // Test 7: Symbolic allele <DEL> -> NormalizationFailed
    // -----------------------------------------------------------------------
    #[test]
    fn symbolic_allele_fails() {
        let records = normalize_vcf_record(100, "A", &["<DEL>"], "0/1");
        assert_eq!(records.len(), 1);
        assert!(matches!(
            records[0].status,
            NormalizationStatus::NormalizationFailed(_)
        ));
        assert_eq!(records[0].gt, "0/1");
    }

    // -----------------------------------------------------------------------
    // Additional: phased GT
    // -----------------------------------------------------------------------
    #[test]
    fn phased_gt_preserved() {
        let records = normalize_vcf_record(100, "A", &["G", "T"], "0|2");
        assert_eq!(records.len(), 2);

        // Split 0: ALT=G, original GT has 0 and 2. Index 1 (G) not present -> 0|.
        assert_eq!(records[0].gt, "0|.");

        // Split 1: ALT=T, original GT has 0 and 2. Index 2 (T) -> 0|1
        assert_eq!(records[1].gt, "0|1");
    }

    // -----------------------------------------------------------------------
    // Additional: hemizygous GT
    // -----------------------------------------------------------------------
    #[test]
    fn hemizygous_gt_split() {
        let records = normalize_vcf_record(100, "A", &["G", "T"], "2");
        assert_eq!(records.len(), 2);

        // Split 0: target is index 1 (G), GT "2" doesn't match -> "."
        assert_eq!(records[0].gt, ".");

        // Split 1: target is index 2 (T), GT "2" matches -> "1"
        assert_eq!(records[1].gt, "1");
    }

    // -----------------------------------------------------------------------
    // Additional: resolve_gt_for_split unit tests
    // -----------------------------------------------------------------------
    #[test]
    fn resolve_gt_ref_hom() {
        assert_eq!(resolve_gt_for_split("0/0", 0, 2), "0/0");
        assert_eq!(resolve_gt_for_split("0/0", 1, 2), "0/0");
    }

    #[test]
    fn resolve_gt_missing() {
        assert_eq!(resolve_gt_for_split("./.", 0, 2), "./.");
    }

    // -----------------------------------------------------------------------
    // trim_alleles unit tests
    // -----------------------------------------------------------------------
    #[test]
    fn trim_simple_snp() {
        let (r, a, p, lead, trail) = trim_alleles("A", "G", 100);
        assert_eq!(r, "A");
        assert_eq!(a, "G");
        assert_eq!(p, 100);
        assert_eq!(lead, 0);
        assert_eq!(trail, 0);
    }

    #[test]
    fn trim_insertion() {
        // REF=A, ALT=AGT -> no shared suffix, no shared prefix beyond anchor
        let (r, a, p, lead, trail) = trim_alleles("A", "AGT", 100);
        assert_eq!(r, "A");
        assert_eq!(a, "AGT");
        assert_eq!(p, 100);
        assert_eq!(lead, 0);
        assert_eq!(trail, 0);
    }

    #[test]
    fn trim_complex_prefix_suffix() {
        // REF=ATCG, ALT=ACCG -> shared suffix G, shared prefix A -> REF=TC, ALT=CC after suffix trim,
        // then prefix trim: T vs C, no shared prefix -> REF=TC, ALT=CC at pos+0
        // Wait, let me recalculate:
        // REF=ATCG, ALT=ACCG
        // Step 1: right-trim. max_suffix = min(4,4)-1 = 3
        //   i=0: G == G -> suffix_len=1
        //   i=1: C == C -> suffix_len=2
        //   i=2: T != C -> stop
        // After suffix trim: REF=AT, ALT=AC
        // Step 2: left-trim. max_prefix = min(2,2)-1 = 1
        //   i=0: A == A -> prefix_len=1
        // After prefix trim: REF=T, ALT=C at pos+1
        let (r, a, p, lead, trail) = trim_alleles("ATCG", "ACCG", 100);
        assert_eq!(r, "T");
        assert_eq!(a, "C");
        assert_eq!(p, 101);
        assert_eq!(lead, 1);
        assert_eq!(trail, 2);
    }

    // -----------------------------------------------------------------------
    // Three-allelic split
    // -----------------------------------------------------------------------
    #[test]
    fn three_allelic_split() {
        let records = normalize_vcf_record(100, "A", &["G", "T", "C"], "1/3");
        assert_eq!(records.len(), 3);

        // Split 0: ALT=G (index 1 matches), GT 1/3 -> 1/.
        assert_eq!(records[0].gt, "1/.");

        // Split 1: ALT=T (index 2), GT 1/3 -> ./.
        assert_eq!(records[1].gt, "./.");

        // Split 2: ALT=C (index 3 matches), GT 1/3 -> ./1
        assert_eq!(records[2].gt, "./1");
    }
}
