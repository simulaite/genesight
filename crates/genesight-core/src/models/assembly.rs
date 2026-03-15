use serde::{Deserialize, Serialize};

/// Human genome reference assembly version.
///
/// DNA raw data files and databases are built against a specific genome
/// assembly. Mismatches between the input file assembly and the database
/// assembly can cause incorrect position-based lookups. This enum tracks
/// which assembly is in use so that warnings can be generated when there
/// is a potential mismatch.
///
/// Phase 1 is informational only — no coordinate LiftOver is performed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GenomeAssembly {
    /// GRCh37 / hg19 (used by most 23andMe and AncestryDNA files)
    GRCh37,
    /// GRCh38 / hg38 (newer assembly, used by some VCF files)
    GRCh38,
    /// Assembly could not be determined from the file or database
    Unknown,
}

impl std::fmt::Display for GenomeAssembly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GenomeAssembly::GRCh37 => write!(f, "GRCh37 (hg19)"),
            GenomeAssembly::GRCh38 => write!(f, "GRCh38 (hg38)"),
            GenomeAssembly::Unknown => write!(f, "Unknown"),
        }
    }
}

impl GenomeAssembly {
    /// Attempt to detect the genome assembly from a single header line.
    ///
    /// Performs a case-insensitive scan for common assembly identifiers:
    /// - GRCh37: "build 37", "hg19", "grch37", "b37"
    /// - GRCh38: "build 38", "hg38", "grch38", "b38"
    ///
    /// Returns `Unknown` if no recognizable identifier is found.
    pub fn from_header_line(s: &str) -> GenomeAssembly {
        let lower = s.to_lowercase();

        // Check GRCh37 identifiers
        if lower.contains("grch37")
            || lower.contains("hg19")
            || lower.contains("build 37")
            || lower.contains("b37")
        {
            return GenomeAssembly::GRCh37;
        }

        // Check GRCh38 identifiers
        if lower.contains("grch38")
            || lower.contains("hg38")
            || lower.contains("build 38")
            || lower.contains("b38")
        {
            return GenomeAssembly::GRCh38;
        }

        GenomeAssembly::Unknown
    }

    /// Check whether two assemblies are compatible.
    ///
    /// `Unknown` is considered compatible with any assembly (since we cannot
    /// determine a mismatch). Two known assemblies are compatible only if
    /// they are equal.
    pub fn is_compatible_with(self, other: GenomeAssembly) -> bool {
        match (self, other) {
            (GenomeAssembly::Unknown, _) | (_, GenomeAssembly::Unknown) => true,
            (a, b) => a == b,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_grch37() {
        assert_eq!(GenomeAssembly::GRCh37.to_string(), "GRCh37 (hg19)");
    }

    #[test]
    fn display_grch38() {
        assert_eq!(GenomeAssembly::GRCh38.to_string(), "GRCh38 (hg38)");
    }

    #[test]
    fn display_unknown() {
        assert_eq!(GenomeAssembly::Unknown.to_string(), "Unknown");
    }

    #[test]
    fn from_header_line_grch37_variants() {
        assert_eq!(
            GenomeAssembly::from_header_line("# build 37"),
            GenomeAssembly::GRCh37
        );
        assert_eq!(
            GenomeAssembly::from_header_line("##reference=hg19"),
            GenomeAssembly::GRCh37
        );
        assert_eq!(
            GenomeAssembly::from_header_line("# This data is on GRCh37"),
            GenomeAssembly::GRCh37
        );
        assert_eq!(
            GenomeAssembly::from_header_line("##contig=<ID=1,assembly=b37>"),
            GenomeAssembly::GRCh37
        );
    }

    #[test]
    fn from_header_line_grch38_variants() {
        assert_eq!(
            GenomeAssembly::from_header_line("# build 38"),
            GenomeAssembly::GRCh38
        );
        assert_eq!(
            GenomeAssembly::from_header_line("##reference=hg38"),
            GenomeAssembly::GRCh38
        );
        assert_eq!(
            GenomeAssembly::from_header_line("##reference=GRCh38"),
            GenomeAssembly::GRCh38
        );
        assert_eq!(
            GenomeAssembly::from_header_line("##contig=<ID=1,assembly=b38>"),
            GenomeAssembly::GRCh38
        );
    }

    #[test]
    fn from_header_line_case_insensitive() {
        assert_eq!(
            GenomeAssembly::from_header_line("# BUILD 37"),
            GenomeAssembly::GRCh37
        );
        assert_eq!(
            GenomeAssembly::from_header_line("# GRCH38 reference"),
            GenomeAssembly::GRCh38
        );
    }

    #[test]
    fn from_header_line_unknown() {
        assert_eq!(
            GenomeAssembly::from_header_line("# some random comment"),
            GenomeAssembly::Unknown
        );
        assert_eq!(
            GenomeAssembly::from_header_line("rsid\tchromosome\tposition\tgenotype"),
            GenomeAssembly::Unknown
        );
    }

    #[test]
    fn compatibility_same_assembly() {
        assert!(GenomeAssembly::GRCh37.is_compatible_with(GenomeAssembly::GRCh37));
        assert!(GenomeAssembly::GRCh38.is_compatible_with(GenomeAssembly::GRCh38));
    }

    #[test]
    fn compatibility_different_assemblies() {
        assert!(!GenomeAssembly::GRCh37.is_compatible_with(GenomeAssembly::GRCh38));
        assert!(!GenomeAssembly::GRCh38.is_compatible_with(GenomeAssembly::GRCh37));
    }

    #[test]
    fn compatibility_unknown_is_always_compatible() {
        assert!(GenomeAssembly::Unknown.is_compatible_with(GenomeAssembly::GRCh37));
        assert!(GenomeAssembly::Unknown.is_compatible_with(GenomeAssembly::GRCh38));
        assert!(GenomeAssembly::Unknown.is_compatible_with(GenomeAssembly::Unknown));
        assert!(GenomeAssembly::GRCh37.is_compatible_with(GenomeAssembly::Unknown));
        assert!(GenomeAssembly::GRCh38.is_compatible_with(GenomeAssembly::Unknown));
    }
}
