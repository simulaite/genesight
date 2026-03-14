use serde::{Deserialize, Serialize};

/// A single genetic variant parsed from a DNA raw data file.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Variant {
    /// The rsID identifier (e.g., "rs4477212"), None for novel variants
    pub rsid: Option<String>,
    /// Chromosome (1-22, X, Y, MT)
    pub chromosome: String,
    /// Genomic position
    pub position: u64,
    /// The genotype call
    pub genotype: Genotype,
    /// Which file format this variant was parsed from
    pub source_format: SourceFormat,
}

/// Represents a genotype call for a variant.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Genotype {
    /// Homozygous genotype (e.g., AA, GG)
    Homozygous(char),
    /// Heterozygous genotype (e.g., AG, CT)
    Heterozygous(char, char),
    /// No call (variant could not be determined)
    NoCall,
    /// Insertion/Deletion
    Indel(String),
}

/// The DNA file format a variant was parsed from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SourceFormat {
    TwentyThreeAndMe,
    AncestryDNA,
    Vcf,
}

impl std::fmt::Display for Genotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Genotype::Homozygous(a) => write!(f, "{a}{a}"),
            Genotype::Heterozygous(a, b) => write!(f, "{a}{b}"),
            Genotype::NoCall => write!(f, "--"),
            Genotype::Indel(s) => write!(f, "{s}"),
        }
    }
}
