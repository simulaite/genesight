use serde::{Deserialize, Serialize};

/// Configuration controlling which databases are queried during annotation.
///
/// Each field corresponds to a database adapter. Setting a field to `false`
/// skips that database entirely, producing no annotations of that type.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationConfig {
    /// Query ClinVar for clinical variant classifications.
    pub clinvar: bool,
    /// Query GWAS Catalog for genome-wide association study hits.
    pub gwas: bool,
    /// Query dbSNP/gnomAD for allele frequency data.
    pub frequencies: bool,
    /// Query PharmGKB for pharmacogenomic annotations.
    pub pharmacogenomics: bool,
}

impl Default for AnnotationConfig {
    fn default() -> Self {
        Self {
            clinvar: true,
            gwas: true,
            frequencies: true,
            pharmacogenomics: true,
        }
    }
}
