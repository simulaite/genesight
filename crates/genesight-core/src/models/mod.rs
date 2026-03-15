pub mod annotation;
pub mod assembly;
pub mod confidence;
pub mod config;
pub mod report;
pub mod variant;

pub use annotation::AnnotatedVariant;
pub use assembly::GenomeAssembly;
pub use confidence::ConfidenceTier;
pub use config::AnnotationConfig;
pub use report::{ConfirmationUrgency, Report, ResultCategory, ScoredResult};
pub use variant::{Genotype, SourceFormat, Variant};
