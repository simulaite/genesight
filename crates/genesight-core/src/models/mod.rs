pub mod annotation;
pub mod confidence;
pub mod report;
pub mod variant;

pub use annotation::AnnotatedVariant;
pub use confidence::ConfidenceTier;
pub use report::{Report, ResultCategory, ScoredResult};
pub use variant::{Genotype, SourceFormat, Variant};
