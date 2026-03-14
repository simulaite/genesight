//! Variant annotation engine.
//!
//! Matches parsed variants against local database entries to produce
//! annotated findings with confidence tiers.

pub mod clinical;
pub mod frequency;
pub mod pharmacogenomics;
pub mod traits;
