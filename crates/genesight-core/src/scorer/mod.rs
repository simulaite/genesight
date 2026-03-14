//! Risk scoring and confidence tier assignment.
//!
//! Every finding MUST be assigned a [`ConfidenceTier`] based on the
//! quality and type of evidence.

pub mod monogenic;
pub mod pharmaco;
pub mod polygenic;
pub mod traits;
