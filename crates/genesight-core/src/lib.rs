//! GeneSight Core Library
//!
//! Privacy-first DNA analysis engine that annotates personal genetic data
//! against local copies of public genome databases.
//!
//! # Architecture
//!
//! - **parser** — Read DNA raw data files (23andMe, AncestryDNA, VCF)
//! - **db** — Query local SQLite databases (ClinVar, SNPedia, GWAS, dbSNP, PharmGKB)
//! - **annotator** — Match variants against database entries
//! - **scorer** — Assign confidence tiers and risk scores
//! - **report** — Generate human-readable reports (Markdown, JSON, HTML)
//! - **models** — Shared types and data structures

pub mod annotator;
pub mod db;
pub mod models;
pub mod parser;
pub mod report;
pub mod scorer;
