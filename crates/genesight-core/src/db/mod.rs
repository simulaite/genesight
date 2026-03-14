//! Database adapters for querying local SQLite copies of genome databases.

pub mod clinvar;
pub mod dbsnp;
pub mod gwas;
pub mod pharmgkb;
pub mod snpedia;

use rusqlite::Connection;

/// Errors that can occur during database operations.
#[derive(Debug, thiserror::Error)]
pub enum DbError {
    #[error("SQLite error: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("Database not found at {path}")]
    NotFound { path: String },
    #[error("Database schema version mismatch: expected {expected}, found {found}")]
    SchemaMismatch { expected: u32, found: u32 },
}

/// Open the GeneSight database at the given path.
pub fn open_database(path: &std::path::Path) -> Result<Connection, DbError> {
    if !path.exists() {
        return Err(DbError::NotFound {
            path: path.display().to_string(),
        });
    }
    let conn = Connection::open(path)?;
    Ok(conn)
}
