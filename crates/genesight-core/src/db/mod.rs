//! Database adapters for querying local SQLite copies of genome databases.

pub mod clinvar;
pub mod dbsnp;
pub mod gwas;
pub mod pharmgkb;
pub mod snpedia;

use crate::models::GenomeAssembly;
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

/// Query the genome assembly version stored in the database metadata.
///
/// Reads the `assembly` key from the `db_metadata` table. Returns
/// `GenomeAssembly::Unknown` if the table does not exist, the key is
/// missing, or the value is not a recognized assembly identifier.
pub fn query_db_assembly(conn: &Connection) -> GenomeAssembly {
    let result: Result<String, _> = conn.query_row(
        "SELECT value FROM db_metadata WHERE key = 'assembly'",
        [],
        |row| row.get(0),
    );

    match result {
        Ok(value) => GenomeAssembly::from_header_line(&value),
        Err(_) => GenomeAssembly::Unknown,
    }
}
