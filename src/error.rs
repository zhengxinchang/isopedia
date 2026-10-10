use thiserror::Error;
use std::io::Error as ioError;


#[derive(Debug, Error)]
pub enum ProfileError {
    #[error("invalid argument: {0}")]
    InvalidArgument(String),

    #[error("{context}: {source}")]
    Io {
        context: String,
        #[source]
        source: std::io::Error,
    },

    #[error("{context}: {source}")]
    Hts {
        context: String,
        #[source]
        source: rust_htslib::errors::Error,
    },

    #[error("invalid record {read_name:?}: {reason}")]
    InvalidRecord {
        read_name: String,
        reason: String,
    },
}