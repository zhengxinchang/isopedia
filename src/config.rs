/// Built-in aliases, applied after uppercasing and removing the CHR prefix.
pub(crate) const CHROMOSOME_ALIASES: &[(&str, &str)] = &[("M", "MT")];

pub(crate) const PROFILE_SCHEMA_VERSION: u32 = 6;

/// PTIR layout without a persisted signature; signature reconstruction uses schema 6.
pub(crate) const ARCHIVE_SCHEMA_VERSION: u32 = 1;
