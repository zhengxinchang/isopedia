pub const MERGED_FILE_NAME: &str = "merged_isoform.dat";
pub const CHROM_FILE_NAME: &str = "chrom.map";
pub const TMPIDX_FILE_NAME: &str = "tmp.idx";
pub const DATASET_INFO_FILE_NAME: &str = "merged_isoform.info";
pub const META_FILE_NAME: &str = "meta.txt";

pub const ORDER: u64 = 250;
pub type KeyType = u64;
pub type ValueType = u64;
pub type NodeIDType = u64;
pub type AddrType = u64;
pub const LRU_CACHE_SIZE: usize = 100000;
pub const MAGIC: u64 = 9236;
pub const MAX_SAMPLE_SIZE: usize = 256;
