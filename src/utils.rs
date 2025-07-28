use ahash::RandomState;
use std::hash::Hash;

pub fn pack_u32(high: u32, low: u32) -> u64 {
    ((high as u64) << 32) | (low as u64)
}

pub fn pack_u64(high: u64, low: u64) -> u64 {
    ((high as u64) << 32) | (low as u64)
}

pub fn u64to2u64(combined: u64) -> (u64, u64) {
    let high = (combined >> 32) as u64;
    let low = combined as u64;
    (high, low)
}

pub fn u64to2u32(combined: u64) -> (u32, u32) {
    let high = (combined >> 32) as u32;
    let low = combined as u32;
    (high, low)
}

pub fn trim_chr_prefix_to_upper(chrom: &str) -> String {
    chrom
        .to_ascii_uppercase()
        .trim_start_matches("CHR")
        .to_string()
}

pub fn pad_chrom_prefix(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

pub fn hash_vec<T: Hash>(value: &T) -> u64 {
    let hasher_builder = RandomState::with_seeds(9336, 5920, 6784, 4496);
    hasher_builder.hash_one(value)
}


pub fn is_overlap(a:&Vec<f64>, b:&Vec<f64>) -> bool {
    if a.is_empty() || b.is_empty() {
        return false;
    }
    let a_start = a[0];
    let a_end = a[a.len() - 1];
    let b_start = b[0];
    let b_end = b[b.len() - 1];

    // Check if the ranges overlap
    (a_start <= b_end) && (b_start <= a_end)
}