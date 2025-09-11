use ahash::RandomState;
use anyhow::Result;
use std::fs;
use std::{fs::File, hash::Hash, io::Read, path::Path};

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

pub fn is_overlap(a: &Vec<f64>, b: &Vec<f64>) -> bool {
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

pub fn calc_cpm(val_u32: &u32, total_u32: &u32) -> f64 {
    if *val_u32 == 0 || *total_u32 == 0 {
        return 0.0;
    }

    (*val_u32 as f64 / *total_u32 as f64) * 1_000_000.0
}

/// Warming up the archive file by reading a portion of the file.
/// Current function only read the first max_bytes, but can be improved by
/// reading the file with imputation strategies
pub fn warmup(path: &Path, max_bytes: usize) -> Result<()> {
    let mut f = File::open(path)?;
    // let file_size = get_file_size(path)?;

    let mut buf = [0u8; 8 * 1024];
    let mut total_read = 0usize;

    while total_read < max_bytes {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break; // EOF
        }
        total_read += n;
    }
    Ok(())
}

pub fn get_total_memory_bytes() -> Option<u64> {
    if let Ok(meminfo) = fs::read_to_string("/proc/meminfo") {
        for line in meminfo.lines() {
            if line.starts_with("MemTotal:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(kb) = parts[1].parse::<u64>() {
                        return Some(kb * 1024);
                    }
                }
            }
        }
    }
    None
}

#[allow(dead_code)]
fn get_file_size(path: &Path) -> Result<u64> {
    let f = File::open(path)?;
    let metadata = f.metadata()?;
    Ok(metadata.len())
}

pub fn parse_breakpoint_str(s: &str) -> Result<(String, u64, String, u64)> {
    let parts: Vec<&str> = s.split(',').collect();

    let left_part = parts[0].split(|c| c == ':').collect::<Vec<&str>>();
    let left_chr = trim_chr_prefix_to_upper(left_part[0]);
    let left_pos: u64 = left_part[1].parse()?;

    let right_part = parts[1].split(|c| c == ':').collect::<Vec<&str>>();
    let right_chr = trim_chr_prefix_to_upper(right_part[0]);
    let right_pos: u64 = right_part[1].parse()?;
    Ok((left_chr, left_pos, right_chr, right_pos))
}

pub fn u64diff2i32(a: u64, b: u64) -> i32 {
    if a >= b {
        let diff = a - b;
        if diff > i32::MAX as u64 {
            i32::MAX
        } else {
            diff as i32
        }
    } else {
        let diff = b - a;
        if diff > i32::MAX as u64 {
            -i32::MAX
        } else {
            -(diff as i32)
        }
    }
}
