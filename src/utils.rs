use ahash::RandomState;
use anyhow::Result;
use log::error;
use std::fs;
use std::path::PathBuf;
use std::{fs::File, hash::Hash, io::Read, path::Path};

use crate::constants::*;

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

pub fn calc_cpm_f32(val_f32: &f32, total_f32: &f32) -> f64 {
    if *val_f32 == 0.0 || *total_f32 == 0.0 {
        return 0.0;
    }

    (*val_f32 as f64 / *total_f32 as f64) * 1_000_000.0
}

/// Warming up the archive file by reading a portion of the file.
/// Current function only read the first max_bytes, but can be improved by
/// reading the file with imputation strategies
pub fn warmup(path: &Path, max_bytes: usize) -> Result<()> {
    if max_bytes == 0 {
        return Ok(());
    }

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
    let left_chr = trim_chr_prefix_to_upper(left_part[0].trim());
    let left_pos: u64 = left_part[1].parse()?;

    let right_part = parts[1].split(|c| c == ':').collect::<Vec<&str>>();
    let right_chr = trim_chr_prefix_to_upper(right_part[0].trim());
    let right_pos: u64 = right_part[1].parse()?;
    Ok((left_chr, left_pos, right_chr, right_pos))
}

pub fn parse_splice_junction_str(s: &str) -> Result<(String, u64, u64)> {
    let parts: Vec<&str> = s.split(|c| c == ':' || c == ',').collect();

    let chr = trim_chr_prefix_to_upper(parts[0]);
    let left_pos: u64 = parts[1].parse()?;
    let right_pos: u64 = parts[2].parse()?;

    Ok((chr, left_pos, right_pos))
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

pub fn check_index_dir(path: &Path) -> bool {
    if !path.is_dir() {
        error!("Index directory {:?} is not a directory", path);
        return false;
    }

    if !path.exists() {
        error!("Index directory {:?} does not exist", path);
        return false;
    }

    if !path.join(MERGED_FILE_NAME).exists() {
        error!(
            "Index directory {:?} does not contain {}",
            path, MERGED_FILE_NAME
        );
        return false;
    }

    if !path.join("bptree_0.idx").exists() {
        error!("Index directory {:?} does not contain the tree files", path);
        return false;
    }

    if !path.join(CHROM_FILE_NAME).exists() {
        error!(
            "Index directory {:?} does not contain {}",
            path, CHROM_FILE_NAME
        );
        return false;
    }

    if !path.join(DATASET_INFO_FILE_NAME).exists() {
        error!(
            "Index directory {:?} does not contain {}",
            path, DATASET_INFO_FILE_NAME
        );
        return false;
    }

    if !path.join(META_FILE_NAME).exists() {
        error!(
            "Index directory {:?} does not contain {}",
            path, META_FILE_NAME
        );
        return false;
    }

    true
}

pub fn line2fields(line: &str) -> Vec<String> {
    // if  \t detcted in the line , use \t to split, else use whitespace to split

    if line.contains('\t') {
        line.trim_end().split('\t').map(|s| s.to_string()).collect()
    } else {
        line.trim_end()
            .split_whitespace()
            .map(|s| s.to_string())
            .collect()
    }
}

pub fn trim_prefix(s: &str, prefix: &Option<&str>) -> String {
    if let Some(p) = prefix {
        s.trim_start_matches(p).to_string()
    } else {
        s.to_string()
    }
}

pub fn add_prefix(s: &str, prefix: &Option<&str>) -> String {
    if let Some(p) = prefix {
        format!("{}{}", p, s)
    } else {
        s.to_string()
    }
}

pub fn calc_confidence(
    evidence_arr: &Vec<u32>,
    total_size: usize,
    sample_total_evidence_vec: &Vec<u32>,
) -> f64 {
    let mut total = 0;
    let mut sorted_evidence = Vec::new();
    let mut evidence_frac_vec = Vec::new();
    let mut max_reads = 0;
    let mut positive_samples: f64 = 0.0;
    for idx in 0..total_size {
        if evidence_arr[idx] > max_reads {
            max_reads = evidence_arr[idx];
        }

        if evidence_arr[idx] > 0 {
            positive_samples += 1.0;
        }

        sorted_evidence.push(evidence_arr[idx]);
        total += evidence_arr[idx];

        // calculate CPM fraction
        if evidence_arr[idx] > 0 {
            let frac = evidence_arr[idx] as f64 / sample_total_evidence_vec[idx] as f64 * 1000000.0;
            evidence_frac_vec.push(frac.ln());
        }
    }

    // let avg_reads = evidence_arr.iter().sum::<u32>() as f64 / (n as f64);

    sorted_evidence.sort_unstable_by(|a, b| b.cmp(a)); // sort in descending order
                                                       // dbg!(&sorted_evidence);

    let mut tmp = 0;
    for (i, &e) in sorted_evidence.iter().enumerate() {
        tmp += (i + 1) as u32 * e;
    }

    let a = 2.0f64 * (tmp as f64) / (total_size as f64 * total as f64);
    let b = ((total_size + 1) as f64) / total_size as f64;

    let gini = a - b;

    if positive_samples == 0.0 {
        return 0.0;
    } else {
        return (positive_samples / total_size as f64)
            * (evidence_frac_vec.iter().sum::<f64>() / evidence_frac_vec.len() as f64).exp()
            * (1.0 - gini);
    }
}

pub fn calc_confidence_f32(
    evidence_arr: &Vec<f32>,
    total_size: usize,
    sample_total_evidence_vec: &Vec<f32>,
) -> f64 {
    let mut total = 0.0;
    let mut sorted_evidence = Vec::new();
    let mut evidence_frac_vec = Vec::new();
    let mut max_reads = 0.0;
    let mut positive_samples: f64 = 0.0;
    for idx in 0..total_size {
        if evidence_arr[idx] > max_reads {
            max_reads = evidence_arr[idx];
        }

        if evidence_arr[idx] > 0.0 {
            positive_samples += 1.0;
        }

        sorted_evidence.push(evidence_arr[idx]);
        total += evidence_arr[idx];

        // calculate CPM fraction
        if evidence_arr[idx] > 0.0 {
            let frac = evidence_arr[idx] as f64 / sample_total_evidence_vec[idx] as f64 * 1000000.0;
            evidence_frac_vec.push(frac.ln());
        }
    }

    // let avg_reads = evidence_arr.iter().sum::<u32>() as f64 / (n as f64);

    sorted_evidence.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)); // sort in descending order
                                                                                                    // dbg!(&sorted_evidence);

    let mut tmp = 0.0;
    for (i, &e) in sorted_evidence.iter().enumerate() {
        tmp += (i + 1) as f32 * e;
    }

    // dbg!(n, tmp);
    let a = 2.0f64 * (tmp as f64) / (total_size as f64 * total as f64);
    let b = ((total_size + 1) as f64) / total_size as f64;

    // dbg!(a, b);
    let gini = a - b;
    // dbg!(gini, avg_reads, max_reads);

    if positive_samples == 0.0 {
        return 0.0;
    } else {
        return (positive_samples / total_size as f64)
            * (evidence_frac_vec.iter().sum::<f64>() / evidence_frac_vec.len() as f64).exp()
            * (1.0 - gini);
    }
}

pub fn add_gz_suffix_if_needed<P: AsRef<Path>>(path: &P) -> PathBuf {
    let mut path = path.as_ref().to_path_buf();
    if let Some(ext) = path.extension() {
        if ext == "gz" {
            path.clone()
        } else {
            let ext_str = ext.to_string_lossy();
            // let path_str = path.file_stem().unwrap().to_string_lossy();
            // dbg!(&ext_str);
            // dbg!(&path_str);
            if ext_str.is_empty() {
                // dbg!(&ext_str);
                path.with_extension("gz")
            } else {
                path.set_extension(format!("{}.gz", ext_str));
                path.clone()
            }
        }
    } else {
        let mut new_path = path.clone();
        new_path.set_extension("gz");
        new_path
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_gz_suffix_if_needed() {
        let path1 = PathBuf::from("file.txt");
        let path2 = PathBuf::from("file.txt.gz");
        let path3 = PathBuf::from("file");
        let path4 = PathBuf::from("file.");

        let new_path1 = add_gz_suffix_if_needed(&path1);
        let new_path2 = add_gz_suffix_if_needed(&path2);
        let new_path3 = add_gz_suffix_if_needed(&path3);
        let new_path4 = add_gz_suffix_if_needed(&path4);

        dbg!(&new_path1);
        dbg!(&new_path2);
        dbg!(&new_path3);
        dbg!(&new_path4);

        // assert_eq!(new_path1, PathBuf::from("file.txt.gz"));
        // assert_eq!(new_path2, PathBuf::from("file.txt.gz"));
        // assert_eq!(new_path3, PathBuf::from("file.gz"));
    }
}

pub fn log_mem_stats() {
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        for line in status.lines() {
            if line.starts_with("VmRSS") || line.starts_with("VmData") || line.starts_with("VmSwap")
            {
                eprintln!("{}", line);
            }
        }
    }
}

pub fn intersect_sorted<T: Ord + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let mut i = 0;
    let mut j = 0;
    let mut result = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                result.push(a[i].clone());
                i += 1;
                j += 1;
            }
        }
    }
    result
}

pub trait GetMemSize {
    fn get_mem_size(&self) -> usize;
}
