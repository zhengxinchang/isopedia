use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::utils::{self, trim_chr_prefix_to_upper};
#[derive(Debug, Clone)]
pub struct SpliceBreakPointPair {
    pub left_chr: String,
    pub left_pos: u64,
    pub right_chr: String,
    pub right_pos: u64,
    pub id: String,
    pub rest_info: String,
}

impl SpliceBreakPointPair {
    pub fn parse_bed_line(s: &str) -> Result<Self> {
        let mut parts = s.split('\t');
        let left_chr = parts
            .next()
            .ok_or_else(|| anyhow!("Missing left_chr"))?
            .to_string();
        let left_chr = trim_chr_prefix_to_upper(&left_chr);

        let left_pos = parts
            .next()
            .ok_or_else(|| anyhow!("Missing left_pos"))?
            .parse()?;
        let right_chr = parts
            .next()
            .ok_or_else(|| anyhow!("Missing right_chr"))?
            .to_string();

        let right_chr = trim_chr_prefix_to_upper(&right_chr);
        let right_pos = parts
            .next()
            .ok_or_else(|| anyhow!("Missing right_pos"))?
            .parse()?;
        let id = parts
            .next()
            .ok_or_else(|| anyhow!("Missing id"))?
            .to_string();
        let rest_info = parts.collect::<Vec<&str>>().join("\t");

        Ok(Self {
            left_chr,
            left_pos,
            right_chr,
            right_pos,
            id,
            rest_info,
        })
    }

    pub fn parse_string(s: &str) -> Result<Self> {
        let (left_chr, left_pos, right_chr, right_pos) = utils::parse_breakpoint_str(s)?;

        Ok(Self {
            left_chr,
            left_pos,
            right_chr,
            right_pos,
            id: "single_query".to_string(),
            rest_info: "".to_string(),
        })
    }

    // pub fn parse_string_sj(s: &str) -> Result<Self> {
    //     let (chr, left_pos, right_pos) = utils::parse_splice_junction_str(s)?;

    //     Ok(Self {
    //         left_chr: chr.clone(),
    //         left_pos,
    //         right_chr: chr,
    //         right_pos,
    //         id: "single_query".to_string(),
    //         rest_info: "".to_string(),
    //     })
    // }

    pub fn is_same_chr(&self) -> bool {
        self.left_chr == self.right_chr
    }

    pub fn sort(&mut self) {
        if self.left_pos > self.right_pos {
            let t = self.left_pos;
            self.left_pos = self.right_pos;
            self.right_pos = t;
        }
    }

    pub fn to_pos_vec(&self) -> Vec<(String, u64)> {
        vec![
            (self.left_chr.clone(), self.left_pos),
            (self.right_chr.clone(), self.right_pos),
        ]
    }
}

pub fn parse_bed_to_breakpoint_pairs<P: AsRef<Path>>(p: P) -> Result<Vec<SpliceBreakPointPair>> {
    let file = File::open(p)?;
    let reader = BufReader::new(file);
    let mut breakpoints: Vec<SpliceBreakPointPair> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue; // Skip empty lines and comments
        }
        let mut bp = SpliceBreakPointPair::parse_bed_line(&line)?;
        bp.sort();
        breakpoints.push(bp);
    }

    Ok(breakpoints)
}
