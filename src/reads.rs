use crate::utils::{self, hash_vec};
use bio_types::strand::ReqStrand;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

impl PartialEq for Strand {
    fn ne(&self, other: &Self) -> bool {
        match self {
            Strand::Plus => match other {
                Strand::Plus => false,
                Strand::Minus => true,
            },
            Strand::Minus => match other {
                Strand::Plus => true,
                Strand::Minus => false,
            },
        }
    }

    fn eq(&self, other: &Self) -> bool {
        core::mem::discriminant(self) == core::mem::discriminant(other)
    }
}

impl Strand {
    pub fn from_string(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }
    pub fn to_string(&self) -> String {
        match self {
            Strand::Plus => "+".to_string(),
            Strand::Minus => "-".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Segment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand, // 0 for +, 1 for -
    pub is_supp: bool,
} // using segment instead of single SJ is more informative and save memory(only one chrom and strand for all SJ)

impl Segment {
    pub fn new(chrom: String, start: u64, end: u64, strand: Strand, is_supp: bool) -> Self {
        let chrom2 = utils::trim_chr_prefix_to_upper(&chrom);

        Self {
            chrom: chrom2,
            start,
            end,
            strand,
            is_supp,
        }
    }

    pub fn cmp_end(&self, other: &Self) -> std::cmp::Ordering {
        if self.end < other.end {
            return std::cmp::Ordering::Less;
        } else if self.end > other.end {
            return std::cmp::Ordering::Greater;
        } else {
            return std::cmp::Ordering::Equal;
        }
    }
}

impl PartialEq for Segment {
    fn ne(&self, other: &Self) -> bool {
        if self.chrom == other.chrom
            && self.start == other.start
            && self.end == other.end
            && self.strand == other.strand
            && self.is_supp == other.is_supp
        {
            return false;
        } else {
            return true;
        }
    }

    fn eq(&self, other: &Self) -> bool {
        !self.ne(other)
    }
}

/// Implement the PartialOrd trait for Segment
/// Only compare the start position of the segment
/// This should be fine if the all segments are from the same cigar string
/// If the segments from supplementary alignment are mixed with the segments
/// from the primary alignment, it might have problem when the supplementary
/// alignment is overlap with the primary alignment.
/// # Arguments
/// * `other` - The other segment to compare
/// # Returns
/// * `Option<std::cmp::Ordering>` - The ordering of the two segments
impl PartialOrd for Segment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if !self.is_supp && other.is_supp {
            // make sure the primary alignment is always before the supplementary alignment
            return Some(std::cmp::Ordering::Less);
        } else if self.is_supp && !other.is_supp {
            return Some(std::cmp::Ordering::Greater);
        } else {
            if self.chrom < other.chrom {
                return Some(std::cmp::Ordering::Less);
            } else if self.chrom > other.chrom {
                return Some(std::cmp::Ordering::Greater);
            } else {
                if self.start < other.start {
                    return Some(std::cmp::Ordering::Less);
                } else if self.start > other.start {
                    return Some(std::cmp::Ordering::Greater);
                } else {
                    Some(self.cmp_end(other))
                }
            }
        }
    }
}

#[derive(Debug)]
#[allow(unused)]
pub struct SingleRead {
    pub chrom: String,
    pub seg_size: u64,
    pub read_len: u64,
    pub read_start: u64,
    pub ref_span: u64, // not sure if it is necessary
    pub mapq: u8,
    pub evidence: u32,
    pub segment_list: Vec<Segment>,
    pub pri_seg_size: u64,
    pub signature: u64,
    pub strand: Strand,
}

impl SingleRead {
    pub fn new(chrom: String, mapq: u8, supp_raw_read: u32, readlen: u64, start: u64) -> Self {
        Self {
            chrom: chrom,
            seg_size: 0,
            segment_list: Vec::new(),
            mapq: mapq,
            evidence: supp_raw_read,
            read_len: readlen,
            read_start: start,
            ref_span: 0,
            pri_seg_size: 0,
            signature: 0,
            strand: Strand::Plus,
        }
    }

    /// Add a segment to the isoform
    /// # Arguments
    /// * `chrom` - The chromosome of the segment
    /// * `start` - The start position of the segment
    /// * `end` - The end position of the segment
    /// * `strand` - The strand of the segment
    /// * `is_supp` - Whether the segment is supported
    /// the strand is a ReqStrand enum, either Forward or Reverse.
    pub fn add_segment(
        &mut self,
        chrom: String,
        start: u64,
        end: u64,
        strand: &ReqStrand,
        is_supp: bool,
    ) {
        let strand = match strand {
            ReqStrand::Forward => Strand::Plus,
            ReqStrand::Reverse => Strand::Minus,
        };

        self.segment_list
            .push(Segment::new(chrom, start, end, strand, is_supp));
        self.seg_size += 1;
        self.pri_seg_size += 1;
    }

    /// Add a segment to the isoform
    /// # Arguments
    /// * `chrom` - The chromosome of the segment
    /// * `start` - The start position of the segment
    /// * `end` - The end position of the segment
    /// * `strand` - The strand of the segment
    /// * `is_supp` - Whether the segment is supported
    /// the strand is a string, either "+" or "-".
    pub fn add_supp_segment(
        &mut self,
        chrom: String,
        start: u64,
        end: u64,
        strand: &str,
        is_supp: bool,
    ) {
        let strand = Strand::from_string(strand);

        self.segment_list
            .push(Segment::new(chrom, start, end, strand, is_supp));
        self.seg_size += 1
    }

    pub fn update_ref_span(&mut self, ref_span: u64) {
        self.ref_span = ref_span;
    }

    pub fn sort_segments(&mut self) {
        self.segment_list.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    pub fn get_pri_splice_junction(&self) -> Vec<(u64, u64)> {
        let mut res = Vec::new();
        for seg in &self.segment_list {
            if !seg.is_supp {
                res.push(seg.start);
                res.push(seg.end);
            }
        }
        // remove the first and last postion
        if res.len() > 2 {
            res.pop();
            res.remove(0);
        }
        // chunk the vector into pairs for splice junction
        res.chunks(2).into_iter().map(|x| (x[0], x[1])).collect()
    }

    pub fn get_supp_segments(&self) -> Vec<Segment> {
        self.segment_list
            .clone()
            .into_iter()
            .filter(|x| x.is_supp)
            .map(|x| x.clone())
            .collect()
    }

    /// Process the isoform
    /// the signature is the hash value of the splice junctions
    pub fn process(&mut self) {
        self.signature = hash_vec(&self.get_pri_splice_junction());
        self.sort_segments();
    }
}

// use crate::extraction::SingleRead;
// // use crate::extractor::aggr_isofrom::Segment;
// // use crate::extractor::aggr_isofrom::Strand;
// use crate::extraction::single_read::Segment;
// use crate::extraction::single_read::Strand;

#[derive(Debug, Clone)]
pub struct ReadDiff {
    pub left: u64,
    pub right: u64,
    pub strand: Strand,
    pub supp_segs: Vec<Segment>,
}

/// Aggregated isoform from single sample in extractor
/// also used in aggregator
#[derive(Debug, Clone)]
pub struct AggrRead {
    pub signature: u64,
    pub chrom: String,
    pub splice_junctions: Vec<(u64, u64)>,
    pub evidence: u32,
    pub read_diffs: Vec<ReadDiff>,
}

impl AggrRead {
    pub fn new(sr: SingleRead) -> AggrRead {
        AggrRead {
            signature: sr.signature,
            chrom: sr.chrom.clone(),
            splice_junctions: sr.get_pri_splice_junction(),
            evidence: 1,
            read_diffs: vec![ReadDiff {
                left: sr.read_start,
                right: sr.ref_span,
                strand: sr.strand,
                supp_segs: sr.get_supp_segments(),
            }],
        }
    }

    pub fn add_read(&mut self, isoform: SingleRead) {
        self.evidence += 1;
        self.read_diffs.push(ReadDiff {
            left: isoform.read_start,
            right: isoform.ref_span,
            strand: isoform.strand,
            supp_segs: isoform.get_supp_segments(),
        });
    }

    /// convert the AggIsoform to a string
    /// format of the string:
    /// 6310607345954554960 1        chr19 35889680-35890651   35889680-35890651:chr19%35889620%35889835%+
    /// TODO:use binary file instead of string
    pub fn to_record(&self) -> String {
        /*
        format: signature \t evidence \t chrom \t splice_junctions \t \t isoform_diffs
         */
        let mut record = String::new();
        record.push_str(&self.signature.to_string());
        record.push('\t');
        record.push_str(&self.evidence.to_string());
        record.push('\t');
        record.push_str(&self.chrom);
        record.push('\t');
        // the common splice junctions
        record.push_str(
            &self
                .splice_junctions
                .iter()
                .map(|(l, r)| format!("{}-{}", l, r))
                .collect::<Vec<String>>()
                .join(","),
        );
        record.push('\t');
        record.push_str(
            &self
                .read_diffs
                .iter()
                .map(|delta| {
                    // isoform diff process
                    format!(
                        "{}-{}:{}:{}",
                        delta.left,
                        delta.right,
                        delta.strand.to_string(),
                        delta
                            .supp_segs
                            .iter()
                            .map(|seg| {
                                // for each reads ,add the supplmentarty segments
                                format!("{}%{}%{}%{}", seg.chrom, seg.start, seg.end, seg.strand)
                            })
                            .collect::<Vec<String>>()
                            .join(",")
                    )
                })
                .collect::<Vec<String>>()
                .join("|"),
        );
        record.push('\n');
        record
    }

    /// create a new AggIsoform from a string
    /// format of the string:
    /// 6310607345954554960 1        chr19 35889680-35890651   35889680-35890651:chr19%35889620%35889835%+
    /// signature           evidence chrom splice_junctions    isoform_diffs
    pub fn from_string(rec: &String) -> AggrRead {
        // dbg!(&rec);
        let mut iter = rec.trim().split('\t');
        // dbg!(&iter.collect::<Vec<&str>>());
        let signature = iter.next().unwrap().parse::<u64>().unwrap();
        let evidence = iter.next().unwrap().parse::<u32>().unwrap();
        let chrom = iter.next().unwrap().to_string();
        let splice_junctions = iter
            .next()
            .unwrap()
            .split(',')
            .map(|s| {
                let mut iter2 = s.split('-');
                (
                    iter2.next().unwrap().parse::<u64>().unwrap(),
                    iter2.next().unwrap().parse::<u64>().unwrap(),
                )
            })
            .collect::<Vec<(u64, u64)>>();
        // dbg!(iter.next().unwrap().split('|').collect());
        // dbg!(&splice_junctions);
        let read_diff = iter
            .next()
            .unwrap()
            .split('|')
            .filter(|s| s.len() > 0)
            .map(|s| {
                let mut iter3 = s.split(':');
                // dbg!(&s);
                let left_right: Vec<u64> = iter3
                    .next()
                    .unwrap()
                    .split('-')
                    .map(|s| s.parse::<u64>().unwrap())
                    .collect();
                let strand = Strand::from_string(iter3.next().unwrap());
                let supp_segs = iter3
                    .next()
                    .unwrap()
                    .split(',')
                    .filter(|s| s.len() > 0)
                    .map(|s| {
                        let mut iter4 = s.split('%');
                        Segment {
                            chrom: iter4.next().unwrap().to_string(),
                            start: iter4.next().unwrap().parse::<u64>().unwrap(),
                            end: iter4.next().unwrap().parse::<u64>().unwrap(),
                            strand: Strand::from_string(iter4.next().unwrap()),
                            is_supp: true,
                        }
                    })
                    .collect::<Vec<Segment>>();
                ReadDiff {
                    left: left_right[0],
                    right: left_right[1],
                    strand: strand,
                    supp_segs,
                }
            });

        AggrRead {
            signature,
            evidence,
            chrom,
            splice_junctions,
            read_diffs: read_diff.collect(),
        }
    }

    pub fn merge(&mut self, other: AggrRead) {
        self.evidence += other.evidence;
        self.read_diffs.extend(other.read_diffs);
    }
}

/// todo: read binary file
pub struct SingleSampleReader {
    pub file_name: String,
    pub offset: u64,
    pub reader: BufReader<GzDecoder<File>>,
    pub curr_rec_string: String,
}

impl SingleSampleReader {
    // TODO: add a function to read a binary file
    pub fn new(file_path: &str) -> SingleSampleReader {
        let file = File::open(file_path).expect("Can not read file...");
        let decoder = GzDecoder::new(file);
        let reader = BufReader::with_capacity(10 * 1024 * 1024, decoder);
        let mut agg_file_reader = SingleSampleReader {
            file_name: file_path.to_string(),
            offset: 0,
            reader: reader,
            curr_rec_string: String::with_capacity(1024),
        };
        agg_file_reader.curr_rec_string.clear();
        agg_file_reader
    }

    pub fn skip_lines(&mut self, skip: u32) {
        for _ in 0..skip {
            self.reader.read_line(&mut self.curr_rec_string).unwrap();
        }
        self.curr_rec_string.clear();
    }

    pub fn next_line_str(&mut self) -> Option<String> {
        // throw a error message if read line failed
        let bytes = self.reader.read_line(&mut self.curr_rec_string);
        match bytes {
            Ok(0) => None,
            Ok(_) => Some(self.curr_rec_string.clone()),
            Err(e) => panic!("Failed to read line from file {}: {} ", self.file_name, e),
        }
    }

    pub fn next_rec(&mut self) -> Option<AggrRead> {
        let bytes = self.reader.read_line(&mut self.curr_rec_string);
        match bytes {
            Ok(0) => None,
            Ok(_) => {
                if self.curr_rec_string.len() == 0 {
                    // dbg!("empty string");
                    None
                } else {
                    let agg_isoform = AggrRead::from_string(&self.curr_rec_string);
                    self.curr_rec_string.clear();
                    Some(agg_isoform)
                }
            }
            Err(e) => panic!("Failed to read line from file {}: {}", self.file_name, e),
        }
    }
}

impl Iterator for SingleSampleReader {
    type Item = AggrRead;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_rec()
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;

    #[test]
    fn test_from_string() {
        let rec = "18421869509402655781\t1\tchr13\t93545258-93545421\t93545258-93545421:+:chr13%93226861%93227618%+,chr13%93830153%94027895%+,chr13%94286347%94403640%+\n".to_string();
        let agg_isoform = AggrRead::from_string(&rec);
        dbg!(agg_isoform);
    }
}
