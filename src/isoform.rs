//! aggregation of metainformation from multiple samples.
//!
//!
//!
// use core::hash;
use std::io::Read;

use flate2::bufread::GzEncoder;
// use itertools::Itertools;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;

use crate::{
    constants::MAX_SAMPLE_SIZE,
    reads::{AggrRead, Segment, Strand},
    utils::{self, hash_vec},
};
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ReadDiffSlim {
    pub left: u64,
    pub right: u64,
    pub supp_seg_vec_offset: u32,
    pub supp_seg_vec_length: u32,
    pub strand: Strand,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum RecordType {
    UnK,
    RefGuided,
}

#[serde_as]
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MergedIsoform {
    pub signature: u64,
    pub total_evidence: u32,
    pub sample_size: u32,
    #[serde_as(as = "[_; 64]")]
    //Records the number of reads (evidence) from each sample that support the merged isoform.
    pub sample_evidence_arr: [u32; MAX_SAMPLE_SIZE],
    #[serde_as(as = "[_; 64]")]
    // Stores the starting index of each sample’s read-level differences (ReadDiffSlim) within the shared isoform_diffs_slim_vec array.
    pub sample_offset_arr: [u32; MAX_SAMPLE_SIZE],
    pub chrom: String,
    pub rec_type: RecordType, // indicate the status of the aggr isoform, also reserved for SV.
    pub splice_junctions_vec: Vec<(u64, u64)>,
    pub isoform_diffs_slim_vec: Vec<ReadDiffSlim>,
    pub supp_segs_vec: Vec<Segment>,
}

impl PartialEq for MergedIsoform {
    fn eq(&self, other: &Self) -> bool {
        self.signature == other.signature
    }
}

impl Eq for MergedIsoform {}

impl PartialOrd for MergedIsoform {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.signature.cmp(&other.signature))
    }
}

impl Ord for MergedIsoform {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.signature.cmp(&other.signature)
    }
}

impl MergedIsoform {
    pub fn init(aggr_isofrom: AggrRead, sample_size: usize, sample_idx: u32) -> MergedIsoform {
        let mut aggr_record = MergedIsoform {
            signature: aggr_isofrom.signature,
            total_evidence: aggr_isofrom.evidence,
            sample_size: sample_size as u32, // the number of samples
            sample_evidence_arr: [0; MAX_SAMPLE_SIZE],
            sample_offset_arr: [0; MAX_SAMPLE_SIZE],
            chrom: aggr_isofrom.chrom.clone(),
            rec_type: RecordType::UnK,
            splice_junctions_vec: Vec::new(),
            isoform_diffs_slim_vec: Vec::new(),
            supp_segs_vec: Vec::new(),
        };

        //  evidence number can be  use to slice the evidence of this sample from aggr_isofrom_delta_slim_arr
        aggr_record.sample_evidence_arr[sample_idx as usize] = aggr_isofrom.evidence;
        // the offset of the sample in the aggr_isofrom_delta_slim_arr
        // first sample is 0.
        aggr_record.sample_offset_arr[sample_idx as usize] = 0;
        aggr_record.splice_junctions_vec = aggr_isofrom.splice_junctions.clone();

        // build supp_segs_vec and isoform_diffs_slim_vec
        let mut supp_segs_vec: Vec<Segment> = Vec::new();
        let aggr_isofrom_delta_slim_arr_tmp: Vec<ReadDiffSlim> = aggr_isofrom
            .read_diffs
            .iter()
            .map(|isoform_delta| {
                supp_segs_vec.extend(isoform_delta.supp_segs.clone()); // build the supp_segs_vec
                ReadDiffSlim {
                    left: isoform_delta.left,
                    right: isoform_delta.right,
                    strand: isoform_delta.strand.clone(),
                    supp_seg_vec_offset: 0,
                    supp_seg_vec_length: isoform_delta.supp_segs.len() as u32,
                }
            })
            .collect();

        // fill the offset of each isoform_delta_slim
        let mut aggr_isofrom_delta_slim_arr: Vec<ReadDiffSlim> =
            Vec::with_capacity(aggr_isofrom_delta_slim_arr_tmp.len());
        let mut curr_offset = 0;
        aggr_isofrom_delta_slim_arr_tmp
            .into_iter()
            .for_each(|mut isoform_delta_slim| {
                isoform_delta_slim.supp_seg_vec_offset = curr_offset;
                curr_offset = curr_offset + isoform_delta_slim.supp_seg_vec_length;
                aggr_isofrom_delta_slim_arr.push(isoform_delta_slim);
            });

        // update supp_segs_vec and isoform_diffs_slim_vec
        aggr_record.isoform_diffs_slim_vec = aggr_isofrom_delta_slim_arr;
        aggr_record.supp_segs_vec = supp_segs_vec;

        aggr_record
    }

    pub fn add(&mut self, aggr_isofrom: AggrRead, sample_idx: u32) {
        // update the total evidence
        self.total_evidence += aggr_isofrom.evidence;

        // self.sample_size += 1;
        // update the sample evidence
        self.sample_evidence_arr[sample_idx as usize] = aggr_isofrom.evidence;
        // update the sample offset
        self.sample_offset_arr[sample_idx as usize] = self.isoform_diffs_slim_vec.len() as u32;

        // build supp_segs_vec and isoform_diffs_slim_vec
        let mut supp_segs_vec: Vec<Segment> = Vec::new();
        let aggr_isofrom_delta_slim_arr_tmp: Vec<ReadDiffSlim> = aggr_isofrom
            .read_diffs
            .iter()
            .map(|isoform_delta| {
                supp_segs_vec.extend(isoform_delta.supp_segs.clone()); // build the supp_segs_vec
                ReadDiffSlim {
                    left: isoform_delta.left,
                    right: isoform_delta.right,
                    strand: isoform_delta.strand.clone(),
                    supp_seg_vec_offset: 0,
                    supp_seg_vec_length: isoform_delta.supp_segs.len() as u32,
                }
            })
            .collect();

        // udpate the offset of each isoform_delta_slim
        // curr_offset is the offset of the last isoform_delta_slim
        let mut curr_offset = self.supp_segs_vec.len() as u32;
        aggr_isofrom_delta_slim_arr_tmp
            .into_iter()
            .for_each(|mut isoform_delta_slim| {
                isoform_delta_slim.supp_seg_vec_offset = curr_offset;
                curr_offset = curr_offset + isoform_delta_slim.supp_seg_vec_length;
                self.isoform_diffs_slim_vec.push(isoform_delta_slim);
            });

        // update supp_segs_vec
        self.supp_segs_vec.extend(supp_segs_vec);
    }

    pub fn get_string(&self) -> String {
        let mut record = String::new();

        record.push_str(&self.signature.to_string());
        record.push('\t');
        record.push_str(&self.chrom);
        record.push('\t');
        record.push_str(self.sample_size.to_string().as_str());
        record.push('\t');
        record.push_str(&self.total_evidence.to_string());
        record.push('\t');
        record.push_str(
            &self
                .sample_evidence_arr
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<String>>()
                .join(":"),
        );
        record.push('\t');
        record.push_str(
            &self
                .sample_offset_arr
                .iter()
                .map(|o| o.to_string())
                .collect::<Vec<String>>()
                .join(":"),
        );
        record.push('\t');
        // record.push_str(&self.chrom);
        record.push('\t');
        record.push_str(stringify!(RecordType::ISOFORM));
        record.push('\t');
        record.push_str(
            &self
                .splice_junctions_vec
                .iter()
                .map(|(l, r)| format!("{}-{}", l, r))
                .collect::<Vec<String>>()
                .join(","),
        );
        record.push('\t');

        record.push_str(
            &self
                .isoform_diffs_slim_vec
                .iter()
                .map(|delta| {
                    format!(
                        "{}-{}:{}:{}:{}",
                        delta.left,
                        delta.right,
                        delta.strand.to_string(),
                        delta.supp_seg_vec_offset,
                        delta.supp_seg_vec_length
                    )
                })
                .collect::<Vec<String>>()
                .join(","),
        );
        record.push('\t');
        record.push_str(
            &self
                .supp_segs_vec
                .iter()
                .map(|seg| format!("{}%{}%{}%{}", seg.chrom, seg.start, seg.end, seg.strand))
                .collect::<Vec<String>>()
                .join(","),
        );

        record
    }

    pub fn gz_encode(&self, buf: &mut Vec<u8>) -> u32 {
        // *buf = bincode::serialize(&self).expect("Can not serialize the AggRecord");
        // dbg!(&buf);
        let tmp_buf = bincode::serialize(&self).expect("Can not serialize the AggRecord");

        let length = GzEncoder::new(tmp_buf.as_slice(), flate2::Compression::default())
            .read_to_end(buf)
            .expect("Can not compress the AggRecord");

        return length as u32;
    }

    pub fn gz_decode(bytes: &[u8]) -> Result<MergedIsoform, bincode::Error> {
        let mut decoder = flate2::bufread::GzDecoder::new(bytes);
        let mut bytes = Vec::new();
        decoder
            .read_to_end(&mut bytes)
            .expect("Can not decompress the AggRecord");

        bincode::deserialize(&bytes[..])
    }

    pub fn get_common_splice_sites(&self) -> Vec<u64> {
        self.splice_junctions_vec
            .iter()
            .flat_map(|&(a, b)| [a, b])
            .collect()
    }

    // return the left and right positions of each isoform diff
    pub fn get_read_ref_span_vec(&self) -> Vec<u64> {
        self.isoform_diffs_slim_vec
            .iter()
            .flat_map(|d| [d.left, d.right])
            .collect()
    }

    pub fn get_positive_count(&self, min_evidence: &u32) -> u32 {
        let mut c = 0;
        for i in 0..self.sample_size as usize {
            if self.sample_evidence_arr[i] >= *min_evidence {
                c += 1;
            }
        }
        c
    }

    pub fn get_sample_evidence_arr(&self) -> Vec<u32> {
        let arr = self.sample_evidence_arr.to_vec();
        arr[0..self.sample_size as usize].to_vec()
    }

    // return the inner exon positions, start from exon2 to exon(n-1)
    pub fn get_inner_exon_positions(&self) -> Vec<(u64, u64)> {
        if self.splice_junctions_vec.len() == 1 {
            return self.splice_junctions_vec.clone();
        } else {
            self.splice_junctions_vec
                .windows(2)
                .map(|window| (window[0].1, window[1].0))
                .collect()
        }
    }

    
    pub fn find_fusion(&self, chrom: &str, pos: u64, flank: u64) -> Vec<u32> {
        let mut fusion_evidence_vec = vec![0u32; MAX_SAMPLE_SIZE];

        // for each sampple
        for (idx, (offset, size)) in self
            .sample_offset_arr
            .iter()
            .zip(self.sample_evidence_arr.iter())
            .enumerate()
        {
            // collect the suppvec of each read
            if *size > 0 {
                // dbg!(idx, offset, size);
                let start = *offset as usize;
                let end = start + *size as usize;
                let diffs = &self.isoform_diffs_slim_vec[start..end];
                for diff in diffs {
                    
                    let supp_vec = &self.supp_segs_vec[diff.supp_seg_vec_offset as usize
                        ..(diff.supp_seg_vec_offset + diff.supp_seg_vec_length) as usize];
                    for seg in supp_vec {
                        // dbg!(seg);
                        // let chrom = utils::pad_chrom_prefix(chrom);
                        if seg.chrom == chrom
                            && (seg.start <= pos + flank && seg.end >= pos.saturating_sub(flank))
                        {
                            fusion_evidence_vec[idx] += 1;
                            // dbg!(self.signature);
                            break;
                        }
                    }
                }
            }
        }
        fusion_evidence_vec
    }
}
