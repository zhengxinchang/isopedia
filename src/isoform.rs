use std::{io::Read, vec};

use flate2::bufread::GzEncoder;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;

use crate::{
    // constants::MAX_SAMPLE_SIZE,
    breakpoints::BreakPointPair,
    dataset_info::DatasetInfo,
    fusion::{FusionAggrReads, FusionSingleRead},
    reads::{AggrRead, Segment, Strand},
    utils::{self, calc_cpm},
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
    // Records the number of reads (evidence) from each sample that support the merged isoform.
    // the evidence number is also the size of read_diffs_slim_vec for this sample.
    pub sample_evidence_arr: Vec<u32>,
    // Stores the starting index of each sample’s read-level differences (ReadDiffSlim) within the shared isoform_diffs_slim_vec array.
    pub sample_offset_arr: Vec<u32>,
    pub chrom: String,
    pub rec_type: RecordType, // indicate the status of the aggr isoform, also reserved for SV.
    pub splice_junctions_vec: Vec<(u64, u64)>,
    pub isoform_reads_slim_vec: Vec<ReadDiffSlim>,
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
            sample_evidence_arr: vec![0; sample_size],
            sample_offset_arr: vec![0; sample_size],
            chrom: aggr_isofrom.chrom.clone(),
            rec_type: RecordType::UnK,
            splice_junctions_vec: Vec::new(),
            isoform_reads_slim_vec: Vec::new(),
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
        aggr_record.isoform_reads_slim_vec = aggr_isofrom_delta_slim_arr;
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
        self.sample_offset_arr[sample_idx as usize] = self.isoform_reads_slim_vec.len() as u32;

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
                self.isoform_reads_slim_vec.push(isoform_delta_slim);
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
                .isoform_reads_slim_vec
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
            .expect("Can not decompress the AggRecord, consider the index is corrupted?");

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
        self.isoform_reads_slim_vec
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

    pub fn get_positive_array(&self, min_evidence: &u32) -> Vec<u32> {
        self.sample_evidence_arr
            .iter()
            .map(|e| if *e >= *min_evidence { 1 } else { 0 })
            .collect()
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

    pub fn find_fusion_by_breakpoints(
        &self,
        chrom: &str,
        pos: u64,
        flank: u64,
        dbinfo: &DatasetInfo,
    ) -> Vec<u32> {
        let mut fusion_evidence_vec = vec![0u32; dbinfo.get_size()];

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
                let diffs = &self.isoform_reads_slim_vec[start..end];
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

    /// Each isoform will be re-scattered at read level, and regrouped by their fusion hash.
    pub fn to_fusion_candidates(&self) -> Option<Vec<FusionAggrReads>> {
        let mut candidates: Vec<FusionSingleRead> = Vec::new();
        let mut aggr_candidates: std::collections::HashMap<
            u64,
            FusionAggrReads,
            rustc_hash::FxBuildHasher,
        > = FxHashMap::default();

        for sample_idx in 0..self.sample_size as usize {
            // process all read from a sample
            let sample_evidence = self.sample_evidence_arr[sample_idx] as usize;
            if sample_evidence == 0 {
                continue; // skip samples with no evidence
            }
            let sample_evidence_offset = self.sample_offset_arr[sample_idx] as usize;

            // process each read in the sample
            let read_diffs = &self.isoform_reads_slim_vec
                [sample_evidence_offset..sample_evidence_offset + sample_evidence];

            for read_diff in read_diffs {
                // skip reads with no supp segments
                if read_diff.supp_seg_vec_length == 0 {
                    continue; // skip reads with no supp segments
                }

                let supp_segments = &self.supp_segs_vec[read_diff.supp_seg_vec_offset as usize
                    ..(read_diff.supp_seg_vec_offset + read_diff.supp_seg_vec_length) as usize];

                // skip reads with supp segments on different chromosomes
                let unique_chr = if let Some(first_chr) = supp_segments.first().map(|o| &o.chrom) {
                    supp_segments.iter().all(|o| &o.chrom == first_chr)
                } else {
                    true
                };

                if !unique_chr {
                    continue; // skip reads with supp segments on different chromosomes
                }

                // create a FusionSingleRead for this read
                let fusion_read = FusionSingleRead::new(
                    self.chrom.clone(),
                    supp_segments.last().unwrap().chrom.clone(),
                    sample_idx as u32,
                    supp_segments,
                    &self.splice_junctions_vec,
                );

                candidates.push(fusion_read);
            }
        }

        for candidate in candidates {
            // check if the candidate already exists in the aggr_candidates map
            if let Some(existing) = aggr_candidates.get_mut(&candidate.fusion_hash) {
                existing.add(&candidate);
            } else {
                // create a new FusionAggrReads from the candidate
                let new_aggr_read = FusionAggrReads::init(&candidate);
                aggr_candidates.insert(candidate.fusion_hash, new_aggr_read);
            };
        }

        if aggr_candidates.is_empty() {
            return None; // no valid candidates found
        } else {
            // convert the aggr_candidates map into a vector
            let fusion_candidates: Vec<FusionAggrReads> = aggr_candidates.into_values().collect();
            Some(fusion_candidates)
        }
    }

    /// Function to get the confidence value of the isoform based on the evidence
    pub fn get_confidence_value(evidence_arr: Vec<u32>, dataset_info: &DatasetInfo) -> f64 {
        let n = dataset_info.get_size();
        let mut total = 0;
        let mut sorted_evidence = Vec::new();
        let mut evidence_frac_vec = Vec::new();
        let mut max_reads = 0;
        let mut positive_samples: f64 = 0.0;
        for idx in 0..n {
            if evidence_arr[idx] > max_reads {
                max_reads = evidence_arr[idx];
            }

            if evidence_arr[idx] > 0 {
                positive_samples += 1.0;
            }

            sorted_evidence.push(evidence_arr[idx]);
            total += evidence_arr[idx];

            if evidence_arr[idx] > 0 {
                let frac = evidence_arr[idx] as f64
                    / dataset_info.sample_total_evidence_vec[idx] as f64
                    * 1000000.0;
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

        // dbg!(n, tmp);
        let a = 2.0f64 * (tmp as f64) / (n as f64 * total as f64);
        let b = ((n + 1) as f64) / n as f64;

        // dbg!(a, b);
        let gini = a - b;
        // dbg!(gini, avg_reads, max_reads);

        if positive_samples == 0.0 {
            return 0.0;
        } else {
            return (positive_samples / n as f64)
                * (evidence_frac_vec.iter().sum::<f64>() / evidence_frac_vec.len() as f64).exp()
                * (1.0 - gini);
        }
    }

    /// get partial report for splice junction searching
    pub fn get_splice_report(
        &self,
        bpp: &BreakPointPair,
        flank: u64,
        dbinfo: &DatasetInfo,
    ) -> Option<String> {
        let mut out_str = Vec::new();

        // support evidence
        out_str.push(self.total_evidence.to_string());

        // cpm of this isoform
        out_str.push(
            calc_cpm(
                &self.total_evidence,
                &(dbinfo.sample_total_evidence_vec.iter().sum::<u32>()),
            )
            .to_string(),
        );

        // index of the SJ in the isoform matched with query
        if let Some((idx, sj)) = self.splice_junctions_vec.iter().enumerate().find(|(_, x)| {
            if (x.0.abs_diff(bpp.left_pos) < flank) && (x.1.abs_diff(bpp.right_pos) < flank) {
                return true;
            } else {
                return false;
            }
        }) {
            // idx of splice junction
            out_str.push((idx + 1).to_string());

            // distance to splice sites
            // dbg!(sj.0, bpp.left_pos, sj.1, bpp.right_pos);

            out_str.push(format!(
                "{},{}",
                utils::u64diff2i32(sj.0, bpp.left_pos),
                utils::u64diff2i32(sj.1, bpp.right_pos)
            ));
        } else {
            return None;
        }

        // for the record only have one SJ, it could be a mono-exon or two exons. it needs to be confirmed here.
        if self.splice_junctions_vec.len() == 1 {
            if self.isoform_reads_slim_vec[0].left == self.splice_junctions_vec[0].0
                && self.isoform_reads_slim_vec[0].right == self.splice_junctions_vec[0].1
            {
                // mono-exon
                out_str.push("1".to_string());
            } else {
                out_str.push("2".to_string());
            }
        } else {
            out_str.push((self.splice_junctions_vec.len() + 1).to_string());
        }

        let mut starts = Vec::with_capacity(10);
        let mut ends = Vec::with_capacity(10);
        let mut sample_vec = Vec::with_capacity(dbinfo.get_size());
        let mut sample_idx = 0usize;

        for (offset, length) in self
            .sample_offset_arr
            .iter()
            .zip(self.sample_evidence_arr.iter())
        {
            // for a particular sample, do:
            if *length > 0 {
                let sub = &self.isoform_reads_slim_vec
                    [*offset as usize..(*offset as usize + *length as usize)];
                let mut sample_read_sub = Vec::with_capacity(*length as usize);
                for read in sub {
                    starts.push(read.left);
                    ends.push(read.right);

                    sample_read_sub.push(format!("{}|{}|{}", read.left, read.right, read.strand));
                }
                let sample_read_sub = sample_read_sub.join(",");
                sample_vec.push(format!(
                    "{}:{}:{}",
                    length,
                    calc_cpm(
                        length,
                        &(dbinfo.sample_total_evidence_vec[sample_idx] as u32)
                    ),
                    sample_read_sub
                ));
            } else {
                sample_vec.push(format!("{}:{}:{}", 0, 0, "NULL"));
            }
            sample_idx = sample_idx + 1;
        }

        //left right most numbers
        out_str.push(starts.iter().min().unwrap().to_string());
        out_str.push(starts.iter().max().unwrap().to_string());
        out_str.push(ends.iter().min().unwrap().to_string());
        out_str.push(ends.iter().max().unwrap().to_string());

        // splice junctions

        out_str.push(
            self.splice_junctions_vec
                .iter()
                .map(|(l, r)| format!("{}-{}", l, r))
                .collect::<Vec<String>>()
                .join(","),
        );

        // format

        out_str.push("COUNT:CPM:START,END,STRAND".to_string());

        out_str.extend_from_slice(&sample_vec);

        Some(out_str.join("\t"))
    }
}

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_get_confidence_value() {
        let evidence_arr = vec![999, 0, 0, 0];
        let mut dataset_info = DatasetInfo::new();
        dataset_info.sample_size = evidence_arr.len();
        dataset_info.sample_total_evidence_vec = vec![10000, 20000, 30000, 40000];

        let confidence = MergedIsoform::get_confidence_value(evidence_arr, &dataset_info);

        dbg!(confidence);
    }
}
