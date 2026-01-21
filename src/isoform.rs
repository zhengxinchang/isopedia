use std::{io::Read, vec};

use crate::{
    // constants::MAX_SAMPLE_SIZE,
    breakpoints::SpliceBreakPointPair,
    cmd::fusion::FusionBreakPointPair,
    dataset_info::DatasetInfo,
    fusion::{FusionAggrReads, FusionSingleRead},
    myio::SampleChip,
    reads::{AggrRead, Segment, Strand},
    utils::{self, calc_cpm},
};
use flate2::bufread::GzEncoder;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ReadDiffSlim {
    pub left: u64,
    pub right: u64,
    pub supp_seg_vec_offset: u32,
    pub supp_seg_vec_length: u32,
    pub info: Vec<(String, String)>,
    pub strand: Strand,
}

impl ReadDiffSlim {
    pub fn to_string(&self) -> String {
        format!(
            "{}-{}:{}:{}:{}:{}",
            self.left,
            self.right,
            self.strand.to_string(),
            self.supp_seg_vec_offset,
            self.supp_seg_vec_length,
            self.info
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect::<Vec<String>>()
                .join(";")
        )
    }

    pub fn to_string_no_offsets(&self) -> String {
        format!(
            "{}|{}|{}|{}",
            self.left,
            self.right,
            self.strand.to_string(),
            self.info
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect::<Vec<String>>()
                .join(";")
        )
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum RecordType {
    Spliced,
    MonoExon,
    UnK,
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
    pub chrom_id: u16,
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
    pub fn init(
        aggr_isofrom: &AggrRead,
        sample_size: usize,
        sample_idx: u32,
        chrom_id: u16,
    ) -> MergedIsoform {
        let mut aggr_record = MergedIsoform {
            signature: aggr_isofrom.signature,
            total_evidence: aggr_isofrom.evidence,
            sample_size: sample_size as u32, // the number of samples
            sample_evidence_arr: vec![0; sample_size],
            sample_offset_arr: vec![0; sample_size],
            chrom: aggr_isofrom.chrom.clone(),
            chrom_id: chrom_id,
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
                    info: isoform_delta.info.clone(),
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
                    info: isoform_delta.info.clone(),
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
                .map(|delta| delta.to_string())
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
        decoder.read_to_end(&mut bytes)?;
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

    pub fn get_sample_offset_arr(&self) -> Vec<u32> {
        self.sample_offset_arr.to_vec()
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

    pub fn get_start_pos(&self) -> u64 {
        self.splice_junctions_vec
            .first()
            .map(|sj| sj.0)
            .unwrap_or(0)
    }

    pub fn is_mono_exonic(&self) -> bool {
        if self.splice_junctions_vec.len() != 1 {
            return false;
        }

        let mut is_mono_exonic = true;

        for (left, right) in self
            .isoform_reads_slim_vec
            .iter()
            .map(|d| (d.left, d.right))
        {
            if left != self.splice_junctions_vec[0].0 || right != self.splice_junctions_vec[0].1 {
                is_mono_exonic = false;
                break;
            }
        }
        is_mono_exonic
    }

    /// Find fusion evidence based on breakpoints
    /// first search one part of the breakpoint, then check the supporting segments to see if they cover the other part of the breakpoint
    pub fn check_fusion_mate_breakpoint(
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
                    // TODO:
                    // optmize to only check the fisrt and last segment
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

    /// Cast the isoform to fusion records based on the breakpoints
    /// format
    /// #1.  #2.    #3.   #4. #5.    #6.  #7.                                   #8.    
    /// chr1,start1,end1,chr2,start2,end2,left_splice_junctions(start-end,...),right_splice_junctions(start-end,...),
    /// #9.              #10...
    /// total_evidence, sample1_evidence,sample2_evidence,...
    pub fn cast_to_fusion_records(
        &self,
        chrom: &str,
        pos: u64,
        flank: u64,
        dbinfo: &DatasetInfo,
        breakpoints: &FusionBreakPointPair,
    ) -> Vec<FusionRecord> {
        let mut fusion_records: Vec<FusionRecord> = Vec::new();

        // for each sampple
        for (sample_idx, (offset, size)) in self
            .sample_offset_arr
            .iter()
            .zip(self.sample_evidence_arr.iter())
            .enumerate()
        {
            // collect the suppvec of each read
            if *size > 0 {
                // dbg!(sample_idx, offset, size);
                let start = *offset as usize;
                let end = start + *size as usize;
                let one_sample_reads = &self.isoform_reads_slim_vec[start..end];

                for single_read in one_sample_reads {
                    // skip reads with no supp segments
                    if single_read.supp_seg_vec_length == 0 {
                        continue; // skip reads with no supp segments
                    }

                    let supp_seg_vec = &self.supp_segs_vec[single_read.supp_seg_vec_offset as usize
                        ..(single_read.supp_seg_vec_offset + single_read.supp_seg_vec_length)
                            as usize];

                    let mut is_ok = false;

                    dbg!(single_read.supp_seg_vec_length);

                    if single_read.supp_seg_vec_length == 1 {
                        let seg = &supp_seg_vec[0];
                        if seg.chrom == chrom
                            && (
                                // sig start in the region
                                (seg.start <= pos + flank && seg.start >= pos.saturating_sub(flank))
                                    ||
                                    // sig end in the region
                                    (seg.end <= pos + flank && seg.end >= pos.saturating_sub(flank))
                            )
                        {
                            // report the fusion record
                            is_ok = true;
                        }
                    } else {
                        // let first_seg = &supp_seg_vec[0];
                        // let last_seg = &supp_seg_vec[supp_seg_vec.len() - 1];

                        // if first_seg.chrom == chrom
                        //     && (
                        //         // sig start in the region
                        //         (first_seg.start <= pos + flank && first_seg.start >= pos.saturating_sub(flank))
                        //             ||
                        //             // sig end in the region
                        //             (first_seg.end <= pos + flank && first_seg.end >= pos.saturating_sub(flank))
                        //     )
                        // {
                        //     is_ok = true;
                        // } else if last_seg.chrom == chrom
                        //     && (
                        //         // sig start in the region
                        //         (last_seg.start <= pos + flank && last_seg.start >= pos.saturating_sub(flank))
                        //             ||
                        //             // sig end in the region
                        //             (last_seg.end <= pos + flank && last_seg.end >= pos.saturating_sub(flank))
                        //     )
                        // {
                        //     is_ok = true;
                        // }

                        // check all segments
                        for seg in supp_seg_vec {
                            if seg.chrom == chrom
                                && (
                                    // sig start in the region
                                    (seg.start <= pos + flank
                                        && seg.start >= pos.saturating_sub(flank))
                                        ||
                                        // sig end in the region
                                        (seg.end <= pos + flank
                                            && seg.end >= pos.saturating_sub(flank))
                                )
                            {
                                // report the fusion record
                                is_ok = true;
                                break;
                            }
                        }
                    }

                    if is_ok {
                        dbg!(supp_seg_vec);

                        let sample_name = dbinfo.get_sample_names()[sample_idx].clone();
                        let fusion_record = FusionRecord::new(
                            breakpoints,
                            single_read,
                            &supp_seg_vec,
                            self,
                            sample_name,
                        );
                        // dbg!(&fusion_record);
                        fusion_records.push(fusion_record);
                    }
                }
            }
        }

        fusion_records
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
    pub fn get_confidence_value(
        evidence_arr: &Vec<u32>,
        total_size: usize,
        sample_total_evidence_vec: &Vec<u32>,
    ) -> f64 {
        utils::calc_confidence(evidence_arr, total_size, sample_total_evidence_vec)
    }

    /// get partial report for splice junction searching
    pub fn get_splice_report(
        &self,
        bpp: &SpliceBreakPointPair,
        flank: u64,
        dbinfo: &DatasetInfo,
    ) -> Option<(Vec<String>, Vec<SampleChip>)> {
        let mut result_buffer = Vec::new();

        let mut sample_chip_vec = Vec::new();

        // support evidence
        result_buffer.push(self.total_evidence.to_string());

        // cpm of this isoform
        result_buffer.push(
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
            result_buffer.push((idx + 1).to_string());

            // distance to splice sites
            // dbg!(sj.0, bpp.left_pos, sj.1, bpp.right_pos);

            result_buffer.push(format!(
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
                result_buffer.push("1".to_string());
            } else {
                result_buffer.push("2".to_string());
            }
        } else {
            result_buffer.push((self.splice_junctions_vec.len() + 1).to_string());
        }

        let mut starts = Vec::with_capacity(10);
        let mut ends = Vec::with_capacity(10);
        let mut sample_vec = Vec::with_capacity(dbinfo.get_size());
        let mut sample_idx = 0usize;

        for (offset, read_count) in self
            .sample_offset_arr
            .iter()
            .zip(self.sample_evidence_arr.iter())
        {
            // for a particular sample, do:
            let mut sample_chip = SampleChip::default(None);

            if *read_count > 0 {
                let sub = &self.isoform_reads_slim_vec
                    [*offset as usize..(*offset as usize + *read_count as usize)];
                let mut sample_read_sub = Vec::with_capacity(*read_count as usize);
                for read in sub {
                    starts.push(read.left);
                    ends.push(read.right);
                    sample_read_sub.push(format!("{}|{}|{}", read.left, read.right, read.strand));
                }
                let sample_read_sub = sample_read_sub.join(","); // example: 2:0.3399716769595925:7668421|7687507|+,7668421|7687489|+

                // sample_chip

                let cpm = calc_cpm(
                    read_count,
                    &(dbinfo.sample_total_evidence_vec[sample_idx] as u32),
                );

                sample_vec.push(format!("{}:{}:{}", &read_count, cpm, sample_read_sub));

                sample_chip.add_item(&read_count.to_string());
                sample_chip.add_item(&cpm.to_string());
                sample_chip.add_item(&sample_read_sub);
            } else {
                sample_vec.push(format!("{}:{}:{}", 0, 0, "NULL"));

                sample_chip.add_item(&"0");
                sample_chip.add_item(&"0");
                sample_chip.add_item(&"NULL");
            }
            sample_idx = sample_idx + 1;
            sample_chip_vec.push(sample_chip);
        }

        //left right most numbers
        result_buffer.push(starts.iter().min().unwrap().to_string());
        result_buffer.push(starts.iter().max().unwrap().to_string());
        result_buffer.push(ends.iter().min().unwrap().to_string());
        result_buffer.push(ends.iter().max().unwrap().to_string());

        result_buffer.push(
            self.splice_junctions_vec
                .iter()
                .map(|(l, r)| format!("{}-{}", l, r))
                .collect::<Vec<String>>()
                .join(","),
        );

        // format

        // result_buffer.push("COUNT:CPM:START,END,STRAND".to_string());

        // result_buffer.extend_from_slice(&sample_vec);

        Some((result_buffer, sample_chip_vec))
    }
}

/// Fusion record for query two regions for example BCR exon 13 and ABL1 exon 2
/// this record is represent a single read fusion evidence
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct FusionRecord {
    pub fusion_hash: u64,
    pub left_chrom: String,
    pub left_start: u64,
    pub left_end: u64,
    pub right_chrom: String,
    pub right_start: u64,
    pub right_end: u64,
    pub left_exons_junctions: Vec<u64>, // exon positions on the left side not splice junctions
    pub right_exons_junctions: Vec<u64>, // exon positions on the right side not splice junctions
    pub sample_name: String,
}

impl FusionRecord {
    pub fn new(
        breakpoints: &FusionBreakPointPair, // query breakpoint info
        read_diff: &ReadDiffSlim,           // read start/end
        supp_segs: &[Segment],              // supplementary segments of the read
        misoform: &MergedIsoform,           // common splice junction isoform
        sample_name: String,
    ) -> FusionRecord {
        // let left_sj_vec = misoform.splice_junctions_vec;
        // let mut all_sj_vec = Vec::new();
        let mut left_exons = Vec::new();
        let mut right_exons = Vec::new();
        // cat left sj and right sj
        left_exons.push(read_diff.left);
        for &(l, r) in misoform.splice_junctions_vec.iter() {
            // all_sj_vec.push(l);
            // all_sj_vec.push(r);
            left_exons.push(l);
            left_exons.push(r);
        }
        left_exons.push(read_diff.right);

        if supp_segs.len() == 1 {
            // mono-exon fusion read
            right_exons.push(supp_segs[0].start);
            right_exons.push(supp_segs[0].end);
            // all_sj_vec.push(supp_segs[0].start);
            // all_sj_vec.push(supp_segs[0].end);
        } else {
            // multi-exon fusion read
            // right_exons.push(supp_segs[0].start);
            for sg in supp_segs.iter() {
                right_exons.push(sg.start);
                right_exons.push(sg.end);
                // all_sj_vec.push(sg.start);
                // all_sj_vec.push(sg.end);
            }
        }

        // let fusion_hash = utils::hash_vec(&all_sj_vec);

        FusionRecord {
            fusion_hash: 0,
            left_chrom: breakpoints.left_chr.clone(),
            left_start: breakpoints.left_start,
            left_end: breakpoints.left_end,
            right_chrom: breakpoints.right_chr.clone(),
            right_start: breakpoints.right_start,
            right_end: breakpoints.right_end,
            left_exons_junctions: left_exons,
            right_exons_junctions: right_exons,
            sample_name: sample_name,
        }
    }

    pub fn get_string(&self) -> String {
        let mut out = String::new();
        out.push_str(&self.left_chrom);
        out.push('\t');
        out.push_str(&self.left_start.to_string());
        out.push('\t');
        out.push_str(&self.left_end.to_string());
        out.push('\t');
        out.push_str(&self.right_chrom);
        out.push('\t');
        out.push_str(&self.right_start.to_string());
        out.push('\t');
        out.push_str(&self.right_end.to_string());
        out.push('\t');
        out.push_str(&(self.left_exons_junctions.len() / 2).to_string());
        out.push('\t');
        out.push_str(&(self.right_exons_junctions.len() / 2).to_string());
        out.push('\t');

        let left_exons_str = self
            .left_exons_junctions
            .chunks(2)
            .map(|w| format!("{}-{}", w[0], w[1]))
            .collect::<Vec<String>>()
            .join(",");
        out.push_str(&left_exons_str);
        out.push('\t');
        let right_exons_str = self
            .right_exons_junctions
            .chunks(2)
            .map(|w| format!("{}-{}", w[0], w[1]))
            .collect::<Vec<String>>()
            .join(",");
        out.push_str(&right_exons_str);
        out.push('\t');
        out.push_str(&self.sample_name);
        out.push_str("\n");
        out
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

        let confidence = MergedIsoform::get_confidence_value(
            &evidence_arr,
            dataset_info.get_size(),
            &dataset_info.sample_total_evidence_vec,
        );

        dbg!(confidence);
    }
}
