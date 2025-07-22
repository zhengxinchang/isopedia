use rustc_hash::FxHashMap;

use crate::{
    gene_index::{self, GeneIntervalTree},
    reads::Segment,
    utils::hash_vec,
};

// isoform derived fusion read event.
// This is a per-read level representation of a supplimentary mapped read.
// could be a potential fusion candidate.
pub struct FusionSingleRead<'a> {
    pub fusion_hash: u64, // hash of the fusion read
    pub chr1: String,
    pub chr2: String,
    pub sample_id: u32,
    pub supp_segments: &'a [Segment],
    pub main_splice_junctions_vec: &'a Vec<(u64, u64)>, // main splice junctions of the read>,
}

impl<'a> FusionSingleRead<'a> {
    pub fn new(
        chr1: String,
        chr2: String,
        sample_id: u32,
        supp_segments: &'a [Segment],
        main_splice_junctions_vec: &'a Vec<(u64, u64)>,
    ) -> Self {
        let mut pos_vec = Vec::new();

        for seg in supp_segments.iter() {
            pos_vec.push(seg.start);
            pos_vec.push(seg.end);
        }
        for pos in main_splice_junctions_vec.iter() {
            pos_vec.push(pos.0);
            pos_vec.push(pos.1);
        }

        pos_vec.sort();

        FusionSingleRead {
            fusion_hash: hash_vec(&pos_vec),
            chr1,
            chr2,
            sample_id,
            supp_segments,
            main_splice_junctions_vec,
        }
    }
}

/// aggregated fusion read event.
/// it has a fusion hash that is derived from the main_splice_junctions_vec and supp_splice_junctions_vec.
pub struct FusionAggrReads {
    pub fusion_hash: u64,
    pub chr1: String,
    pub chr2: String,
    pub sample_evidence: FxHashMap<u32, u32>, // sample_id -> evidence count
    pub tootal_evidences: u32,
    pub main_splice_junctions_vec: Vec<(u64, u64)>,
    pub supp_splice_junctions_vec: Vec<(u64, u64)>,
    pub left_matched_gene: String,
    pub right_matched_gene: String,
    pub left_matched_splice_junctions: Vec<u64>,
    pub right_matched_splice_junctions: Vec<u64>,
}

impl FusionAggrReads {
    pub fn init(fusion_single_read: &FusionSingleRead) -> Self {
        let supp_splice_junctions_vec: Vec<(u64, u64)> = fusion_single_read
            .supp_segments
            .iter()
            .map(|seg| (seg.start as u64, seg.end as u64))
            .collect();

        FusionAggrReads {
            fusion_hash: fusion_single_read.fusion_hash,
            chr1: fusion_single_read.chr1.clone(),
            chr2: fusion_single_read.chr2.clone(),
            sample_evidence: FxHashMap::from_iter(vec![(fusion_single_read.sample_id, 1)]),
            tootal_evidences: 1,
            main_splice_junctions_vec: fusion_single_read.main_splice_junctions_vec.clone(),
            supp_splice_junctions_vec: supp_splice_junctions_vec,
            left_matched_gene: String::new(),
            right_matched_gene: String::new(),
            left_matched_splice_junctions: Vec::new(),
            right_matched_splice_junctions: Vec::new(),
        }
    }

    pub fn add(&mut self, other: &FusionSingleRead) {
        self.tootal_evidences += 1;
        *self.sample_evidence.entry(other.sample_id).or_insert(0) += 1
    }

    pub fn match_gene(&mut self, gene_indexing: &GeneIntervalTree, flank: u64) -> bool {
        let mut is_good = true;

        let flat_main_splices = self
            .main_splice_junctions_vec
            .iter()
            .flat_map(|x| vec![x.0, x.1])
            .collect::<Vec<u64>>();

        let results = gene_indexing.match2(self.chr1.as_str(), &flat_main_splices, flank);
        if let Some((chrom, splice_sites)) = results {
            self.left_matched_gene = chrom;
            self.left_matched_splice_junctions = splice_sites;
        } else {
            is_good = false;
        }

        let flat_supp_splices = self
            .supp_splice_junctions_vec
            .iter()
            .flat_map(|x| vec![x.0, x.1])
            .collect::<Vec<u64>>();
        let results = gene_indexing.match2(self.chr2.as_str(), &flat_supp_splices, flank);

        if let Some((chrom, splice_sites)) = results {
            self.right_matched_gene = chrom;
            self.right_matched_splice_junctions = splice_sites;
        } else {
            is_good = false;
        }

        is_good
    }

    pub fn get_string(&self, total_samples: usize) -> String {
        todo!()
    }
}
