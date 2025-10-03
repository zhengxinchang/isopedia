use std::vec;

use crate::{
    dataset_info::DatasetInfo,
    fusion,
    gene_index::{CandidateMatchStatus, GeneInterval, GeneIntervalTree},
    io::{Line, SampleChip},
    output::{FusionBrkPtTableOut, FusionDiscoveryTableOut, GeneralTableOutput},
    reads::Segment,
    utils::{hash_vec, is_overlap},
};
use ahash::HashSet;
use anyhow::Result;
use rustc_hash::FxHashMap;

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
#[derive(Debug, Clone)]
pub struct FusionAggrReads {
    pub fusion_hash: u64,
    pub chr1: String,
    pub chr2: String,
    pub sample_evidence: FxHashMap<u32, u32>, // sample_id -> evidence count
    pub tootal_evidences: u32,
    pub is_ambiguous: bool, // if the fusion is ambiguous, i.e. multiple genes matched
    pub main_splice_junctions_vec: Vec<(u64, u64)>,
    pub supp_splice_junctions_vec: Vec<(u64, u64)>,
    pub left_matched_gene: GeneInterval,
    pub right_matched_gene: GeneInterval,
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
            is_ambiguous: false, // initially not ambiguous
            main_splice_junctions_vec: fusion_single_read.main_splice_junctions_vec.clone(),
            supp_splice_junctions_vec: supp_splice_junctions_vec,
            left_matched_gene: GeneInterval::default(),
            right_matched_gene: GeneInterval::default(),
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

        if let Some((mached_gene1, splice_sites, status)) = results {
            self.left_matched_gene = mached_gene1;
            self.left_matched_splice_junctions = splice_sites;
            if status == CandidateMatchStatus::Ambiguous {
                self.is_ambiguous = true; // if the match is ambiguous, set the flag
            }
        } else {
            is_good = false;
        }

        let flat_supp_splices = self
            .supp_splice_junctions_vec
            .iter()
            .flat_map(|x| vec![x.0, x.1])
            .collect::<Vec<u64>>();
        let results = gene_indexing.match2(self.chr2.as_str(), &flat_supp_splices, flank);

        if let Some((matched_gene2, splice_sites, status)) = results {
            self.right_matched_gene = matched_gene2;
            self.right_matched_splice_junctions = splice_sites;
            if status == CandidateMatchStatus::Ambiguous {
                self.is_ambiguous = true; // if the match is ambiguous, set the flag
            }
        } else {
            is_good = false;
        }

        if self.left_matched_gene.gene_id == self.right_matched_gene.gene_id {
            is_good = false; // both genes are the same, not a fusion
        }

        is_good
    }

    pub fn get_gene_hash(&self) -> (String, String, String) {
        let mut v = vec![
            self.left_matched_gene.gene_name.clone(),
            self.right_matched_gene.gene_name.clone(),
        ];
        v.sort();
        return (
            v.join("_"),
            self.left_matched_gene.gene_name.clone(),
            self.right_matched_gene.gene_name.clone(),
        );
    }
}

#[derive(Debug, Clone)]
pub struct FusionCluster {
    pub gene_hash: String,
    pub gene_combination: HashSet<String>,
    pub gene1: String,
    pub gene2: String,
    pub chr1: String,
    pub chr2: String,
    pub evidence_by_sample: Vec<FxHashMap<u32, u32>>,
    pub total_evidences: u32,
    pub order_a_main_left_most_splice_sites: Vec<u64>,
    pub order_a_main_right_most_splice_sites: Vec<u64>,
    pub order_a_supp_left_most_splice_sites: Vec<u64>,
    pub order_a_supp_right_most_splice_sites: Vec<u64>,
    pub order_b_main_left_most_splice_sites: Vec<u64>,
    pub order_b_main_right_most_splice_sites: Vec<u64>,
    pub order_b_supp_left_most_splice_sites: Vec<u64>,
    pub order_b_supp_right_most_splice_sites: Vec<u64>,
    pub order_a_main_left_mean: f64,
    pub order_a_main_right_mean: f64,
    pub order_a_supp_left_mean: f64,
    pub order_a_supp_right_mean: f64,
    pub order_b_main_left_mean: f64,
    pub order_b_main_right_mean: f64,
    pub order_b_supp_left_mean: f64,
    pub order_b_supp_right_mean: f64,
    pub order_a_main_left_std: f64,
    pub order_a_main_right_std: f64,
    pub order_a_supp_left_std: f64,
    pub order_a_supp_right_std: f64,
    pub order_b_main_left_std: f64,
    pub order_b_main_right_std: f64,
    pub order_b_supp_left_std: f64,
    pub order_b_supp_right_std: f64,
}

impl FusionCluster {
    pub fn new(fusion_aggr: &FusionAggrReads) -> Self {
        let (gene_hash, gene1, gene2) = fusion_aggr.get_gene_hash();

        let evidence_by_sample = vec![fusion_aggr.sample_evidence.clone()];
        let gene_natural_order = format!("{}_{}", gene1, gene2);

        let mut cluster = FusionCluster {
            gene_hash: gene_hash,
            gene_combination: HashSet::from_iter(vec![gene_natural_order]),
            gene1,
            gene2,
            chr1: fusion_aggr.chr1.clone(),
            chr2: fusion_aggr.chr2.clone(),
            evidence_by_sample,
            total_evidences: fusion_aggr.tootal_evidences,
            order_a_main_left_most_splice_sites: Vec::new(),
            order_a_main_right_most_splice_sites: Vec::new(),
            order_a_supp_left_most_splice_sites: Vec::new(),
            order_a_supp_right_most_splice_sites: Vec::new(),
            order_b_main_left_most_splice_sites: Vec::new(),
            order_b_main_right_most_splice_sites: Vec::new(),
            order_b_supp_left_most_splice_sites: Vec::new(),
            order_b_supp_right_most_splice_sites: Vec::new(),

            order_a_main_left_mean: 0.0,
            order_a_main_right_mean: 0.0,
            order_a_supp_left_mean: 0.0,
            order_a_supp_right_mean: 0.0,

            order_b_main_left_mean: 0.0,
            order_b_main_right_mean: 0.0,
            order_b_supp_left_mean: 0.0,
            order_b_supp_right_mean: 0.0,

            order_a_main_left_std: 0.0,
            order_a_main_right_std: 0.0,
            order_a_supp_left_std: 0.0,
            order_a_supp_right_std: 0.0,

            order_b_main_left_std: 0.0,
            order_b_main_right_std: 0.0,
            order_b_supp_left_std: 0.0,
            order_b_supp_right_std: 0.0,
        };

        cluster.order_a_main_left_most_splice_sites.push(
            fusion_aggr
                .left_matched_splice_junctions
                .first()
                .cloned()
                .unwrap_or_default(),
        );
        cluster.order_a_main_right_most_splice_sites.push(
            fusion_aggr
                .left_matched_splice_junctions
                .last()
                .cloned()
                .unwrap_or_default(),
        );
        cluster.order_a_supp_left_most_splice_sites.push(
            fusion_aggr
                .right_matched_splice_junctions
                .first()
                .cloned()
                .unwrap_or_default(),
        );
        cluster.order_a_supp_right_most_splice_sites.push(
            fusion_aggr
                .right_matched_splice_junctions
                .last()
                .cloned()
                .unwrap_or_default(),
        );

        cluster
    }

    pub fn add(&mut self, fusion_aggr: &FusionAggrReads) {
        let (_, gene1, gene2) = fusion_aggr.get_gene_hash();
        let gene_natural_order = format!("{}_{}", gene1, gene2);
        self.total_evidences += fusion_aggr.tootal_evidences;
        self.evidence_by_sample
            .push(fusion_aggr.sample_evidence.clone());
        self.gene_combination.insert(gene_natural_order);

        if (self.gene1 == gene1) && (self.gene2 == gene2) {
            self.order_a_main_left_most_splice_sites.push(
                fusion_aggr
                    .left_matched_splice_junctions
                    .first()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_a_main_right_most_splice_sites.push(
                fusion_aggr
                    .left_matched_splice_junctions
                    .last()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_a_supp_left_most_splice_sites.push(
                fusion_aggr
                    .right_matched_splice_junctions
                    .first()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_a_supp_right_most_splice_sites.push(
                fusion_aggr
                    .right_matched_splice_junctions
                    .last()
                    .cloned()
                    .unwrap_or_default(),
            );
        } else {
            //oerder b
            self.order_b_main_left_most_splice_sites.push(
                fusion_aggr
                    .left_matched_splice_junctions
                    .first()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_b_main_right_most_splice_sites.push(
                fusion_aggr
                    .left_matched_splice_junctions
                    .last()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_b_supp_left_most_splice_sites.push(
                fusion_aggr
                    .right_matched_splice_junctions
                    .first()
                    .cloned()
                    .unwrap_or_default(),
            );
            self.order_b_supp_right_most_splice_sites.push(
                fusion_aggr
                    .right_matched_splice_junctions
                    .last()
                    .cloned()
                    .unwrap_or_default(),
            );
        }
    }

    pub fn evaluate(&mut self) -> bool {
        // calcuate the mean and std of the splice sites

        self.order_a_main_left_mean = self
            .order_a_main_left_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_a_main_left_most_splice_sites.len() as f64;
        self.order_a_main_right_mean = self
            .order_a_main_right_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_a_main_right_most_splice_sites.len() as f64;
        self.order_a_supp_left_mean = self
            .order_a_supp_left_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_a_supp_left_most_splice_sites.len() as f64;
        self.order_a_supp_right_mean = self
            .order_a_supp_right_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_a_supp_right_most_splice_sites.len() as f64;

        let left_std: Vec<f64> = self
            .order_a_main_left_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_a_main_left_mean).powi(2))
            .collect();
        let right_std: Vec<f64> = self
            .order_a_main_right_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_a_main_right_mean).powi(2))
            .collect();
        let supp_left_std: Vec<f64> = self
            .order_a_supp_left_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_a_supp_left_mean).powi(2))
            .collect();
        let supp_right_std: Vec<f64> = self
            .order_a_supp_right_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_a_supp_right_mean).powi(2))
            .collect();

        if !left_std.is_empty() {
            self.order_a_main_left_std =
                (left_std.iter().sum::<f64>() / left_std.len() as f64).sqrt();
        }
        if !right_std.is_empty() {
            self.order_a_main_right_std =
                (right_std.iter().sum::<f64>() / right_std.len() as f64).sqrt();
        }
        if !supp_left_std.is_empty() {
            self.order_a_supp_left_std =
                (supp_left_std.iter().sum::<f64>() / supp_left_std.len() as f64).sqrt();
        }
        if !supp_right_std.is_empty() {
            self.order_a_supp_right_std =
                (supp_right_std.iter().sum::<f64>() / supp_right_std.len() as f64).sqrt();
        }

        // order b
        self.order_b_main_left_mean = self
            .order_b_main_left_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_b_main_left_most_splice_sites.len() as f64;
        self.order_b_main_right_mean = self
            .order_b_main_right_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_b_main_right_most_splice_sites.len() as f64;
        self.order_b_supp_left_mean = self
            .order_b_supp_left_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_b_supp_left_most_splice_sites.len() as f64;
        self.order_b_supp_right_mean = self
            .order_b_supp_right_most_splice_sites
            .iter()
            .copied()
            .map(|x| x as f64)
            .sum::<f64>()
            / self.order_b_supp_right_most_splice_sites.len() as f64;
        let left_std: Vec<f64> = self
            .order_b_main_left_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_b_main_left_mean).powi(2))
            .collect();
        let right_std: Vec<f64> = self
            .order_b_main_right_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_b_main_right_mean).powi(2))
            .collect();
        let supp_left_std: Vec<f64> = self
            .order_b_supp_left_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_b_supp_left_mean).powi(2))
            .collect();
        let supp_right_std: Vec<f64> = self
            .order_b_supp_right_most_splice_sites
            .iter()
            .map(|&x| (x as f64 - self.order_b_supp_right_mean).powi(2))
            .collect();
        if !left_std.is_empty() {
            self.order_b_main_left_std =
                (left_std.iter().sum::<f64>() / left_std.len() as f64).sqrt();
        }
        if !right_std.is_empty() {
            self.order_b_main_right_std =
                (right_std.iter().sum::<f64>() / right_std.len() as f64).sqrt();
        }
        if !supp_left_std.is_empty() {
            self.order_b_supp_left_std =
                (supp_left_std.iter().sum::<f64>() / supp_left_std.len() as f64).sqrt();
        }
        if !supp_right_std.is_empty() {
            self.order_b_supp_right_std =
                (supp_right_std.iter().sum::<f64>() / supp_right_std.len() as f64).sqrt();
        }

        let mut is_good = false;
        // if order_a_main overlaop with order_b_supp, then it is not a good fusion
        let main_a = vec![self.order_a_main_left_mean, self.order_a_main_right_mean];
        let main_b = vec![self.order_b_main_left_mean, self.order_b_main_right_mean];
        let supp_a = vec![self.order_a_supp_left_mean, self.order_a_supp_right_mean];
        let supp_b = vec![self.order_b_supp_left_mean, self.order_b_supp_right_mean];
        if is_overlap(&main_a, &supp_b) && is_overlap(&main_b, &supp_a) {
            is_good = true;
        }

        is_good
    }

    pub fn generate_record_line(
        &self,
        index_info: &DatasetInfo,
        fusiondiscovery_out: &mut FusionDiscoveryTableOut,
    ) -> Result<()> {
        // let mut out_str = String::new();
        let total_samples = index_info.get_size();

        let mut merged_sample_evidence = FxHashMap::default();

        for sample_evidence in self.evidence_by_sample.iter() {
            for (sample_id, count) in sample_evidence.iter() {
                *merged_sample_evidence.entry(*sample_id).or_insert(0) += *count;
            }
        }
        let mut out_line = Line::new();

        let mut sample_count_str = String::new();
        for sample_id in 0..total_samples {
            if let Some(count) = merged_sample_evidence.get(&(sample_id as u32)) {
                sample_count_str.push_str(&format!("{}\t", count));
                let sample_chip = SampleChip::new(None, vec![count.to_string()]);
                out_line.add_sample(sample_chip);
            } else {
                sample_count_str.push_str("0\t");
                out_line.add_sample(SampleChip::new(None, vec!["0".to_string()]));
            }
        }

        out_line.add_field(&self.gene1);
        out_line.add_field(&self.gene2);
        out_line.add_field(&self.total_evidences.to_string());
        out_line.add_field(&total_samples.to_string());
        out_line.add_field(&(self.gene_combination.len() > 1).to_string());
        out_line.add_field(&format!(
            "{}:{}-{}",
            self.chr1, self.order_a_main_left_mean as u32, self.order_a_main_right_mean as u32
        ));
        out_line.add_field(&format!(
            "{}:{}-{}",
            self.chr2, self.order_a_supp_left_mean as u32, self.order_a_supp_right_mean as u32
        ));
        out_line.add_field(&format!(
            "{}:{}-{}",
            self.chr2, self.order_b_main_left_mean as u32, self.order_b_main_right_mean as u32
        ));
        out_line.add_field(&format!(
            "{}:{}-{}",
            self.chr1, self.order_b_supp_left_mean as u32, self.order_b_supp_right_mean as u32
        ));

        fusiondiscovery_out.add_line(&out_line)?;

        Ok(())
    }
}
