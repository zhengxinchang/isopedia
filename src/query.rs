use std::os::unix::process;

use log::{debug, info};
use noodles_core::position;

use crate::{
    bptree::BPForest,
    cmd::isoform::AnnIsoCli,
    em::{TxAbundance, MSJC},
    gtf::Transcript,
    isoformarchive::ArchiveCache,
    tmpidx::MergedIsoformOffsetPtr,
};

// gruop tx into classes based on shared splice juncitons
// input, list of Transcript struct
// function, iter_groups() -> QueriesManager
// sites-> tx id
// tx id -> sits

// start with tx_id, create a grouped_tx, then find all splice releated tx to add to the group,
// iterately until no more tx can be added to the group
// then create a new group with remaining tx

pub struct ChromGroupedTxManager {
    chrom: String,
    group_queries: Vec<GroupedTx>,
    sample_size: usize,
    // init_group_len: usize,
    // merged_group_len: usize,
    // _position_to_group_idx: std::collections::HashMap<u64, usize>,
    // offset_to_group_idx: std::collections::HashMap<usize, usize>, // merge gruop if two groups have same offset.
}

impl ChromGroupedTxManager {
    pub fn new(chrom: &str, sample_size: usize) -> Self {
        ChromGroupedTxManager {
            chrom: chrom.to_string(),
            group_queries: Vec::new(),
            sample_size,
            // _position_to_group_idx: std::collections::HashMap::new(),
            // offset_to_group_idx: std::collections::HashMap::new(),
            // init_group_len: 0,
            // merged_group_len: 0,
        }
    }

    pub fn get_mem_size(&self) -> usize {
        let size = std::mem::size_of_val(&self.chrom) + std::mem::size_of_val(&self.group_queries);
        // + std::mem::size_of_val(&self._position_to_group_idx);
        size
    }

    pub fn add_transcript_by_chrom(&mut self, tx_v: &Vec<Transcript>) {
        let mut merged_intervals: Vec<(u64, u64, Vec<&Transcript>)> = Vec::new();

        let curr_tx = tx_v
            .first()
            .unwrap_or_else(|| panic!("No transcript found for chromosome {}", self.chrom));

        let mut curr_mrg_tx = (curr_tx.start, curr_tx.end, vec![curr_tx]);

        for tx in tx_v.iter().skip(1) {
            if tx.start <= curr_mrg_tx.1 {
                // overlap
                curr_mrg_tx.1 = curr_mrg_tx.1.max(tx.end);
                curr_mrg_tx.2.push(tx);
            } else {
                // no overlap, push current and start new
                merged_intervals.push(curr_mrg_tx);
                curr_mrg_tx = (tx.start, tx.end, vec![tx]);
            }
        }

        // push last
        merged_intervals.push(curr_mrg_tx);

        debug!(
            "Chromosome {}: formed {} grouped transcripts from {} transcripts",
            self.chrom,
            merged_intervals.len(),
            tx_v.len()
        );

        // convert each merged group into GroupedTx

        for (_start, _end, txs) in merged_intervals.iter() {
            let grouped_tx = GroupedTx::from_grouped_txs(self.sample_size, txs);
            self.group_queries.push(grouped_tx);
        }
    }

    pub fn query_groups(
        &mut self,
        bpforest: &mut BPForest,
        cli: &AnnIsoCli,
        archive_cache: &mut ArchiveCache,
    ) {
        // query and get results for each group
        let mut processed_cnt = 0usize;
        let total_groups = self.group_queries.len();
        for grouped_tx in self.group_queries.iter_mut() {
            processed_cnt += 1;
            if processed_cnt % 1000 == 0 || processed_cnt == total_groups {
                info!(
                    "Chromosome {}: processing {}/{} groups",
                    self.chrom, processed_cnt, total_groups
                );
            }

            let queries = grouped_tx.get_quieries(&self.chrom);
            // query bpforest
            let mut all_res = bpforest.search2_all_match_nofsm(&queries, cli.flank, cli.lru_size);

            // sort the inner vecs by offset to make sure the same offset are together
            for res_vec in all_res.iter_mut() {
                res_vec.sort_by_key(|x| x.offset);
            }

            // debug!("got {} results", all_res.iter().map(|v| v.len()).sum::<usize>());
            grouped_tx.update_results(&all_res, archive_cache);

            grouped_tx.em(cli);

            grouped_tx.post_cleanup();
        }
    }
}

pub struct GroupedTx {
    pub positions: Vec<u64>, // deduped positions sort by genomic coordinate
    pub tx_abundances: Vec<TxAbundance>,
    pub msjcs: Vec<MSJC>,
}

impl GroupedTx {
    pub fn from_grouped_txs(sample_size: usize, txs: &Vec<&Transcript>) -> Self {
        let mut pos_set = std::collections::BTreeSet::new();

        let mut tx_abundances = Vec::new();
        for (i, tx) in txs.iter().enumerate() {
            for sj in tx.get_splice_junction_pairs().iter() {
                pos_set.insert(sj.0);
                pos_set.insert(sj.1);
            }
            let tx_abd = TxAbundance::new(i, sample_size, tx);
            tx_abundances.push(tx_abd);
        }

        let mut positions: Vec<u64> = pos_set.into_iter().collect();
        positions.sort();

        GroupedTx {
            positions,
            tx_abundances: tx_abundances,
            msjcs: Vec::new(),
        }
    }

    pub fn post_cleanup(&mut self) {
        self.positions.clear();
        self.msjcs.clear();
        self.positions.shrink_to_fit();
        self.msjcs.shrink_to_fit();
    }

    pub fn get_mem_size(&self) -> usize {
        let mut size = std::mem::size_of_val(&self.positions);

        for tx_ab in self.tx_abundances.iter() {
            size += tx_ab.get_mem_size();
        }

        for msjc in self.msjcs.iter() {
            size += msjc.get_mem_size();
        }

        size
    }

    pub fn get_quieries(&self, chrom: &str) -> Vec<(String, u64)> {
        let mut queries: Vec<(String, u64)> = Vec::new();

        for &pos in self.positions.iter() {
            queries.push((chrom.to_string(), pos));
        }

        queries
    }

    pub fn update_results(
        &mut self,
        all_res: &Vec<Vec<MergedIsoformOffsetPtr>>,
        archive_cache: &mut ArchiveCache,
    ) {
        if all_res.len() != self.positions.len() {
            panic!("Number of results does not match number of positions");
        }

        // for each transcript, commpare all results to find 1) fsm 2) misoform offsets
        for (tx_idx, tx_abd) in self.tx_abundances.iter_mut().enumerate() {
            // compare the misoform offsets that match all splice junctions of the transcript and add them into fsm, then find the misoform offsets
            // that match at least one splice junction and add them into msjcs.
            // msjc must deduped, if its already exist, just add the tx_abundance reference.
            // the tx_abundance.sj_pairs is a Vec<(u64,u64)>, so we can use that to find the matching results.
            let mut fsm_misoforms: Vec<MergedIsoformOffsetPtr> = Vec::new();
            let mut msjc_map: std::collections::HashMap<usize, MergedIsoformOffsetPtr> =
                std::collections::HashMap::new();

            for sj_pair in tx_abd.sj_pairs.iter() {
                let sj_start = sj_pair.0;
                let sj_end = sj_pair.1;

                // find the index of the position in self.positions
                let pos_idx = self.positions.iter().position(|&p| p == sj_start).unwrap();
                let res_vec = &all_res[pos_idx];
            }
        }
    }

    pub fn em(&mut self, cli: &AnnIsoCli) {
        // the tx_abundances and msjcs should be already prepared
        // then run em
        let em_iter = cli.em_iter;
        for _ in 0..em_iter {
            // E step
            for msjc in self.msjcs.iter_mut() {
                msjc.e_step(&self.tx_abundances);
            }

            // M step
            for tx_abd in self.tx_abundances.iter_mut() {
                tx_abd.m_step(&self.msjcs);
            }
        }
    }
}
