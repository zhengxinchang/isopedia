use crate::{
    bptree::BPForest,
    em::{TxAbundance, MSJC},
    gtf::TranscriptMini,
    tmpidx::MergedIsoformOffsetPtr,
};

// gruop tx into classes based on shared splice juncitons
// input, list of Transcript struct
// function, iter_groups() -> QueriesManager
pub struct GroupedTxManager {
    chrom: String,
    group_queries: Vec<GroupedTx>,
    init_group_len: usize,
    merged_group_len: usize,
    _position_to_group_idx: std::collections::HashMap<u64, usize>,
    offset_to_group_idx: std::collections::HashMap<usize, usize>, // merge gruop if two groups have same offset.
}

impl GroupedTxManager {
    pub fn new(chrom: &str) -> Self {
        GroupedTxManager {
            chrom: chrom.to_string(),
            group_queries: Vec::new(),
            _position_to_group_idx: std::collections::HashMap::new(),

            offset_to_group_idx: std::collections::HashMap::new(),
            init_group_len: 0,
            merged_group_len: 0,
        }
    }

    pub fn add_transcript(&mut self, tx: &crate::gtf::Transcript) {
        todo!()
        // add transcript to existing group or create a new group. if a transcript shares splice juncitons with existing group, add to that group, else create a new group
    }

    pub fn query_groups(&self, bpforest: &BPForest) {
        todo!()
        // query each group against the bpforest

        // add resuts back to each group, if multiple groups have same offset, merge them

        // dump offests for each group into disk but record the offset and length in memory for later retrieval
    }
}

pub struct GroupedTx {
    pub positions: Vec<u64>, // deduped positions sort by genomic coordinate
    pub query_offsets: Vec<Vec<usize>>,
    pub dedup_offsets: Vec<MergedIsoformOffsetPtr>,
    pub length: usize,
    pub offset_cnt: usize,
    pub tx_minis: Vec<TranscriptMini>,

    pub TxAbundances: Vec<TxAbundance>,
    pub MSJCs: Vec<MSJC>,
}

impl GroupedTx {
    pub fn new() -> Self {
        todo!()
    }

    pub fn get_quieries(&self, chrom: &str) -> Vec<(String, u64)> {
        todo!()
    }

    pub fn update_results(&mut self, all_res: &Vec<Vec<MergedIsoformOffsetPtr>>) {
        todo!()
    }

    pub fn em(&mut self) -> Vec<TxAbundance> {
        todo!()
    }
}
