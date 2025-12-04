use crate::{
    bptree::BPForest,
    cmd::isoform::AnnIsoCli,
    dataset_info::DatasetInfo,
    global_stats::GlobalStats,
    gtf::Transcript,
    io::{Line, SampleChip},
    isoform::MergedIsoform,
    isoformarchive::ArchiveCache,
    results::TableOutput,
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, intersect_sorted},
};
use ahash::HashMap;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    io::{Read, Seek, Write},
    path::PathBuf,
    sync::Mutex,
};

pub struct ChromGroupedTxManager {
    chrom: String,
    group_queries: Vec<GroupedTx>,
    sample_size: usize,
}

impl ChromGroupedTxManager {
    pub fn new(chrom: &str, sample_size: usize) -> Self {
        ChromGroupedTxManager {
            chrom: chrom.to_string(),
            group_queries: Vec::new(),
            sample_size,
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

    pub fn process_tx_groups(
        &mut self,
        bpforest: &mut BPForest,
        cli: &AnnIsoCli,
        archive_cache: &mut ArchiveCache,
        dbinfo: &DatasetInfo,
        tmp_tx_manager: &mut Mutex<TmpOutputManager>,
        global_stats: &mut GlobalStats,
    ) {
        let total_groups = self.group_queries.len();
        let processed_cnt = Mutex::new(0usize);
        let global_stats_mutex = Mutex::new(global_stats);

        // Make archive_cache and bpforest thread-safe
        let archive_cache_mutex = Mutex::new(archive_cache);
        let bpforest_mutex = Mutex::new(bpforest);

        self.group_queries.par_iter_mut().for_each(|grouped_tx| {
            let queries = grouped_tx.get_quieries(&self.chrom);

            // Lock resources when needed
            let mut all_res = {
                let mut bp = bpforest_mutex.lock().unwrap();
                bp.search2_all_match_only(&queries, cli.flank, cli.cached_nodes)
            };

            for res_vec in all_res.iter_mut() {
                res_vec.sort_by_key(|x| x.offset);
            }

            {
                let mut cache = archive_cache_mutex.lock().unwrap();
                grouped_tx.update_results(&all_res, &mut *cache, dbinfo, cli);
            }

            grouped_tx.prepare_em();
            grouped_tx.em(cli);
            grouped_tx.post_cleanup(&cli);

            // Update global stats,
            {
                let mut stats = global_stats_mutex.lock().unwrap();
                for txbd in grouped_tx.tx_abundances.iter() {
                    stats.update_fsm_total(&txbd.fsm_abundance);
                    stats.update_em_total(&txbd.abundance_cur);
                }
            }

            // Update counter and log
            {
                let mut cnt = processed_cnt.lock().unwrap();
                *cnt += 1;
                if *cnt % 10 == 0 || *cnt == total_groups {
                    info!(
                        "Chromosome {}: processing {}/{} groups",
                        self.chrom, *cnt, total_groups
                    );
                }
            }

            // write to table output

            {
                let mut tmp_manager = tmp_tx_manager.lock().unwrap();
                tmp_manager.dump_grouped_tx(grouped_tx);
            }
        });
    }
}

pub struct GroupedTx {
    pub positions: Vec<u64>, // deduped positions sort by genomic coordinate
    pub position_tx_abd_map: std::collections::HashMap<u64, usize>, // position -> tx_abundance ids
    pub tx_abundances: Vec<TxAbundance>,
    pub msjcs: HashMap<u64, MSJC>,
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

        let mut position_tx_abd_map = std::collections::HashMap::new();

        for i in 0..positions.len() {
            position_tx_abd_map.insert(positions[i], i);
        }

        GroupedTx {
            positions,
            position_tx_abd_map: position_tx_abd_map,
            tx_abundances: tx_abundances,
            msjcs: HashMap::default(),
        }
    }

    pub fn post_cleanup(&mut self, _cli: &AnnIsoCli) {
        self.positions.clear();
        self.msjcs.clear();
        self.positions.shrink_to_fit();
        self.msjcs.shrink_to_fit();
        for tx_abd in self.tx_abundances.iter_mut() {
            tx_abd.msjc_ids.clear();
            tx_abd.msjc_ids.shrink_to_fit();
            tx_abd.sj_pairs.clear();
            tx_abd.sj_pairs.shrink_to_fit();
            tx_abd.abundance_prev.clear();
            tx_abd.abundance_prev.shrink_to_fit();
        }
    }

    pub fn get_mem_size(&self) -> usize {
        let mut size = std::mem::size_of_val(&self.positions);

        for tx_ab in self.tx_abundances.iter() {
            size += tx_ab.get_mem_size();
        }

        for msjc in self.msjcs.values() {
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
        dbinfo: &DatasetInfo,
        cli: &AnnIsoCli,
    ) {
        if all_res.len() != self.positions.len() {
            panic!("Number of results does not match number of positions");
        }

        // dbg!(&all_res[0]);

        let mut msjc_map_global: HashMap<u64, MSJC> = HashMap::default();

        if cli.verbose {
            info!(
                "Updating results for grouped tx with {} transcripts and {} positions",
                self.tx_abundances.len(),
                self.positions.len()
            );
        }

        // for each transcript, commpare all results to find 1) fsm 2) misoform offsets
        for (tx_idx, tx_abd) in self.tx_abundances.iter_mut().enumerate() {
            // compare the misoform offsets that match all splice junctions of the transcript and add them into fsm, then find the misoform offsets
            // that match at least one splice junction and add them into msjcs.
            // msjc must deduped, if its already exist, just add the tx_abundance reference.
            // the tx_abundance.sj_pairs is a Vec<(u64,u64)>, so we can use that to find the matching results.
            let mut fsm_misoforms: Vec<&Vec<MergedIsoformOffsetPtr>> = Vec::new();

            // at least mapping with one sj to be considered as msjc
            let mut partial_msjc_map: HashMap<u64, &MergedIsoformOffsetPtr> = HashMap::default();

            for sj_pair in tx_abd.sj_pairs.iter() {
                let sj_start_idx = self
                    .position_tx_abd_map
                    .get(&sj_pair.0)
                    .unwrap_or_else(|| panic!("Can not get position index for {}", sj_pair.0));
                let sj_end_idx = self
                    .position_tx_abd_map
                    .get(&sj_pair.1)
                    .unwrap_or_else(|| panic!("Can not get position index for {}", sj_pair.1));

                fsm_misoforms.push(&all_res[*sj_start_idx]);
                fsm_misoforms.push(&all_res[*sj_end_idx]);
                quick_find_partial_msjc_exclude_terminals(
                    &all_res[*sj_start_idx],
                    &all_res[*sj_end_idx],
                    &mut partial_msjc_map,
                );
            }

            let fsm_msjc_offsets = find_fsm(fsm_misoforms);

            // remove fsm offsets from msjc_map
            for fsm in fsm_msjc_offsets.iter() {
                partial_msjc_map.remove(&fsm.offset);
            }

            if cli.verbose {
                info!(
                    "Transcript {}: found {} FSM misoforms and {} MSJC misoforms",
                    tx_idx,
                    fsm_msjc_offsets.len(),
                    partial_msjc_map.len()
                );
            }

            // add fsm misoforms
            for fsm in fsm_msjc_offsets.into_iter() {
                // info!("adding fsm misoform offset {:?}", fsm);
                let fsm_msjc_rec = archive_cache.read_bytes(&fsm);
                tx_abd.add_fsm_misoform(&fsm_msjc_rec);
            }

            // add msjcs

            for (offset, msjc_ptr) in partial_msjc_map.into_iter() {
                if msjc_map_global.contains_key(&offset) {
                    let existing_msjc = msjc_map_global.get_mut(&offset).unwrap();
                    // let mut existing_msjc = archive_cache.read_bytes(&existing_msjc_ptr);
                    let local_id = existing_msjc.add_txabundance(tx_abd);
                    tx_abd.msjc_ids.push((offset, local_id));
                } else {
                    let mut msjc = MSJC::new(
                        offset,
                        dbinfo.get_size(),
                        &archive_cache.read_bytes(&msjc_ptr),
                    );
                    let local_id = msjc.add_txabundance(tx_abd);
                    tx_abd.msjc_ids.push((offset, local_id));
                    // self.msjcs.push(msjc);
                    msjc_map_global.insert(offset, msjc);
                }
            }
        }

        self.msjcs = msjc_map_global;
    }

    pub fn prepare_em(&mut self) {
        for msjc in self.msjcs.values_mut() {
            msjc.prepare_em(self.tx_abundances.len());
        }

        // clean up the position_tx_abd_map to save memory
        self.position_tx_abd_map.clear();
        self.position_tx_abd_map.shrink_to_fit();
    }

    pub fn em(&mut self, cli: &AnnIsoCli) {
        // the tx_abundances and msjcs should be already prepared
        // then run em

        if cli.verbose {
            info!(
                "Starting EM for grouped tx with {} transcripts and {} msjcs, memory usage: {:.2} MB",
                self.tx_abundances.len(),
                self.msjcs.len(),
                self.get_mem_size() as f64 / 1024.0 / 1024.0
            );
        }

        let start_time = std::time::Instant::now();

        for i in 0..cli.em_iter {
            if cli.verbose {
                info!("====== Starting EM iteration {} ======", i + 1);
            }

            // E step
            for msjc in self.msjcs.values_mut() {
                msjc.e_step(&self.tx_abundances);
            }

            // M step
            for tx_abd in self.tx_abundances.iter_mut() {
                tx_abd.m_step(&self.msjcs, cli);
            }

            // check convergence
            let mut all_converged = true;
            for tx_abd in self.tx_abundances.iter_mut() {
                let is_converged = tx_abd.check_convergence(cli.em_converge, cli);
                if !is_converged {
                    all_converged = false;
                }
            }
            if all_converged {
                if cli.verbose {
                    info!("All transcripts converged at iteration {}", i + 1);
                }
                break;
            }
        }

        let duration = start_time.elapsed();

        if cli.verbose {
            info!(
                "EM finished for grouped tx with {} transcripts in {:?}",
                self.tx_abundances.len(),
                duration
            );
        }

        // dbg!(&self.tx_abundances);
    }
}

pub fn quick_find_partial_msjc_exclude_terminals<'a>(
    a: &'a Vec<MergedIsoformOffsetPtr>,
    b: &'a Vec<MergedIsoformOffsetPtr>,
    collection: &mut HashMap<u64, &'a MergedIsoformOffsetPtr>,
) {
    // let mut res: Vec<MergedIsoformOffsetPtr> = Vec::new();

    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        if a[i].offset == b[j].offset {
            // make sure n_splice_site are more than 0, if its 0 ,then its the read terminal
            if a[i].n_splice_sites == 0 {
                i += 1;
                j += 1;
                continue;
            }

            collection.insert(a[i].offset, &a[i]);

            i += 1;
            j += 1;
        } else if a[i].offset < b[j].offset {
            i += 1;
        } else {
            j += 1;
        }
    }
}

pub fn find_fsm(vecs: Vec<&Vec<MergedIsoformOffsetPtr>>) -> Vec<MergedIsoformOffsetPtr> {
    if vecs.is_empty() {
        return vec![];
    }

    let n_splice_sites = vecs.len(); // each splice junction has two positions
    let mut result = vecs[0].clone();

    // dbg!(n_splice_sites);
    // dbg!(&result);
    // result.dedup();

    for v in vecs.iter().skip(1) {
        result = intersect_sorted(&result, v);

        // dbg!(result.len());

        if result.is_empty() {
            break;
        }
    }

    // filter out the MergedIsoformOffsetPtr dose not exeactly match n_splice
    result.retain(|x| x.n_splice_sites as usize == n_splice_sites);

    result
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TxAbundance {
    pub id: usize, // idx in the grouped_tx
    pub orig_idx: u32,
    pub orig_tx_id: String,
    pub orig_gene_id: String,
    pub orig_tx_len: u32,
    pub orig_n_exon: u32,
    pub orig_chrom: String,
    pub orig_start: u64,
    pub orig_end: u64,
    pub orig_attrs: String,
    pub sj_pairs: Vec<(u64, u64)>, // splice junction positions
    pub fsm_abundance: Vec<f32>,
    pub abundance_cur: Vec<f32>, // sample size length
    pub abundance_prev: Vec<f32>,
    sample_size: usize,
    inv_pt: f32, // scaling factor based solely on the complexity of transcript, same role as effective length,
    // is_ok: Vec<bool>,
    pub msjc_ids: Vec<(u64, usize)>, // (msjc id, index in msjc)
}

impl TxAbundance {
    pub fn new(txid: usize, sample_size: usize, tx: &Transcript) -> TxAbundance {
        // calclulate pt
        let nsj = tx.splice_junc.len() + 10; // J + 1
        TxAbundance {
            id: txid,
            orig_idx: tx.origin_idx,
            orig_tx_id: tx.tx_id.clone(),
            orig_gene_id: tx.gene_id.clone(),
            orig_tx_len: tx.get_transcript_seq_length() as u32,
            orig_n_exon: tx.get_exon_count() as u32,
            orig_chrom: tx.chrom.clone(),
            orig_start: tx.start,
            orig_end: tx.end,
            orig_attrs: tx.get_attributes(),
            sj_pairs: tx.get_splice_junction_pairs(),
            fsm_abundance: vec![0.0; sample_size],
            abundance_cur: vec![1.0; sample_size],
            abundance_prev: vec![1.0; sample_size],
            sample_size,
            inv_pt: 1.0 / nsj as f32,
            // is_ok: vec![true; sample_size],
            msjc_ids: Vec::new(),
        }
    }

    pub fn add_fsm_misoform(&mut self, misoform: &MergedIsoform) {
        for sid in 0..self.sample_size {
            self.fsm_abundance[sid] += misoform.sample_evidence_arr[sid] as f32;
        }
    }

    pub fn get_mem_size(&self) -> usize {
        let size = std::mem::size_of_val(&self.id)
            + std::mem::size_of_val(&0f32) * self.fsm_abundance.len()
            + std::mem::size_of_val(&0f32) * self.abundance_cur.len()
            + std::mem::size_of_val(&0f32) * self.abundance_prev.len()
            + std::mem::size_of_val(&self.sample_size)
            + std::mem::size_of_val(&self.inv_pt)
            // + std::mem::size_of_val(&self.is_ok)
            + std::mem::size_of_val(&self.msjc_ids);
        size
    }

    pub fn m_step(&mut self, msjc_vec: &HashMap<u64, MSJC>, cli: &AnnIsoCli) {
        if cli.verbose {
            info!(
                "Updating TxAbundance id {},tx_id {:?},msjc length: {}",
                self.id,
                self.orig_tx_id,
                self.msjc_ids.len()
            );
        }

        let sample_size = self.sample_size;

        for sid in 0..sample_size {
            let mut total_abd = 0.0f32;

            for (msjc_id, tx_local_id) in self.msjc_ids.iter() {
                // HashMap 仍用安全访问，真正热点在 Vec 上
                let msjc = msjc_vec
                    .get(msjc_id)
                    .unwrap_or_else(|| panic!("MSJC id {} not found", msjc_id));

                let idx = *tx_local_id * sample_size + sid;

                unsafe {
                    // cov_vec[sid]
                    let cov = *msjc.cov_vec.get_unchecked(sid);
                    // gamma_vec[idx]
                    let gamma = *msjc.gamma_vec.get_unchecked(idx);
                    total_abd += cov * gamma;
                }
            }

            unsafe {
                // abundance_prev[sid] = abundance_cur[sid];
                let prev_ptr = self.abundance_prev.get_unchecked_mut(sid);
                let cur_ptr = self.abundance_cur.get_unchecked_mut(sid);
                *prev_ptr = *cur_ptr;
                *cur_ptr = total_abd;
            }

            if cli.verbose {
                unsafe {
                    let prev = *self.abundance_prev.get_unchecked(sid);
                    let cur = *self.abundance_cur.get_unchecked(sid);
                    let fsm = *self.fsm_abundance.get_unchecked(sid);

                    info!(
                        "TxAbundance id {}, tx_id {:?}, sample {} prev {:.6} cur {:.6} fsm {:.6}",
                        self.id, self.orig_tx_id, sid, prev, cur, fsm
                    );
                }
            }
        }
    }

    pub fn check_convergence(&mut self, tol: f32, cli: &AnnIsoCli) -> bool {
        // check convergence
        // if relative change for all samples are below tol, then converged
        let mut is_converged = true;
        for sid in 0..self.sample_size {
            let prev = self.abundance_prev[sid];
            let cur = self.abundance_cur[sid];
            let diff = (cur - prev).abs();
            let denom = if prev.abs() < 1e-6 { 1.0 } else { prev.abs() };
            let rel_diff = diff / denom;
            if cli.verbose {
                info!(
                    "TxAbundance id {}, tx_id {:?}, sample {} prev {:.6} cur {:.6},diff {:.6} rel_diff {:.6}",
                    self.id, self.orig_tx_id, sid, prev, cur, diff, rel_diff
                );
            }

            if rel_diff > tol {
                is_converged = false;
            }
        }
        is_converged
    }

    pub fn encode(&self) -> Vec<u8> {
        bincode::serialize(self).expect("Can not serialize the grouped_tx")
    }

    pub fn decode(data: &[u8]) -> Self {
        bincode::deserialize(data).expect("Can not deserialize the grouped_tx")
    }

    pub fn get_positive_samples(&self, min_read: f32) -> usize {
        let mut count = 0;
        for (abd1, abd2) in self.abundance_cur.iter().zip(&self.fsm_abundance) {
            let total = abd1 + abd2;
            if total >= min_read {
                count += 1;
            }
        }
        count
    }

    pub fn to_output_line(
        &self,
        global_stats: &GlobalStats,
        dbinfo: &DatasetInfo,
        cli: &AnnIsoCli,
    ) -> Line {
        let mut out_line = Line::new();

        out_line.add_field(&self.orig_chrom);
        out_line.add_field(&self.orig_start.to_string());
        out_line.add_field(&self.orig_end.to_string());
        out_line.add_field(&self.orig_tx_len.to_string());
        out_line.add_field(&self.orig_n_exon.to_string());
        out_line.add_field(&self.orig_tx_id.to_string());
        out_line.add_field(&self.orig_gene_id.to_string());

        // confidence scores

        let u32_arr = self
            .fsm_abundance
            .iter()
            .zip(&self.abundance_cur)
            .map(|(&fsm, em)| fsm + em)
            .collect::<Vec<f32>>();

        let confscore =
            utils::calc_confidence_f32(&u32_arr, dbinfo.get_size(), &global_stats.fsm_em_total);

        out_line.add_field(&confscore.to_string());

        let positive_samples = self.get_positive_samples(cli.min_read as f32);

        if positive_samples != 0 {
            out_line.add_field("yes");
        } else {
            out_line.add_field("no");
        }

        out_line.add_field(&format!("{}/{}", positive_samples, dbinfo.get_size()));

        out_line.add_field(&self.orig_attrs);

        // process each sample

        for sid in 0..dbinfo.get_size() {
            let fsm = self.fsm_abundance[sid];
            let em = self.abundance_cur[sid];
            let total = fsm + em;
            let total_cov = global_stats.fsm_em_total[sid];

            let fsm_cpm = utils::calc_cpm_f32(&fsm, &total_cov);
            let em_cpm = utils::calc_cpm_f32(&em, &total_cov);
            let total_cpm = utils::calc_cpm_f32(&total, &total_cov);

            let sample = SampleChip::new(
                Some(dbinfo.get_sample_names()[sid].clone()),
                vec![format!(
                    "{}:{}:{}:{}:{}:{}",
                    total_cpm, total, fsm_cpm, fsm, em_cpm, em
                )],
            );
            out_line.add_sample(sample);
        }

        out_line
    }
}

pub struct MSJC {
    id: u64,
    sample_size: usize,
    cov_vec: Vec<f32>,
    txids: Vec<usize>,
    gamma_vec: Vec<f32>, // length is same as txids
}

impl MSJC {
    pub fn new(id: u64, sample_size: usize, misoform: &MergedIsoform) -> MSJC {
        let cov = misoform
            .sample_evidence_arr
            .iter()
            .map(|&x| x as f32)
            .collect();

        // info!("Creating MSJC id {} sample_size {} cov_vec {:?}", id, sample_size, &cov);
        MSJC {
            id,
            sample_size,
            cov_vec: cov,
            txids: Vec::new(),
            gamma_vec: Vec::new(),
        }
    }

    pub fn add_txabundance(&mut self, txabd: &TxAbundance) -> usize {
        self.txids.push(txabd.id);
        self.txids.len() - 1
    }

    pub fn prepare_em(&mut self, tx_len: usize) {
        self.gamma_vec = vec![0.0; self.sample_size * tx_len];
    }

    pub fn get_mem_size(&self) -> usize {
        let size = std::mem::size_of_val(&self.id)
            + std::mem::size_of_val(&self.sample_size)
            + std::mem::size_of_val(&self.cov_vec)
            + std::mem::size_of_val(&self.txids)
            + std::mem::size_of_val(&self.gamma_vec);
        size
    }

    pub fn e_step2(&mut self, txabds: &Vec<TxAbundance>) {
        let mut total = vec![0.0f32; self.sample_size];

        // pass 1: total[s]
        for &txid in self.txids.iter() {
            let txabd = &txabds[txid];
            for sid in 0..self.sample_size {
                total[sid] += txabd.abundance_cur[sid] * txabd.inv_pt;
            }
        }

        for (tid, &txid) in self.txids.iter().enumerate() {
            // dbg!("inv_pt of txid {}: {}", txid, txabds[txid].inv_pt);
            let txabd = &txabds[txid];
            let base = tid * self.sample_size;
            for sid in 0..self.sample_size {
                let numer = txabd.abundance_cur[sid] * txabd.inv_pt;
                // let denom = total[sid];
                let idx = base + sid;
                if total[sid] > 0.0 {
                    self.gamma_vec[idx] = numer / total[sid];
                } else {
                    self.gamma_vec[idx] = 0.0;
                }
            }
        }
        // info!("updated gamma_vec after E step: {:?}", self.gamma_vec);
    }

    pub fn e_step(&mut self, txabds: &Vec<TxAbundance>) {
        let sample_size = self.sample_size;
        let txabds_slice = txabds.as_slice();

        // total[s] = Σ_t α_t[s] * inv_pt_t
        let mut total = vec![0.0f32; sample_size];

        // pass 1: 计算 total
        for &txid in self.txids.iter() {
            unsafe {
                let txabd = txabds_slice.get_unchecked(txid);
                let inv_pt = txabd.inv_pt;
                let abd_slice = txabd.abundance_cur.as_slice();

                for sid in 0..sample_size {
                    let t = total.get_unchecked_mut(sid);
                    *t += *abd_slice.get_unchecked(sid) * inv_pt;
                }
            }
        }

        // pass 2: 计算 gamma
        for (tid, &txid) in self.txids.iter().enumerate() {
            unsafe {
                let txabd = txabds_slice.get_unchecked(txid);
                let inv_pt = txabd.inv_pt;
                let abd_slice = txabd.abundance_cur.as_slice();

                let base = tid * sample_size;

                for sid in 0..sample_size {
                    let total_s = *total.get_unchecked(sid);
                    let idx = base + sid;
                    let gamma_slot = self.gamma_vec.get_unchecked_mut(idx);

                    if total_s > 0.0 {
                        *gamma_slot = *abd_slice.get_unchecked(sid) * inv_pt / total_s;
                    } else {
                        *gamma_slot = 0.0;
                    }
                }
            }
        }
    }
}

type orig_idx = u32;
type offset = u64;
type length = u64;

/// A temporary output manager to dump grouped_tx results before collect all informations
pub struct TmpOutputManager {
    pub curr_offset: u64,
    pub curr_processed_idx: usize,
    pub offsets: Vec<(orig_idx, offset, length)>, // (orig_tx_idx, offset, length)
    pub file_path: PathBuf,
    pub file_handle: std::fs::File,
}

impl TmpOutputManager {
    pub fn new(file_path: &PathBuf) -> Self {
        let file_handle = std::fs::File::create(file_path).unwrap_or_else(|_| {
            panic!(
                "Could not create temporary output file {}",
                file_path.display()
            )
        });
        TmpOutputManager {
            curr_offset: 0,
            curr_processed_idx: 0,
            offsets: Vec::new(),
            file_path: file_path.clone(),
            file_handle,
        }
    }

    pub fn dump_grouped_tx(&mut self, grouped_tx: &GroupedTx) {
        for tx_abd in grouped_tx.tx_abundances.iter() {
            let bytes = tx_abd.encode();
            let length = bytes.len() as u64;
            self.file_handle.write_all(&bytes).unwrap_or_else(|_| {
                panic!(
                    "Could not write to temporary output file {}",
                    self.file_path.display()
                )
            });
            self.offsets
                .push((tx_abd.orig_idx, self.curr_offset, length));
            self.curr_offset += length;
        }
    }

    pub fn finish(&mut self) {
        self.offsets.sort_by_key(|x| x.0);
        self.file_handle.sync_all().unwrap_or_else(|_| {
            panic!(
                "Could not sync temporary output file {}",
                self.file_path.display()
            )
        });
        self.file_handle = std::fs::File::open(&self.file_path).unwrap_or_else(|_| {
            panic!(
                "Could not open temporary output file {}",
                self.file_path.display()
            )
        });
        self.curr_processed_idx = 0;
    }

    pub fn get_next_txabundance(&mut self) -> Option<TxAbundance> {
        if self.curr_processed_idx >= self.offsets.len() {
            return None;
        }

        let (_orig_idx, offset, length) = self.offsets[self.curr_processed_idx];
        self.file_handle
            .seek(std::io::SeekFrom::Start(offset))
            .unwrap_or_else(|_| {
                panic!(
                    "Could not seek in temporary output file {}",
                    self.file_path.display()
                )
            });
        let mut buffer = vec![0u8; length as usize];
        self.file_handle
            .read_exact(&mut buffer)
            .unwrap_or_else(|e| {
                panic!(
                    "Could not read from temporary output file {}, with offset {}, length {}: {}",
                    self.file_path.display(),
                    offset,
                    length,
                    e
                )
            });
        let tx_abd = TxAbundance::decode(&buffer);
        self.curr_processed_idx += 1;
        Some(tx_abd)
    }
}

impl Iterator for TmpOutputManager {
    type Item = TxAbundance;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_txabundance()
    }
}
