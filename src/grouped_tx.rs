use crate::{
    bptree::BPForest,
    cmd::isoform::AnnIsoCli,
    dataset_info::DatasetInfo,
    global_stats::GlobalStats,
    gtf::Transcript,
    isoform::MergedIsoform,
    isoformarchive::ArchiveCache,
    results::TableOutput,
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, intersect_sorted, GetMemSize},
};
use ahash::HashMap;
use anyhow::Result;
use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    path::PathBuf,
};
pub struct ChromGroupedTxManager {
    pub chrom: String,
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

    pub fn clear(&mut self) {
        self.group_queries.clear();
        self.group_queries.shrink_to_fit();
    }

    pub fn add_transcript_by_chrom(&mut self, tx_v: &Vec<Transcript>, cli: &AnnIsoCli) {
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

        for (idx, (_start, _end, txs)) in merged_intervals.iter().enumerate() {
            let mut grouped_tx =
                GroupedTx::from_grouped_txs(self.sample_size, txs, idx as u32, cli);
            grouped_tx.id = idx as u32;
            self.group_queries.push(grouped_tx);
        }
    }

    pub fn process_tx_groups(
        &mut self,
        bpforest: &mut BPForest,
        cli: &AnnIsoCli,
        archive_cache: &mut ArchiveCache,
        dbinfo: &DatasetInfo,
        tmp_tx_manager: &mut TmpOutputManager,
        global_stats: &mut GlobalStats,
    ) {
        let total_groups = self.group_queries.len();
        let mut processed_cnt = 0;

        let cli_arc = std::sync::Arc::new(cli.clone());

        for (chunk_idx, chunk) in self.group_queries.chunks_mut(cli.em_chunk_size).enumerate() {
            let start = std::time::Instant::now();

            for grouped_tx in chunk.iter_mut() {
                processed_cnt += 1;
                if processed_cnt % 100 == 0 || processed_cnt == total_groups {
                    info!(
                        "Chromosome {}: processing {}/{} groups",
                        self.chrom, processed_cnt, total_groups
                    );
                }

                let queries = grouped_tx.get_quieries(&self.chrom);

                let mut all_res = bpforest.all_candidate_match_search_flank(
                    &queries,
                    cli.flank,
                    cli.cached_nodes,
                );

                for res_vec in all_res.iter_mut() {
                    res_vec.sort_by_key(|x| x.offset);
                }

                let all_res_mono_exonic = bpforest.search_mono_exons(
                    &self.chrom,
                    &grouped_tx.get_quieries_mono_exon(),
                    cli.flank,
                    cli.cached_nodes,
                );

                // info!("{:?}",all_res_mono_exonic);

                grouped_tx.update_results(
                    &all_res,
                    &all_res_mono_exonic,
                    archive_cache,
                    dbinfo,
                    cli,
                );
                drop(all_res);
                drop(all_res_mono_exonic);

                grouped_tx.prepare_em(cli);
            }

            chunk.par_iter_mut().for_each(|grouped_tx| {
                grouped_tx.em(cli_arc.as_ref());

                grouped_tx.post_cleanup(cli);
            });

            for grouped_tx in chunk.iter_mut() {
                for txbd in grouped_tx.tx_abundances.iter() {
                    global_stats.update_fsm_total(&txbd.fsm_abundance);
                    global_stats.update_em_total(&txbd.abundance_cur);
                }
                tmp_tx_manager.dump_grouped_tx(grouped_tx);
            }

            let duration = start.elapsed();

            debug!(
                "EM for chunk {} of size {} took {:?}",
                chunk_idx,
                chunk.len(),
                duration
            );
        }
    }
}

#[derive(Clone)]
pub struct GroupedTx {
    pub id: u32,
    pub positions: Vec<u64>, // deduped positions sort by genomic coordinate
    pub positions_mono_exonic: Vec<(u64, u64)>, //  positions for mono exonic transcripts
    pub position_tx_abd_map: FxHashMap<u64, usize>, // position -> tx_abundance ids
    pub tx_abundances: Vec<TxAbundance>,
    pub msjcs: Vec<MSJC>,
}
impl GroupedTx {
    pub fn from_grouped_txs(
        sample_size: usize,
        txs: &Vec<&Transcript>,
        id: u32,
        cli: &AnnIsoCli,
    ) -> Self {
        let mut pos_set = std::collections::BTreeSet::new();
        let mut pos_set_mono_exonic = Vec::new();

        let mut tx_abundances = Vec::new();
        for (i, tx) in txs.iter().enumerate() {
            if tx.is_mono_exonic {
                // add logic for mono exonic transcripts
                pos_set_mono_exonic.push((tx.start, tx.end));
            } else {
                for sj in tx.get_splice_junction_pairs().iter() {
                    pos_set.insert(sj.0);
                    pos_set.insert(sj.1);
                }
            }
            // no difference in creating TxAbundance for mono exonic transcripts or not
            let tx_abd = TxAbundance::new(i, sample_size, tx, cli);
            tx_abundances.push(tx_abd);
        }

        let mut positions: Vec<u64> = pos_set.into_iter().collect();
        positions.sort();

        let mut position_tx_abd_map = FxHashMap::default();

        for i in 0..positions.len() {
            position_tx_abd_map.insert(positions[i], i);
        }

        GroupedTx {
            id,
            positions,
            positions_mono_exonic: pos_set_mono_exonic,
            position_tx_abd_map: position_tx_abd_map,
            tx_abundances: tx_abundances,
            msjcs: Vec::new(),
        }
    }

    pub fn post_cleanup(&mut self, _cli: &AnnIsoCli) {
        for msjc in self.msjcs.iter_mut() {
            msjc.gamma_data.clear();
            msjc.gamma_data.shrink_to_fit();
            msjc.cov_vec.clear();
            msjc.cov_vec.shrink_to_fit();
            msjc.txids.clear();
            msjc.txids.shrink_to_fit();
            msjc.totals_buffer.clear();
            msjc.totals_buffer.shrink_to_fit();
            msjc.nonzero_sample_indices.clear();
            msjc.nonzero_sample_indices.shrink_to_fit();
        }

        self.positions.clear();
        self.positions.shrink_to_fit();

        self.position_tx_abd_map.clear();
        self.position_tx_abd_map.shrink_to_fit();

        self.msjcs.clear();
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

    pub fn get_quieries(&self, chrom: &str) -> Vec<(String, u64)> {
        let mut queries: Vec<(String, u64)> = Vec::new();

        for &pos in self.positions.iter() {
            queries.push((chrom.to_string(), pos));
        }

        queries
    }

    pub fn get_quieries_mono_exon(&self) -> Vec<(u64, u64)> {
        self.positions_mono_exonic.clone()
    }

    pub fn update_results(
        &mut self,
        all_res: &Vec<Vec<MergedIsoformOffsetPtr>>,
        all_res_mono_exonic: &Vec<Vec<MergedIsoformOffsetPtr>>,
        archive_cache: &mut ArchiveCache,
        dbinfo: &DatasetInfo,
        cli: &AnnIsoCli,
    ) {
        if all_res.len() != self.positions.len() {
            panic!("Number of results does not match number of positions");
        }

        if all_res_mono_exonic.len() != self.positions_mono_exonic.len() {
            panic!("Number of mono exonic results does not match number of mono exonic positions");
        }

        let mut msjc_map_global: HashMap<u64, usize> = HashMap::default();

        if cli.verbose {
            info!(
                "Updating results for grouped tx with {} transcripts and {} positions",
                self.tx_abundances.len(),
                self.positions.len()
            );
        }

        // for each transcript, commpare all results to find 1) fsm 2) misoform offsets
        let mut local_mono_exonic_idx = 0;
        for (tx_idx, tx_abd) in self.tx_abundances.iter_mut().enumerate() {
            // mono exonic need special query
            if tx_abd.is_mono_exonic {
                let res_vec = &all_res_mono_exonic[local_mono_exonic_idx];

                if res_vec.is_empty() {
                    local_mono_exonic_idx += 1;
                    continue;
                }

                for misoform_ptr in res_vec.iter() {
                    let misoform = archive_cache.read_bytes(&misoform_ptr);

                    if misoform.is_mono_exonic() {
                        // only consider mono exonic misoforms records here

                        let msjc_idx = if let Some(&idx) = msjc_map_global.get(&misoform_ptr.offset)
                        {
                            idx
                        } else {
                            let msjc = MSJC::new(
                                misoform_ptr.offset,
                                dbinfo.get_size(),
                                &archive_cache.read_bytes(&misoform_ptr),
                            );
                            self.msjcs.push(msjc);
                            let new_idx = self.msjcs.len() - 1;
                            msjc_map_global.insert(misoform_ptr.offset, new_idx);
                            new_idx
                        };

                        assert!(msjc_idx < self.msjcs.len());

                        // check if fsm

                        let msjc = &self.msjcs[msjc_idx];

                        if msjc.check_mono_exon_fsm(&tx_abd, cli) {
                            tx_abd.add_fsm_misoform(&misoform);
                        } else {
                            if msjc.check_msjc_belong_to_tx(&tx_abd, cli) {
                                let tx_local_id_in_this_msjc =
                                    self.msjcs[msjc_idx].add_txabundance(tx_abd);
                                tx_abd.add_msjc_mono_exonic(msjc_idx, tx_local_id_in_this_msjc);
                            }
                        }
                    } else {
                        continue;
                    }
                }

                local_mono_exonic_idx += 1;
            } else {
                let mut fsm_misoforms_candidates: Vec<&Vec<MergedIsoformOffsetPtr>> = Vec::new();
                // at least mapping with one sj to be considered as msjc
                let mut partial_msjc_map: HashMap<u64, &MergedIsoformOffsetPtr> =
                    HashMap::default();

                for sj_pair in tx_abd.sj_pairs.iter() {
                    let sj_start_idx = self
                        .position_tx_abd_map
                        .get(&sj_pair.0)
                        .unwrap_or_else(|| panic!("Can not get position index for {}", sj_pair.0));
                    let sj_end_idx = self
                        .position_tx_abd_map
                        .get(&sj_pair.1)
                        .unwrap_or_else(|| panic!("Can not get position index for {}", sj_pair.1));

                    fsm_misoforms_candidates.push(&all_res[*sj_start_idx]);
                    fsm_misoforms_candidates.push(&all_res[*sj_end_idx]);
                    quick_find_partial_msjc_exclude_read_st_end(
                        &all_res[*sj_start_idx],
                        &all_res[*sj_end_idx],
                        &mut partial_msjc_map,
                    );
                }

                let fsm_msjc_offsets = find_fsm(fsm_misoforms_candidates, tx_abd.sj_pairs.len());

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
                    let fsm_msjc_rec = archive_cache.read_bytes(&fsm);

                    tx_abd.add_fsm_misoform(&fsm_msjc_rec);
                }

                // add msjcs

                for (offset, msjc_ptr) in partial_msjc_map.into_iter() {
                    let msjc_idx = if let Some(&msjc_idx) = msjc_map_global.get(&offset) {
                        msjc_idx
                    } else {
                        let msjc = MSJC::new(
                            offset,
                            dbinfo.get_size(),
                            &archive_cache.read_bytes(&msjc_ptr),
                        );
                        self.msjcs.push(msjc);
                        let grouped_tx_msjc_idx = self.msjcs.len() - 1;
                        msjc_map_global.insert(offset, grouped_tx_msjc_idx);
                        grouped_tx_msjc_idx
                    };

                    let msjc = &mut self.msjcs[msjc_idx];
                    let (is_all_matched, first_match_pos, matched_count) =
                        msjc.splice_junctions_aligned_to_tx(&tx_abd.sj_pairs, cli.flank);

                    // consider adjust the match logic here and compare the correlations with ground truth
                    if is_all_matched {
                        if msjc.splice_junctions_vec.len() >= 1 {
                            if matched_count >= 1 {
                                let tx_local_id_in_this_msjc = msjc.add_txabundance(tx_abd);
                                tx_abd.add_msjc(
                                    msjc_map_global[&offset],
                                    tx_local_id_in_this_msjc,
                                    first_match_pos as usize,
                                    matched_count,
                                );
                            }
                        }
                    }
                }
            }
        }

        // self.msjcs = msjc_map_global;
    }

    pub fn prepare_em(&mut self, cli: &AnnIsoCli) {
        for msjc in self.msjcs.iter_mut() {
            msjc.prepare_em(self.tx_abundances.len());
            assert!(msjc.nonzero_count() > 0);
        }

        for txabd in self.tx_abundances.iter_mut() {
            txabd.prepare_em(&self.msjcs, cli);
        }

        // clean up the position_tx_abd_map to save memory
        self.position_tx_abd_map.clear();
        self.position_tx_abd_map.shrink_to_fit();
    }

    pub fn em(&mut self, cli: &AnnIsoCli) {
        if cli.verbose {
            info!(
                "Starting EM for grouped tx with {} transcripts and {} msjcs, memory usage: {:.2} MB",
                self.tx_abundances.len(),
                self.msjcs.len(),
                self.get_mem_size() as f64 / 1024.0 / 1024.0
            );
        }

        let start_time = std::time::Instant::now();

        for i in 0..cli.em_max_iter {
            if cli.verbose {
                info!("====== Starting EM iteration {} ======", i + 1);
            }

            // E step

            // info!("EM iteration {}: E step", i + 1);
            self.msjcs.par_iter_mut().for_each(|msjc| {
                msjc.e_step(&self.tx_abundances);
                // dbg!(msjc.)
            });

            // M step

            self.tx_abundances
                .par_iter_mut()
                .with_min_len(100)
                .for_each(|tx_abd| {
                    tx_abd.m_step(&self.msjcs, cli);
                });

            let all_converged = self
                .tx_abundances
                .par_iter_mut()
                .map(|tx_abd| tx_abd.check_convergence(cli.em_conv_min_diff, cli))
                .all(|x| x);

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

    }
}

pub fn quick_find_partial_msjc_exclude_read_st_end<'a>(
    a: &'a Vec<MergedIsoformOffsetPtr>,
    b: &'a Vec<MergedIsoformOffsetPtr>,
    collection: &mut HashMap<u64, &'a MergedIsoformOffsetPtr>,
) {

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

pub fn find_fsm(
    vecs: Vec<&Vec<MergedIsoformOffsetPtr>>,
    sj_pairs_n: usize,
) -> Vec<MergedIsoformOffsetPtr> {
    if vecs.is_empty() {
        return vec![];
    }

    let n_splice_sites = sj_pairs_n * 2; // each splice junction has two positions
    let mut result = vecs[0].clone();

    for v in vecs.iter().skip(1) {
        result = intersect_sorted(&result, v);

        if result.is_empty() {
            break;
        }
    }

    // filter out the MergedIsoformOffsetPtr dose not exeactly match n_splice
    result.retain(|x| x.n_splice_sites as usize == n_splice_sites);

    result
}

#[derive(Debug, Serialize, Deserialize, Clone)]
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
    pub is_mono_exonic: bool,
    sample_size: usize,
    inv_pt: f32, // scaling factor based solely on the complexity of transcript, same role as effective length,
    // is_ok: Vec<bool>,
    pub msjc_ids: Vec<(usize, usize)>, // (msjc id, index in msjc)
    pub is_need_em: bool,
    pub alive_samples_bitmap: Vec<u8>, // bitmap to indicate which samples are alive for em
    pub covered_sj_bits: Vec<u8>,      // bitmap to indicate which splice junctions are covered
}

impl TxAbundance {
    pub fn new(txid: usize, sample_size: usize, tx: &Transcript, cli: &AnnIsoCli) -> TxAbundance {
        // calclulate pt
        let nsj = tx.splice_junc.len() + cli.em_effective_len_coef; // J + 1
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
            msjc_ids: Vec::new(),
            is_need_em: true,
            alive_samples_bitmap: vec![0; (sample_size + 7) / 8], // initialize all samples as alive
            is_mono_exonic: tx.is_mono_exonic,
            covered_sj_bits: vec![0; (tx.get_splice_junction_pairs().len() + 7) / 8], // each splice junction has two positions
        }
    }

    pub fn add_fsm_misoform(&mut self, misoform: &MergedIsoform) {
        for sid in 0..self.sample_size {
            self.fsm_abundance[sid] += misoform.sample_evidence_arr[sid] as f32;
        }
    }

    pub fn add_msjc(
        &mut self,
        msjc_idx: usize,
        tx_local_id: usize,
        first_match_pos: usize,
        matched_count: usize,
    ) {
        self.msjc_ids.push((msjc_idx, tx_local_id));

        // mark covered splice junctions
        // first_match_pos is the index in self.sj_pairs where the first matched splice junction is located
        for i in 0..matched_count {
            let sj_idx = first_match_pos + i;
            self.covered_sj_bits[sj_idx / 8] |= 1 << (sj_idx % 8);
        }
    }

    pub fn add_msjc_mono_exonic(&mut self, msjc_idx: usize, tx_local_id: usize) {
        self.msjc_ids.push((msjc_idx, tx_local_id));
    }

    pub fn is_all_sj_covered(&self, cli: &AnnIsoCli) -> bool {
        // must match all splice junctions including the first and last bits
        if cli.only_fully_covered_tx {
            let n_sj = self.sj_pairs.len();
            let n_bytes = (n_sj + 7) / 8;
            for i in 0..n_bytes {
                let byte = self.covered_sj_bits[i];
                if i == n_bytes - 1 {
                    // last byte
                    let remaining_bits = n_sj % 8;
                    let mask = if remaining_bits == 0 {
                        0xFF
                    } else {
                        (1 << remaining_bits) - 1
                    };
                    if byte & mask != mask {
                        return false;
                    }
                } else {
                    if byte != 0xFF {
                        return false;
                    }
                }
            }
            true
        } else {
            // dont consider first and last splice junctions

            let n_sj = self.sj_pairs.len();
            if n_sj <= 2 {
                return true;
            } else {
                let mut all_covered = true;
                for sj_idx in 1..(n_sj - 1) {
                    let byte_idx = sj_idx / 8;
                    let bit_idx = sj_idx % 8;
                    if (self.covered_sj_bits[byte_idx] & (1 << bit_idx)) == 0 {
                        all_covered = false;
                        break;
                    }
                }
                all_covered
            }
        }
    }

    pub fn prepare_em(&mut self, msjc_vec: &Vec<MSJC>, cli: &AnnIsoCli) {
        // check all associated MSJCs, and mark the samples that have zero coverage in all MSJCs
        // set the self.alive_samples_bitmap accordingly
        // iterately process each msjc, set the bitmap

        // update the inv_pt if its a mono exonic transcript
        if self.is_mono_exonic {
            if self.msjc_ids.is_empty() {
                debug!("Warning: mono exonic transcript {} has no associated msjc to calculate effective length, keep inv_pt as 1.0", self.orig_tx_id);
            } else {
                let mut avg_msjc_len = 0.0;
                let mut count = 0;

                for (msjc_idx, _tx_local_id) in self.msjc_ids.iter() {
                    // info!(
                    //     "xsafe 1"
                    // );
                    let msjc = unsafe { msjc_vec.get_unchecked(*msjc_idx) };
                    avg_msjc_len += msjc.get_effective_length() as f32;
                    count += 1;
                }
                if count > 0 && avg_msjc_len > 0.0 {
                    avg_msjc_len /= count as f32;
                    self.inv_pt = 1.0 / (avg_msjc_len as f32);
                } else {
                    debug!("Warning: mono exonic transcript {} has no associated msjc to calculate effective length, keep inv_pt as 1.0", self.orig_tx_id);
                }
            }
        } else {
            if !self.is_all_sj_covered(cli) {
                // dont need em
                self.is_need_em = false;

                unsafe {
                    std::ptr::write_bytes(self.abundance_cur.as_mut_ptr(), 0, self.sample_size);
                }
            } else {
                for (msjc_idx, _tx_local_id) in self.msjc_ids.iter() {

                    let msjc = unsafe { msjc_vec.get_unchecked(*msjc_idx) };
                    for &sid in msjc.nonzero_sample_indices.iter() {
                        let byte_idx = sid / 8;
                        let bit_idx = sid % 8;
                        self.alive_samples_bitmap[byte_idx] |= 1 << bit_idx;
                    }
                }
            }
        }
    }

    pub fn m_step(&mut self, msjc_vec: &Vec<MSJC>, cli: &AnnIsoCli) {
        if self.msjc_ids.is_empty() {
            // no msjc support, set abundance to 0
            if self.is_need_em {
                unsafe {
                    std::ptr::write_bytes(self.abundance_cur.as_mut_ptr(), 0, self.sample_size);
                }
                self.is_need_em = false;
                return;
            } else {
                return;
            }
        }

        if cli.verbose {
            info!(
                "Updating TxAbundance id {},tx_id {:?},msjc length: {}",
                self.id,
                self.orig_tx_id,
                self.msjc_ids.len()
            );
        }

        let sample_size = self.sample_size;

        unsafe {
            for sid in 0..sample_size {
                *self.abundance_prev.get_unchecked_mut(sid) =
                    *self.abundance_cur.get_unchecked(sid);
            }

            std::ptr::write_bytes(self.abundance_cur.as_mut_ptr(), 0, sample_size);

            for (msjc_idx, tx_local_id) in self.msjc_ids.iter() {
                let msjc = msjc_vec.get_unchecked(*msjc_idx);
                let k = msjc.nonzero_count();

                let mut i = 0;
                let end = k & !3;

                while i < end {
                    let sid0 = *msjc.nonzero_sample_indices.get_unchecked(i);
                    let sid1 = *msjc.nonzero_sample_indices.get_unchecked(i + 1);
                    let sid2 = *msjc.nonzero_sample_indices.get_unchecked(i + 2);
                    let sid3 = *msjc.nonzero_sample_indices.get_unchecked(i + 3);

                    let cov0 = *msjc.cov_vec.get_unchecked(i);
                    let cov1 = *msjc.cov_vec.get_unchecked(i + 1);
                    let cov2 = *msjc.cov_vec.get_unchecked(i + 2);
                    let cov3 = *msjc.cov_vec.get_unchecked(i + 3);

                    let base = *tx_local_id * k;
                    let gamma0 = *msjc.gamma_data.get_unchecked(base + i);
                    let gamma1 = *msjc.gamma_data.get_unchecked(base + i + 1);
                    let gamma2 = *msjc.gamma_data.get_unchecked(base + i + 2);
                    let gamma3 = *msjc.gamma_data.get_unchecked(base + i + 3);

                    *self.abundance_cur.get_unchecked_mut(sid0) += cov0 * gamma0;
                    *self.abundance_cur.get_unchecked_mut(sid1) += cov1 * gamma1;
                    *self.abundance_cur.get_unchecked_mut(sid2) += cov2 * gamma2;
                    *self.abundance_cur.get_unchecked_mut(sid3) += cov3 * gamma3;

                    i += 4;
                }

                while i < k {
                    let sid = *msjc.nonzero_sample_indices.get_unchecked(i);
                    let cov = *msjc.cov_vec.get_unchecked(i);
                    let gamma = *msjc.gamma_data.get_unchecked(*tx_local_id * k + i);

                    *self.abundance_cur.get_unchecked_mut(sid) += cov * gamma;
                    i += 1;
                }
            }
        }
    }

    pub fn check_convergence(&mut self, tol: f32, cli: &AnnIsoCli) -> bool {
        // 这里必须check 所有的postive的样的abundance，而不是所有，否则，初始化的abundance都是1.0，realdiff都是1，永远不会converge。

        // only check the alive samples
        let mut is_converged = true;
        for sid in 0..self.sample_size {
            let byte_idx = sid / 8;
            let bit_idx = sid % 8;
            if (self.alive_samples_bitmap[byte_idx] & (1 << bit_idx)) != 0 {
                let prev = self.abundance_prev[sid];
                let cur = self.abundance_cur[sid];
                let diff = (cur - prev).abs();
                let denom = if prev.abs() < 1e-6 { 1.0 } else { prev.abs() };
                let rel_diff = diff / denom;
                if cli.verbose {
                    info!(
                        "TxAbundance id {}, tx_id {:?}, sample {} fsm {} em prev {:.6} em cur {:.6},diff {:.6} rel_diff {:.6}",
                        self.id, self.orig_tx_id, sid, self.fsm_abundance[sid], prev, cur, diff, rel_diff
                    );
                }
                if rel_diff > tol {
                    is_converged = false;
                }
            }
        }

        is_converged
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        TxAbundanceView::encode(self)
    }
}

pub struct TxAbundanceView {
    _orig_idx: u32,
    orig_tx_id: Vec<u8>,
    orig_gene_id: Vec<u8>,
    orig_tx_len: u32,
    orig_n_exon: u32,
    orig_chrom: Vec<u8>,
    orig_start: u64,
    orig_end: u64,
    orig_attrs: Vec<u8>,
    fsm_abundance: Vec<f32>,
    em_abundance: Vec<f32>,
}

impl TxAbundanceView {
    pub fn encode(txabd: &TxAbundance) -> Vec<u8> {
        let mut buf = Vec::new();

        buf.extend_from_slice(&txabd.orig_idx.to_le_bytes());
        buf.extend_from_slice(&txabd.orig_tx_len.to_le_bytes());
        buf.extend_from_slice(&txabd.orig_n_exon.to_le_bytes());
        buf.extend_from_slice(&txabd.orig_start.to_le_bytes());
        buf.extend_from_slice(&txabd.orig_end.to_le_bytes());

        // txid,geneid,chrom,attr
        buf.extend_from_slice(&(txabd.orig_tx_id.len() as u32).to_le_bytes());
        buf.extend_from_slice(txabd.orig_tx_id.as_bytes());
        buf.extend_from_slice(&(txabd.orig_gene_id.len() as u32).to_le_bytes());
        buf.extend_from_slice(txabd.orig_gene_id.as_bytes());
        buf.extend_from_slice(&(txabd.orig_chrom.len() as u32).to_le_bytes());
        buf.extend_from_slice(txabd.orig_chrom.as_bytes());
        buf.extend_from_slice(&(txabd.orig_attrs.len() as u32).to_le_bytes());
        buf.extend_from_slice(txabd.orig_attrs.as_bytes());

        // fsm_abundance,em_abundance

        buf.extend_from_slice(&(txabd.fsm_abundance.len() as u32).to_le_bytes());
        for &v in &txabd.fsm_abundance {
            buf.extend_from_slice(&v.to_le_bytes());
        }

        buf.extend_from_slice(&(txabd.abundance_cur.len() as u32).to_le_bytes());
        for &v in &txabd.abundance_cur {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut pos = 0;

        let orig_idx = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_tx_len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_n_exon = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_start = u64::from_le_bytes(data[pos..pos + 8].try_into()?);
        pos += 8;
        let orig_end = u64::from_le_bytes(data[pos..pos + 8].try_into()?);
        pos += 8;

        // txid,geneid,chrom,attr
        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_tx_id = &data[pos..pos + len as usize];
        pos += len as usize;
        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_gene_id = &data[pos..pos + len as usize];
        pos += len as usize;
        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_chrom = &data[pos..pos + len as usize];
        pos += len as usize;

        // attr
        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        let orig_attrs = &data[pos..pos + len as usize];

        pos += len as usize;

        // fsm em abundance
        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;

        let mut fsm_abundance = Vec::new();
        let mut em_abundance = Vec::new();

        for _ in 0..len {
            fsm_abundance.push(f32::from_le_bytes(data[pos..pos + 4].try_into()?));
            pos += 4;
        }

        let len = u32::from_le_bytes(data[pos..pos + 4].try_into()?);
        pos += 4;
        for _ in 0..len {
            em_abundance.push(f32::from_le_bytes(data[pos..pos + 4].try_into()?));
            pos += 4;
        }

        Ok(TxAbundanceView {
            _orig_idx: orig_idx,
            orig_tx_id: orig_tx_id.to_vec(),
            orig_gene_id: orig_gene_id.to_vec(),
            orig_tx_len,
            orig_n_exon,
            orig_chrom: orig_chrom.to_vec(),
            orig_start,
            orig_end,
            orig_attrs: orig_attrs.to_vec(),
            fsm_abundance,
            em_abundance,
        })
    }

    pub fn get_positive_samples(&self, min_read: f32) -> usize {
        let mut count = 0;
        for (abd1, abd2) in self.em_abundance.iter().zip(&self.fsm_abundance) {
            let total = abd1 + abd2;
            if total >= min_read {
                count += 1;
            }
        }
        count
    }

    pub fn write_line_directly(
        &self,
        global_stats: &GlobalStats,
        dbinfo: &DatasetInfo,
        cli: &AnnIsoCli,
        tableout: &mut TableOutput,
    ) -> Result<()> {
        tableout.write_bytes(&self.orig_chrom)?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(self.orig_start.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(self.orig_end.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(self.orig_tx_len.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(self.orig_n_exon.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(&self.orig_tx_id)?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(&self.orig_gene_id)?;
        tableout.write_bytes(b"\t")?;
        // confidence scores
        let u32_arr = self
            .fsm_abundance
            .iter()
            .zip(&self.em_abundance)
            .map(|(&fsm, em)| fsm + em)
            .collect::<Vec<f32>>();

        let confscore =
            utils::calc_confidence_f32(&u32_arr, dbinfo.get_size(), &global_stats.fsm_em_total);

        tableout.write_bytes(confscore.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;

        let positive_samples = self.get_positive_samples(cli.min_read as f32);
        if positive_samples != 0 {
            tableout.write_bytes(b"yes")?;
        } else {
            tableout.write_bytes(b"no")?;
        }
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(cli.min_read.to_string().as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(format!("{}/{}", positive_samples, dbinfo.get_size()).as_bytes())?;
        tableout.write_bytes(b"\t")?;
        tableout.write_bytes(&self.orig_attrs)?;
        tableout.write_bytes(b"\t")?;
        // let format_str = tableout.format_str.clone();
        // tableout.write_bytes(format_str.as_bytes())?;
        tableout.write_format_str()?;
        for sid in 0..dbinfo.get_size() {
            let fsm = self.fsm_abundance[sid];
            let em = match self.em_abundance[sid] >= cli.min_em_abundance {
                true => self.em_abundance[sid],
                false => 0.0,
            };

            let total_abd = fsm + em;
            let total_cov = global_stats.fsm_em_total[sid];

            let fsm_cpm = utils::calc_cpm_f32(&fsm, &total_cov);
            let em_cpm = utils::calc_cpm_f32(&em, &total_cov);
            let total_cpm = utils::calc_cpm_f32(&total_abd, &total_cov);

            tableout.write_bytes(b"\t")?;
            tableout.write_bytes(total_cpm.to_string().as_bytes())?;
            tableout.write_bytes(b":")?;
            tableout.write_bytes(total_abd.to_string().as_bytes())?;
            tableout.write_bytes(b":")?;
            tableout.write_bytes(fsm_cpm.to_string().as_bytes())?;
            tableout.write_bytes(b":")?;
            tableout.write_bytes(fsm.to_string().as_bytes())?;
            tableout.write_bytes(b":")?;
            tableout.write_bytes(em_cpm.to_string().as_bytes())?;
            tableout.write_bytes(b":")?;
            tableout.write_bytes(em.to_string().as_bytes())?;

        }
        tableout.write_bytes(b"\n")?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct MSJC {
    #[allow(unused)]
    id: u64,
    #[allow(unused)]
    sample_size: usize,

    nonzero_sample_indices: Vec<usize>, // length: K (K <= sample_size), values are sample indices with non-zero coverage
    cov_vec: Vec<f32>,                  // length: K, coverage values for non-zero samples

    // Transcript IDs
    txids: Vec<usize>, // length: T

    // Gamma matrix in a flattened format
    // gamma[tx_local_id][nonzero_idx] = gamma_data[tx_local_id * K + nonzero_idx]
    gamma_data: Vec<f32>, // length: T * K

    // Buffer to hold totals during E-step
    totals_buffer: Vec<f32>, // lenghth: K

    splice_junctions_vec: Vec<(u64, u64)>,
}

impl MSJC {
    pub fn new(id: u64, sample_size: usize, misoform: &MergedIsoform) -> MSJC {
        let mut nonzero_indices = Vec::new();
        let mut cov_vec = Vec::new();

        for (idx, &cov) in misoform.sample_evidence_arr.iter().enumerate() {
            if cov > 0 {
                nonzero_indices.push(idx);
                cov_vec.push(cov as f32);
            }
        }

        MSJC {
            id,
            sample_size,
            nonzero_sample_indices: nonzero_indices,
            cov_vec,
            txids: Vec::new(),
            gamma_data: Vec::new(),
            totals_buffer: Vec::new(),
            splice_junctions_vec: misoform.splice_junctions_vec.clone(),
        }
    }

    pub fn add_txabundance(&mut self, txabd: &TxAbundance) -> usize {
        self.txids.push(txabd.id);
        self.txids.len() - 1
    }

    pub fn get_effective_length(&self) -> u32 {
        // start -end
        if self.splice_junctions_vec.is_empty() {
            0
        } else {
            let start = self
                .splice_junctions_vec
                .first()
                .expect("MSJC can not get start of splice junction vec")
                .0;
            let end = self
                .splice_junctions_vec
                .last()
                .expect("MSJC can not get end of splice junction vec")
                .1;
            (end - start) as u32
        }
    }

    pub fn prepare_em(&mut self, tx_count: usize) {
        let k = self.nonzero_sample_indices.len();
        self.gamma_data = vec![0.0; k * tx_count];
        self.totals_buffer = vec![0.0; k];
    }

    pub fn check_mono_exon_fsm(&self, txabd: &TxAbundance, cli: &AnnIsoCli) -> bool {
        if self.splice_junctions_vec.is_empty() {
            return false;
        }

        if self.splice_junctions_vec.len() > 1 {
            return false;
        }

        let is_fsm = true;

        if self
            .splice_junctions_vec
            .first()
            .unwrap()
            .0
            .abs_diff(txabd.orig_start)
            > cli.flank as u64
        {
            return false;
        }

        if self
            .splice_junctions_vec
            .last()
            .unwrap()
            .1
            .abs_diff(txabd.orig_end)
            > cli.flank as u64
        {
            return false;
        }

        is_fsm
    }

    pub fn check_msjc_belong_to_tx(&self, txabd: &TxAbundance, _cli: &AnnIsoCli) -> bool {
        // the MJSC must with in the transcript region
        if self.splice_junctions_vec.is_empty() {
            return false;
        }

        if self.splice_junctions_vec.len() > 1 {
            return false;
        }

        if self.splice_junctions_vec.first().unwrap().0 < txabd.orig_start {
            return false;
        }

        if self.splice_junctions_vec.last().unwrap().1 > txabd.orig_end {
            return false;
        }

        true
    }

    #[inline(always)]
    pub fn nonzero_count(&self) -> usize {
        self.nonzero_sample_indices.len()
    }

    #[inline(always)]
    fn e_step(&mut self, txabds: &Vec<TxAbundance>) {
        let k = self.nonzero_sample_indices.len();
        let txabds_slice = txabds.as_slice();
        unsafe {
            for i in 0..k {
                *self.totals_buffer.get_unchecked_mut(i) = 0.0;
            }

            for &txid in self.txids.iter() {
                let txabd = txabds_slice.get_unchecked(txid);

                if txabd.is_need_em == false {
                    continue;
                }

                let inv_pt = txabd.inv_pt;
                let abd_slice = txabd.abundance_cur.as_slice();

                let mut i = 0;
                let end = k & !3;

                while i < end {
                    let sid0 = *self.nonzero_sample_indices.get_unchecked(i);

                    let sid1 = *self.nonzero_sample_indices.get_unchecked(i + 1);

                    let sid2 = *self.nonzero_sample_indices.get_unchecked(i + 2);

                    let sid3 = *self.nonzero_sample_indices.get_unchecked(i + 3);

                    *self.totals_buffer.get_unchecked_mut(i) +=
                        *abd_slice.get_unchecked(sid0) * inv_pt;

                    *self.totals_buffer.get_unchecked_mut(i + 1) +=
                        *abd_slice.get_unchecked(sid1) * inv_pt;

                    *self.totals_buffer.get_unchecked_mut(i + 2) +=
                        *abd_slice.get_unchecked(sid2) * inv_pt;

                    *self.totals_buffer.get_unchecked_mut(i + 3) +=
                        *abd_slice.get_unchecked(sid3) * inv_pt;

                    i += 4;
                }

                // info!("step3");
                while i < k {
                    // info!("step3.1 i={}", i);
                    let sid = *self.nonzero_sample_indices.get_unchecked(i);
                    // info!("step3.2 i={}", i);
                    *self.totals_buffer.get_unchecked_mut(i) +=
                        *abd_slice.get_unchecked(sid) * inv_pt;
                    i += 1;
                }
            }

            for (tx_local_id, &txid) in self.txids.iter().enumerate() {
                let txabd = txabds_slice.get_unchecked(txid);
                let inv_pt = txabd.inv_pt;
                let abd_slice = txabd.abundance_cur.as_slice();

                let mut i = 0;
                let end = k & !3;

                while i < end {
                    // info!("step4.1 i={}", i);
                    let sid0 = *self.nonzero_sample_indices.get_unchecked(i);
                    // info!("step4.2 i={}", i+1);
                    let sid1 = *self.nonzero_sample_indices.get_unchecked(i + 1);
                    // info!("step4.3 i={}", i+2);
                    let sid2 = *self.nonzero_sample_indices.get_unchecked(i + 2);
                    // info!("step4.4 i={}", i+3);
                    let sid3 = *self.nonzero_sample_indices.get_unchecked(i + 3);

                    let compute_gamma = |nonzero_idx: usize, sid: usize| {
                        // info!("step4.5 i={} sid={}", i, sid);
                        let total = *self.totals_buffer.get_unchecked(nonzero_idx);
                        if total > 0.0 {
                            *abd_slice.get_unchecked(sid) * inv_pt / total
                        } else {
                            0.0
                        }
                    };

                    let base = tx_local_id * k;
                    // info!("step4.6 i={}", i);
                    *self.gamma_data.get_unchecked_mut(base + i) = compute_gamma(i, sid0);
                    // info!("step4.7 i={}", i+1);
                    *self.gamma_data.get_unchecked_mut(base + i + 1) = compute_gamma(i + 1, sid1);
                    // info!("step4.8 i={}", i+2);
                    *self.gamma_data.get_unchecked_mut(base + i + 2) = compute_gamma(i + 2, sid2);
                    // info!("step4.9 i={}", i+3);
                    *self.gamma_data.get_unchecked_mut(base + i + 3) = compute_gamma(i + 3, sid3);

                    i += 4;
                }

                while i < k {
                    // info!("step5.1 i={}", i);
                    let sid = *self.nonzero_sample_indices.get_unchecked(i);
                    // info!("step5.2 i={}", i);
                    let total = *self.totals_buffer.get_unchecked(i);
                    let gamma = if total > 0.0 {
                        // info!("step5.3 i={}", i);
                        *abd_slice.get_unchecked(sid) * inv_pt / total
                    } else {
                        0.0
                    };
                    // info!("step5.4 i={}", i);
                    *self.gamma_data.get_unchecked_mut(tx_local_id * k + i) = gamma;
                    i += 1;
                }
            }
        }
    }

    pub fn splice_junctions_aligned_to_tx(
        &self,
        reference_sjs: &Vec<(u64, u64)>,
        flank: u64,
    ) -> (bool, i32, usize) {
        let mut matched_count = 0usize;
        let mut first_match_pos = -1i32;
        let mut is_all_matched = true;

        if self.splice_junctions_vec.len() > reference_sjs.len() {
            is_all_matched = false;
            return (is_all_matched, first_match_pos, matched_count);
        }

        let mut qidx = 0;
        let mut sidx = 0;

        loop {
            if qidx >= reference_sjs.len() || sidx >= self.splice_junctions_vec.len() {
                if sidx < self.splice_junctions_vec.len() {
                    is_all_matched = false;
                }
                break;
            }

            let query_sj = &reference_sjs[qidx];
            let isoform_sj = &self.splice_junctions_vec[sidx];

            if (isoform_sj.0.abs_diff(query_sj.0) <= flank)
                && (isoform_sj.1.abs_diff(query_sj.1) <= flank)
            {
                matched_count += 1;
                if first_match_pos == -1 {
                    first_match_pos = qidx as i32;
                }
                qidx += 1;
                sidx += 1;
            } else {
                if first_match_pos == -1 {
                    qidx += 1;
                } else {
                    // skip sjs
                    is_all_matched = false;
                    break;
                }
            }
        }

        (is_all_matched, first_match_pos, matched_count)
    }
}

type OrigIdx = u32;
// type Offset = u64;
// type Length = u64;

/// write to temporary output manager with sharding,
/// two files will be generated, one is the data file, another is the index file
pub struct TmpOutputManager {
    // pub curr_offset: u64,
    pub curr_processed_idx: usize,
    pub shards_offsets_vec: Vec<Vec<OrigIdx>>,
    pub records: Vec<TxAbundance>,
    pub file_base_path: PathBuf,
    pub curr_shard_idx: usize,
    pub shard_capacity: usize,
    pub heap: BinaryHeap<Reverse<HeapEntry>>, // (orig_idx, shard_idx)
    pub data_readers: Vec<BufReader<File>>,
    pub shard_read_indices: Vec<usize>,
    pub data_length: [u8; 8],
    pub data_record: Vec<u8>,
}

impl TmpOutputManager {
    pub fn new(file_path: &PathBuf, cli: &AnnIsoCli) -> Self {
        TmpOutputManager {
            // curr_offset: 0,
            curr_processed_idx: 0,
            shards_offsets_vec: Vec::new(),
            records: Vec::new(),
            file_base_path: file_path.clone(),
            curr_shard_idx: 0,
            shard_capacity: cli.output_tmp_shard_counts,
            heap: BinaryHeap::new(),
            data_readers: Vec::new(),
            shard_read_indices: Vec::new(),
            data_length: [0; 8],
            data_record: Vec::new(),
        }
    }

    pub fn dump_grouped_tx(&mut self, grouped_tx: &mut GroupedTx) {
        // for tx_abd in grouped_tx.tx_abundances.iter() {
        //     // self.records.push(tx_abd.clone());

        // }

        let txs: Vec<TxAbundance> = std::mem::take(&mut grouped_tx.tx_abundances);
        self.records.extend(txs);

        if self.records.len() >= self.shard_capacity {
            // self.curr_offset = 0;
            // sort the records by orig_idx
            self.records.sort_by_key(|tx_abd| tx_abd.orig_idx);
            // write to file
            let shard_file_path = self
                .file_base_path
                .with_extension(format!("tmp.{}", self.curr_shard_idx));

            // create a buffer writer
            let mut file_h = std::fs::File::create(&shard_file_path).unwrap_or_else(|_| {
                panic!(
                    "Could not create temporary output file {}",
                    shard_file_path.display()
                )
            });
            let mut bufwriter = BufWriter::new(&mut file_h);

            let mut curr_shard_offsets = Vec::new();

            info!(
                "Writing shard {} with {} records to {}",
                self.curr_shard_idx,
                self.records.len(),
                shard_file_path.display()
            );

            for tx_abd in self.records.iter() {
                let bytes = tx_abd.to_bytes();
                let length = bytes.len() as u64;
                bufwriter
                    .write_all(&length.to_le_bytes())
                    .unwrap_or_else(|_| {
                        panic!(
                            "Could not write byte lengthto temporary output file {}",
                            shard_file_path.display()
                        )
                    });
                bufwriter.write_all(&bytes).unwrap_or_else(|_| {
                    panic!(
                        "Could not write to byte data temporary output file {}",
                        shard_file_path.display()
                    )
                });
                curr_shard_offsets.push(tx_abd.orig_idx);
                // self.curr_offset += length + 8;
            }

            self.shards_offsets_vec.push(curr_shard_offsets);
            bufwriter.flush().unwrap_or_else(|_| {
                panic!(
                    "Could not flush to temporary output file {}",
                    shard_file_path.display()
                )
            });
            self.curr_shard_idx += 1;
            self.records.clear();
        }
    }

    pub fn finish(&mut self) {
        // write remaining records

        if self.records.len() > 0 {
            // self.curr_offset = 0;
            // sort the records by orig_idx
            self.records.sort_by_key(|tx_abd| tx_abd.orig_idx);
            // write to file
            let shard_file_path = self
                .file_base_path
                .with_extension(format!("tmp.{}", self.curr_shard_idx));

            // create a buffer writer
            let mut file_h = std::fs::File::create(&shard_file_path).unwrap_or_else(|_| {
                panic!(
                    "Could not create temporary output file {}",
                    shard_file_path.display()
                )
            });
            let mut bufwriter = BufWriter::new(&mut file_h);

            let mut curr_shard_offsets = Vec::new();

            info!(
                "Writing final shard {} with {} records to {}",
                self.curr_shard_idx,
                self.records.len(),
                shard_file_path.display()
            );

            for tx_abd in self.records.iter() {
                let bytes = tx_abd.to_bytes();
                let length = bytes.len() as u64;
                bufwriter
                    .write_all(&length.to_le_bytes())
                    .unwrap_or_else(|_| {
                        panic!(
                            "Could not write byte length to temporary output file {}",
                            shard_file_path.display()
                        )
                    });
                bufwriter.write_all(&bytes).unwrap_or_else(|_| {
                    panic!(
                        "Could not write to byte data temporary output file {}",
                        shard_file_path.display()
                    )
                });
                curr_shard_offsets.push(tx_abd.orig_idx);
                // self.curr_offset += length + 8;
            }

            self.shards_offsets_vec.push(curr_shard_offsets);
            bufwriter.flush().unwrap_or_else(|_| {
                panic!(
                    "Could not flush to temporary output file {}",
                    shard_file_path.display()
                )
            });

            self.curr_shard_idx += 1;
            self.records.clear();
        }

        // prepare mmap readers
        for shard_idx in 0..self.curr_shard_idx {
            let shard_file_path = self
                .file_base_path
                .with_extension(format!("tmp.{}", shard_idx));
            let file_h = File::open(&shard_file_path).unwrap_or_else(|_| {
                panic!(
                    "Could not open temporary output file {} for mmap",
                    shard_file_path.display()
                )
            });
            //             info!(
            //     "xsafe 5"
            // );
            // make bufreader for each shard
            let bufreader = BufReader::new(file_h);
            self.data_readers.push(bufreader);
        }
        // 初始化每个 shard 的读取索引
        self.shard_read_indices = vec![0; self.curr_shard_idx];

        // init the heap - 将每个 shard 的第一个元素加入堆
        for shard_idx in 0..self.curr_shard_idx {
            if !self.shards_offsets_vec[shard_idx].is_empty() {
                let orig_idx = self.shards_offsets_vec[shard_idx][0];
                self.heap.push(Reverse(HeapEntry {
                    orig_idx,
                    shard_idx,
                }));
            }
        }

        info!(
            "K-way merge initialized with {} shards, heap size: {}",
            self.curr_shard_idx,
            self.heap.len()
        );
    }

    pub fn clean_up(&mut self) -> Result<()> {
        // remove temporary files
        for shard_idx in 0..self.curr_shard_idx {
            let shard_file_path = self
                .file_base_path
                .with_extension(format!("tmp.{}", shard_idx));
            if std::fs::remove_file(&shard_file_path).is_err() {
                warn!(
                    "Could not remove temporary output file {}",
                    shard_file_path.display()
                );
            }
        }
        
        Ok(())
    }
}

impl Iterator for TmpOutputManager {
    type Item = TxAbundanceView;

    fn next(&mut self) -> Option<Self::Item> {
        // 从堆中取出最小的元素
        let Reverse(entry) = self.heap.pop()?;

        let HeapEntry { shard_idx, .. } = entry;

        let bufreader = &mut self.data_readers[shard_idx];
        // let mut data_length = vec![0u8; 8]; // u64
        self.data_length.fill(0);
        bufreader
            .read_exact(&mut self.data_length)
            .expect("Can't read length of tx_abd_view record");

        self.data_record
            .resize(u64::from_le_bytes(self.data_length) as usize, 0u8);
        bufreader
            .read_exact(&mut self.data_record)
            .expect("Can not read tx_abd_view record");

        let tx_abd_view =
            TxAbundanceView::from_bytes(&self.data_record).expect("Can't read TxabundanceView");

        self.shard_read_indices[shard_idx] += 1;
        let next_idx = self.shard_read_indices[shard_idx];

        if next_idx < self.shards_offsets_vec[shard_idx].len() {
            let next_orig_idx = self.shards_offsets_vec[shard_idx][next_idx];
            self.heap.push(Reverse(HeapEntry {
                orig_idx: next_orig_idx,
                shard_idx,
            }));
        }

        Some(tx_abd_view)
    }
}

#[derive(Eq, PartialEq)]
pub struct HeapEntry {
    orig_idx: OrigIdx,
    shard_idx: usize,
    // offset: Offset,
    // length: Length,
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.orig_idx
            .cmp(&other.orig_idx)
            .then_with(|| self.shard_idx.cmp(&other.shard_idx))
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl GetMemSize for TxAbundance {
    fn get_mem_size(&self) -> usize {
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
}

impl GetMemSize for GroupedTx {
    fn get_mem_size(&self) -> usize {
        let mut size = std::mem::size_of_val(&self.positions);

        for tx_ab in self.tx_abundances.iter() {
            size += tx_ab.get_mem_size();
        }

        for msjc in self.msjcs.iter() {
            size += msjc.get_mem_size();
        }

        size
    }
}

impl GetMemSize for MSJC {
    fn get_mem_size(&self) -> usize {
        std::mem::size_of::<u64>()
            + std::mem::size_of::<usize>()
            + self.nonzero_sample_indices.len() * std::mem::size_of::<usize>()
            + self.cov_vec.len() * std::mem::size_of::<f32>()
            + self.txids.len() * std::mem::size_of::<usize>()
            + self.gamma_data.len() * std::mem::size_of::<f32>()
            + self.totals_buffer.len() * std::mem::size_of::<f32>()
    }
}

impl GetMemSize for ChromGroupedTxManager {
    fn get_mem_size(&self) -> usize {
        let mut size = std::mem::size_of_val(&self.chrom);

        for grouped_tx in self.group_queries.iter() {
            size += grouped_tx.get_mem_size();
        }
        // + std::mem::size_of_val(&self._position_to_group_idx);
        size
    }
}
