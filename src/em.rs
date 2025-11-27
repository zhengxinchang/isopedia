use bincode::de;

use crate::{gtf::Transcript, isoform::MergedIsoform};

#[derive(Debug)]
pub struct TxAbundance {
    pub id: usize,
    pub sj_pairs: Vec<(u64, u64)>, // splice junction positions
    pub fsm_abundance: Vec<f32>,
    pub abundance_cur: Vec<f32>, // sample size length
    pub abundance_prev: Vec<f32>,
    sample_size: usize,
    inv_pt: f32, // scaling factor based solely on the complexity of transcript, same role as effective length,
    is_ok: Vec<bool>,
    msjc_ids: Vec<(usize, usize)>, // (msjc id, index in msjc)
}

impl TxAbundance {
    pub fn new(txid: usize, sample_size: usize, tx: &Transcript) -> TxAbundance {
        // calclulate pt
        let nsj = tx.splice_junc.len() + 10; // J + 1
        TxAbundance {
            id: txid,
            sj_pairs: tx.get_splice_junction_pairs(),
            fsm_abundance: vec![0.0; sample_size],
            abundance_cur: vec![1.0; sample_size],
            abundance_prev: vec![1.0; sample_size],
            sample_size,
            inv_pt: 1.0 / nsj as f32,
            is_ok: vec![true; sample_size],
            msjc_ids: Vec::new(),
        }
    }

    pub fn add_msjc(&mut self, msjc: &mut MSJC) {
        let local_id = msjc.add_txabundance(self);
        self.msjc_ids.push((msjc.id, local_id));
    }

    pub fn add_fsm_misoform(&mut self, misoform: &MergedIsoform) {
        self.fsm_abundance = misoform
            .sample_evidence_arr
            .iter()
            .map(|&x| x as f32)
            .collect();
    }

    pub fn get_mem_size(&self) -> usize {
        let size = std::mem::size_of_val(&self.id)
            + std::mem::size_of_val(&0f32) * self.fsm_abundance.len()
            + std::mem::size_of_val(&0f32) * self.abundance_cur.len()
            + std::mem::size_of_val(&0f32) * self.abundance_prev.len()
            + std::mem::size_of_val(&self.sample_size)
            + std::mem::size_of_val(&self.inv_pt)
            + std::mem::size_of_val(&self.is_ok)
            + std::mem::size_of_val(&self.msjc_ids);
        size
    }

    pub fn m_step(&mut self, msjc_vec: &Vec<MSJC>) {
        for sid in 0..self.abundance_cur.len() {
            let mut total_abd = 0.0f32;
            for (msjc_id, tx_local_id) in self.msjc_ids.iter() {
                let msjc = &msjc_vec[*msjc_id as usize];
                let idx = *tx_local_id * self.sample_size + sid;
                total_abd += msjc.cov_vec[sid] * msjc.gamma_vec[idx];
            }
            self.abundance_prev[sid] = self.abundance_cur[sid];
            self.abundance_cur[sid] = total_abd;
        }
    }
}

pub struct MSJC {
    id: usize,
    sample_size: usize,
    cov_vec: Vec<f32>,
    txids: Vec<usize>,
    gamma_vec: Vec<f32>, // length is same as txids
}

impl MSJC {
    pub fn new(id: usize, sample_size: usize, misoform: &MergedIsoform) -> MSJC {
        MSJC {
            id,
            sample_size,
            cov_vec: misoform
                .sample_evidence_arr
                .iter()
                .map(|&x| x as f32)
                .collect(),
            txids: Vec::new(),
            gamma_vec: Vec::new(),
        }
    }

    pub fn add_txabundance(&mut self, txabd: &TxAbundance) -> usize {
        self.txids.push(txabd.id);
        self.txids.len() - 1
    }

    pub fn prepare_em(&mut self) {
        self.gamma_vec = vec![0.0; self.sample_size * self.txids.len()];
    }

    pub fn get_mem_size(&self) -> usize {
        let size = std::mem::size_of_val(&self.id)
            + std::mem::size_of_val(&self.sample_size)
            + std::mem::size_of_val(&self.cov_vec)
            + std::mem::size_of_val(&self.txids)
            + std::mem::size_of_val(&self.gamma_vec);
        size
    }

    pub fn e_step(&mut self, txabds: &Vec<TxAbundance>) {
        let mut total = vec![0.0f32; self.sample_size];

        // pass 1: compute total[s]
        for &txid in &self.txids {
            let txabd = &txabds[txid];
            // let inv_pt = 1.0 / txabd.pt;
            for sid in 0..self.sample_size {
                total[sid] += txabd.abundance_cur[sid] * txabd.inv_pt;
            }
        }

        // pass 2: compute gamma
        for (tid, &txid) in self.txids.iter().enumerate() {
            let txabd = &txabds[txid];
            // let inv_pt = 1.0 / txabd.pt;
            let base = tid * self.sample_size;
            for sid in 0..self.sample_size {
                let denom = total[sid];
                let idx = base + sid;
                if denom > 0.0 {
                    let numer = txabd.abundance_cur[sid] * txabd.inv_pt;
                    self.gamma_vec[idx] = numer / denom;
                } else {
                    self.gamma_vec[idx] = 0.0;
                }
            }
        }
    }
}
