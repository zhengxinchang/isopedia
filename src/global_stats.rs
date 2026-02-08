use crate::{cmd::isoform::AnnIsoCli, grouped_tx::TxAbundanceView};

pub struct GlobalStats {
    pub fsm_tx_abd_total: Vec<f32>,    // totoal fsm transcripts per sample
    pub em_tx_abd_total: Vec<f32>,     // total ism transcripts per sample
    pub fsm_em_tx_abd_total: Vec<f32>, // total fsm+ism transcripts per sample
    pub sample_posi_tx_count_fsm: Vec<usize>, // size == sample size, each element shows how may fsm postive transcripts in the sample
    pub sample_posi_tx_count_em: Vec<usize>, // size == sample size, each element shows how may ism postive transcripts in the sample
    pub sample_posi_tx_count_fsm_em: Vec<usize>, // size == sample size, each element shows how may fsm+ism postive transcripts in the sample
}

impl GlobalStats {
    pub fn new(n_sample: usize) -> Self {
        GlobalStats {
            fsm_tx_abd_total: vec![0.0; n_sample],
            em_tx_abd_total: vec![0.0; n_sample],
            fsm_em_tx_abd_total: vec![0.0; n_sample],
            sample_posi_tx_count_fsm: vec![0; n_sample],
            sample_posi_tx_count_em: vec![0; n_sample],
            sample_posi_tx_count_fsm_em: vec![0; n_sample],
        }
    }

    pub fn update_fsm_tx_abd_total(&mut self, fsm_vec: &Vec<f32>) {
        for (i, count) in fsm_vec.iter().enumerate() {
            self.fsm_tx_abd_total[i] += *count;
            self.fsm_em_tx_abd_total[i] += *count;
        }
    }

    pub fn update_em_tx_abd_total(&mut self, em_vec: &Vec<f32>) {
        for (i, count) in em_vec.iter().enumerate() {
            self.em_tx_abd_total[i] += *count;
            self.fsm_em_tx_abd_total[i] += *count;
        }
    }

    pub fn update_sample_level_stats(&mut self, txview: &TxAbundanceView, cli: &AnnIsoCli) {
        txview
            .fsm_abundance
            .iter()
            .enumerate()
            .for_each(|(i, abd)| {
                if *abd >= cli.min_read as f32 {
                    self.sample_posi_tx_count_fsm[i] += 1;
                }

                if txview.em_abundance[i] >= cli.min_read as f32 {
                    self.sample_posi_tx_count_em[i] += 1;
                }

                if *abd >= cli.min_read as f32 || txview.em_abundance[i] >= cli.min_read as f32 {
                    self.sample_posi_tx_count_fsm_em[i] += 1;
                }
            });
    }
    pub fn get_fsm_tx_by_sample_idx(&self, idx: usize) -> usize {
        self.sample_posi_tx_count_fsm_em[idx]
    }
    pub fn get_em_tx_by_sample_idx(&self, idx: usize) -> usize {
        self.sample_posi_tx_count_em[idx]
    }
    pub fn get_fsm_em_tx_by_sample_idx(&self, idx: usize) -> usize {
        self.sample_posi_tx_count_fsm_em[idx]
    }
}
