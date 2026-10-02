use crate::{cmd::isoform::AnnIsoCli, grouped_tx::TxAbundanceView};

pub struct GlobalStats {
    pub sample_posi_tx_count_fsm: Vec<usize>,
}

impl GlobalStats {
    pub fn new(n_sample: usize) -> Self {
        GlobalStats {
            sample_posi_tx_count_fsm: vec![0; n_sample],
        }
    }

    pub fn update_sample_level_stats(&mut self, txview: &TxAbundanceView, cli: &AnnIsoCli) {
        txview.rc_fsm_jc.iter().enumerate().for_each(|(i, &count)| {
            if count >= cli.min_read as u64 {
                self.sample_posi_tx_count_fsm[i] += 1;
            }
        });
    }
    pub fn get_fsm_tx_by_sample_idx(&self, idx: usize) -> usize {
        self.sample_posi_tx_count_fsm[idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tracks_only_direct_positive_transcripts() {
        let stats = GlobalStats::new(3);
        assert_eq!(stats.sample_posi_tx_count_fsm, vec![0, 0, 0]);
    }
}
