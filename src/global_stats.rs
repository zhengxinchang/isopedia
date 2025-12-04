pub struct GlobalStats {
    pub fsm_total: Vec<f32>,    // totoal fsm transcripts per sample
    pub em_total: Vec<f32>,     // total ism transcripts per sample
    pub fsm_em_total: Vec<f32>, // total fsm+ism transcripts per sample
}

impl GlobalStats {
    pub fn new(n_sample: usize) -> Self {
        GlobalStats {
            fsm_total: vec![0.0; n_sample],
            em_total: vec![0.0; n_sample],
            fsm_em_total: vec![0.0; n_sample],
        }
    }

    pub fn update_fsm_total(&mut self, fsm_vec: &Vec<f32>) {
        for (i, count) in fsm_vec.iter().enumerate() {
            self.fsm_total[i] += *count;
            self.fsm_em_total[i] += *count;
        }
    }

    pub fn update_em_total(&mut self, em_vec: &Vec<f32>) {
        for (i, count) in em_vec.iter().enumerate() {
            self.em_total[i] += *count;
            self.fsm_em_total[i] += *count;
        }
    }

    pub fn get_positive_samples(&self, min_read: u32) -> usize {
        let mut count = 0;
        for i in 0..self.fsm_em_total.len() {
            if self.fsm_em_total[i] >= min_read as f32 {
                count += 1;
            }
        }
        count
    }
}
