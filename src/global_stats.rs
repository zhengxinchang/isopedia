pub struct GlobalStats {
    pub fsm_total: Vec<f32>,     // totoal fsm transcripts per sample
    pub ism_total: Vec<f32>,     // total ism transcripts per sample
    pub fsm_ism_total: Vec<f32>, // total fsm+ism transcripts per sample
}

impl GlobalStats {
    pub fn new(n_sample: usize) -> Self {
        GlobalStats {
            fsm_total: vec![0.0; n_sample],
            ism_total: vec![0.0; n_sample],
            fsm_ism_total: vec![0.0; n_sample],
        }
    }

    pub fn update_fsm_total(&mut self, fsm_vec: &Vec<f32>) {
        for (i, count) in fsm_vec.iter().enumerate() {
            self.fsm_total[i] += *count;
            self.fsm_ism_total[i] += *count;
        }
    }

    pub fn update_ism_total(&mut self, ism_vec: &Vec<f32>) {
        for (i, count) in ism_vec.iter().enumerate() {
            self.ism_total[i] += *count;
            self.fsm_ism_total[i] += *count;
        }
    }
}

// use crate::{cmd::isoform::AnnIsoCli, dataset_info::DatasetInfo, isoform::MergedIsoform, utils};
// use anyhow::Result;

// pub struct Runtime {
//     n_sample: usize,
//     pub total: u64,
//     pub fsm_hit: u64,
//     pub ism_hit: u64,
//     pub missed_hit: u64,
//     pub sample_wide_positive_transcript_fsm: Vec<u64>,
//     pub sample_wide_positive_transcript_ism: Vec<u64>,
//     pub sample_wide_positive_transcript_both: Vec<u64>,
//     fsm_acc_positive_sample_vec_by_min_read: Vec<usize>, // positive sample count that defined by min-read, aggregated multiple record
//     fsm_acc_sample_read_info_vec: Vec<String>, // per-sample info string, aggregated multiple record
//     fsm_acc_sample_evidence_arr: Vec<u32>, // per-sample evidence strings, aggregated multiple recor
//     ism_acc_positive_sample_vec_by_min_read: Vec<usize>, // positive sample count that defined by min-read, aggregated multiple record
//     ism_acc_sample_read_info_vec: Vec<String>, // per-sample info string, aggregated multiple record
//     ism_acc_sample_evidence_arr: Vec<u32>, // per-sample evidence strings, aggregated multiple record
// }

// impl Runtime {
//     pub fn init(n_sample: usize) -> Self {
//         Runtime {
//             n_sample,
//             total: 0,
//             fsm_acc_positive_sample_vec_by_min_read: vec![0; n_sample],
//             fsm_acc_sample_read_info_vec: vec!["NULL".to_string(); n_sample],
//             fsm_acc_sample_evidence_arr: vec![0; n_sample],
//             ism_acc_positive_sample_vec_by_min_read: vec![0; n_sample],
//             ism_acc_sample_read_info_vec: vec!["NULL".to_string(); n_sample],
//             ism_acc_sample_evidence_arr: vec![0; n_sample],
//             fsm_hit: 0,
//             ism_hit: 0,
//             missed_hit: 0,
//             sample_wide_positive_transcript_fsm: vec![0; n_sample],
//             sample_wide_positive_transcript_ism: vec![0; n_sample],
//             sample_wide_positive_transcript_both: vec![0; n_sample],
//         }
//     }

//     pub fn update_fsm_record(&mut self, record: &MergedIsoform) -> Result<()> {
//         for (sample_idx, evidence) in record.sample_evidence_arr.iter().enumerate() {
//             self.fsm_acc_sample_evidence_arr[sample_idx] += evidence;
//             // self.positive_transcript_count_by_sample_fsm[sample_idx] += 1;
//         }

//         Ok(())
//     }

//     pub fn update_ism_record(&mut self, cov_vec: &Vec<u32>) -> Result<()> {
//         if self.n_sample != cov_vec.len() {
//             return Err(anyhow::anyhow!(
//                 "Mismatch sample size when updating ism record"
//             ));
//         }

//         for (sample_idx, evidence) in cov_vec.iter().enumerate() {
//             self.ism_acc_sample_evidence_arr[sample_idx] += evidence;
//         }

//         Ok(())
//     }

//     pub fn add_one_fsm_hit(&mut self) {
//         self.fsm_hit += 1;
//     }

//     pub fn add_one_ism_hit(&mut self) {
//         self.ism_hit += 1;
//     }

//     pub fn add_one_missed_hit(&mut self) {
//         self.missed_hit += 1;
//     }

//     pub fn add_one_total(&mut self) {
//         self.total += 1;
//     }

//     pub fn get_fsm_positive_count_by_min_read(&self, cli: &AnnIsoCli) -> usize {
//         let mut count = 0;
//         for evidence in self.fsm_acc_sample_evidence_arr.iter() {
//             if *evidence >= cli.min_read {
//                 count += 1;
//             }
//         }
//         count
//     }

//     pub fn get_fsm_ism_positive_count_by_min_read(&self, cli: &AnnIsoCli) -> usize {
//         let mut count = 0;
//         for i in 0..self.n_sample {
//             let combined_evidence =
//                 self.fsm_acc_sample_evidence_arr[i] + self.ism_acc_sample_evidence_arr[i];
//             if combined_evidence >= cli.min_read {
//                 count += 1;
//             }
//         }
//         count
//     }

//     pub fn get_fsm_confidence(&self, dataset_info: &DatasetInfo) -> f64 {
//         MergedIsoform::get_confidence_value(
//             &self.fsm_acc_sample_evidence_arr,
//             dataset_info.get_size(),
//             &dataset_info.sample_total_evidence_vec,
//         )
//     }

//     pub fn get_fsm_ism_confidence(&self, dataset_info: &DatasetInfo) -> f64 {
//         let mut combined_arr: Vec<u32> = vec![0; self.n_sample];
//         for i in 0..self.n_sample {
//             combined_arr[i] =
//                 self.fsm_acc_sample_evidence_arr[i] + self.ism_acc_sample_evidence_arr[i];
//         }

//         MergedIsoform::get_confidence_value(
//             &combined_arr,
//             dataset_info.get_size(),
//             &dataset_info.sample_total_evidence_vec,
//         )
//     }

//     pub fn get_fsm_read_info_vec(&self) -> &Vec<String> {
//         &self.fsm_acc_sample_read_info_vec
//     }

//     pub fn get_ism_read_info_vec(&self) -> &Vec<String> {
//         &self.ism_acc_sample_read_info_vec
//     }

//     pub fn get_fsm_ism_read_info_vec(&self) -> Vec<String> {
//         let mut combined_vec: Vec<String> = vec!["NULL".to_string(); self.n_sample];
//         for i in 0..self.n_sample {
//             if self.fsm_acc_sample_read_info_vec[i] == "NULL"
//                 && self.ism_acc_sample_read_info_vec[i] == "NULL"
//             {
//                 combined_vec[i] = "NULL".to_string();
//             } else if self.fsm_acc_sample_read_info_vec[i] == "NULL" {
//                 combined_vec[i] = self.ism_acc_sample_read_info_vec[i].clone();
//             } else if self.ism_acc_sample_read_info_vec[i] == "NULL" {
//                 combined_vec[i] = self.fsm_acc_sample_read_info_vec[i].clone();
//             } else {
//                 combined_vec[i] = format!(
//                     "{},{}",
//                     self.fsm_acc_sample_read_info_vec[i], self.ism_acc_sample_read_info_vec[i]
//                 );
//             }
//         }
//         combined_vec
//     }

//     pub fn get_sample_chip_data(
//         &self,
//         dataset_info: &DatasetInfo,
//         cli: &AnnIsoCli,
//     ) -> Vec<(f64, u32, String)> {
//         let mut cov_vec = Vec::with_capacity(self.n_sample);
//         let mut cpm_vec = Vec::with_capacity(self.n_sample);
//         let mut read_info_vec = Vec::with_capacity(self.n_sample);

//         match cli.asm {
//             true => {
//                 let read_info = self.get_fsm_ism_read_info_vec();

//                 for i in 0..self.n_sample {
//                     let combined_evidence =
//                         self.fsm_acc_sample_evidence_arr[i] + self.ism_acc_sample_evidence_arr[i];
//                     let cpm = utils::calc_cpm(
//                         &combined_evidence,
//                         &dataset_info.sample_total_evidence_vec[i],
//                     );
//                     cov_vec.push(combined_evidence);
//                     cpm_vec.push(cpm);
//                     // combine read info
//                     let mut info = "TIER=FSM+ISM,".to_string();
//                     match cli.info {
//                         true => {
//                             info.push_str(&read_info[i]);
//                             read_info_vec.push(info);
//                         }
//                         false => {
//                             read_info_vec.push("NULL".to_string());
//                         }
//                     }
//                 }
//             }
//             false => {
//                 let read_info = self.get_fsm_read_info_vec();

//                 for i in 0..self.n_sample {
//                     let evidence = self.fsm_acc_sample_evidence_arr[i];
//                     let cpm =
//                         utils::calc_cpm(&evidence, &dataset_info.sample_total_evidence_vec[i]);
//                     cov_vec.push(evidence);
//                     cpm_vec.push(cpm);
//                     // read info
//                     let mut info = "TIER=FSM,".to_string();
//                     match cli.info {
//                         true => {
//                             info.push_str(&read_info[i]);
//                             read_info_vec.push(info);
//                         }
//                         false => {
//                             read_info_vec.push("NULL".to_string());
//                         }
//                     }
//                 }
//             }
//         }

//         let mut result_vec: Vec<(f64, u32, String)> = Vec::with_capacity(self.n_sample);
//         for i in 0..self.n_sample {
//             result_vec.push((cpm_vec[i], cov_vec[i], read_info_vec[i].clone()));
//         }

//         result_vec
//     }

//     pub fn fsm_trigger_add_read_info(&mut self, record: &MergedIsoform) -> Result<()> {
//         let single_sample_evidence_arr = record.get_sample_evidence_arr();
//         for (i, ofs) in record.get_sample_offset_arr().iter().enumerate() {
//             let ofs = *ofs as usize;
//             let length = single_sample_evidence_arr[i] as usize;
//             // get readdiffslim
//             if length > 0 {
//                 let read_info_str = record.isoform_reads_slim_vec[ofs..ofs + length]
//                     .iter()
//                     .map(|delta| delta.to_string_no_offsets())
//                     .collect::<Vec<String>>()
//                     .join(",");
//                 if self.fsm_acc_sample_read_info_vec[i] == "NULL" {
//                     self.fsm_acc_sample_read_info_vec[i] = read_info_str;
//                 } else {
//                     self.fsm_acc_sample_read_info_vec[i].push(',');
//                     self.fsm_acc_sample_read_info_vec[i].push_str(&read_info_str);
//                     // format!("{},{}", self.acc_sample_read_info_vec[i], read_info_str);
//                 }
//             }
//         }
//         Ok(())
//     }

//     pub fn ism_trigger_add_read_info(&mut self, record: &MergedIsoform) -> Result<()> {
//         let single_sample_evidence_arr = record.get_sample_evidence_arr();
//         for (i, ofs) in record.get_sample_offset_arr().iter().enumerate() {
//             let ofs = *ofs as usize;
//             let length = single_sample_evidence_arr[i] as usize;
//             // get readdiffslim
//             if length > 0 {
//                 let read_info_str = record.isoform_reads_slim_vec[ofs..ofs + length]
//                     .iter()
//                     .map(|delta| delta.to_string_no_offsets())
//                     .collect::<Vec<String>>()
//                     .join(",");
//                 if self.ism_acc_sample_read_info_vec[i] == "NULL" {
//                     self.ism_acc_sample_read_info_vec[i] = read_info_str;
//                 } else {
//                     self.ism_acc_sample_read_info_vec[i].push(',');
//                     self.ism_acc_sample_read_info_vec[i].push_str(&read_info_str);
//                     // format!("{},{}", self.acc_sample_read_info_vec[i], read_info_str);
//                 }
//             }
//         }
//         Ok(())
//     }

//     pub fn log_stats(&mut self) {
//         // 这里是检验每个query是否在每个样本中hit了，如果有，则添加到对应的
//         // sample_wide* 中
//         // 因为有的query，既可以是FSM，也可以是ISM，所以 某个样本的 FSM hit + ISM hit
//         // 可能会大于1
//         // 这里需要注意鉴别的是 FSM 和 ISM 是分别统计的

//         for sample_idx in 0..self.n_sample {
//             if self.fsm_acc_sample_evidence_arr[sample_idx] > 0 {
//                 self.sample_wide_positive_transcript_fsm[sample_idx] += 1;
//             }
//             if self.ism_acc_sample_evidence_arr[sample_idx] > 0 {
//                 self.sample_wide_positive_transcript_ism[sample_idx] += 1;
//             }

//             if self.fsm_acc_sample_evidence_arr[sample_idx] > 0
//                 || self.ism_acc_sample_evidence_arr[sample_idx] > 0
//             {
//                 self.sample_wide_positive_transcript_both[sample_idx] += 1;
//             }
//         }
//     }
//     pub fn reset(&mut self) {
//         self.fsm_acc_positive_sample_vec_by_min_read = vec![0; self.n_sample];
//         self.ism_acc_positive_sample_vec_by_min_read = vec![0; self.n_sample];
//         self.fsm_acc_sample_read_info_vec = vec!["NULL".to_string(); self.n_sample];
//         self.ism_acc_sample_read_info_vec = vec!["NULL".to_string(); self.n_sample];
//         self.fsm_acc_sample_evidence_arr = vec![0; self.n_sample];
//         self.ism_acc_sample_evidence_arr = vec![0; self.n_sample];
//     }
// }
