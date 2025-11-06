// maintain an function that inks fragmented reads into query
// a reusable matrix structure, that records the row and column max.
//

use crate::cmd::isoform::AnnIsoCli;
use crate::isoformarchive::ArchiveCache;
use crate::runtime::Runtime;
use crate::tmpidx::MergedIsoformOffsetPtr;
use ahash::HashSet;
pub struct Assembler {
    pub matix: Vec<u32>, // flattened matrix
    pub sample_count: usize,
    pub splice_junction_count: usize, // each junciton is a node
    pub record_buf: Vec<u8>,
}

impl Assembler {
    pub fn init(n_sample: usize) -> Self {
        //          sj1 sj2 sj3 ...
        // sample 1
        // sample 2

        let matrix = Vec::new();

        Assembler {
            matix: matrix, // flattened matrix
            sample_count: n_sample,
            splice_junction_count: 0,
            record_buf: Vec::new(),
        }
    }

    /// reset the assembler state for reusable
    pub fn reset(&mut self) {
        // reset the matrix demensions
        // self.sample_count = 0;
        self.splice_junction_count = 0;
        // clear the matrix content
        self.matix.clear();
    }

    fn add_at(&mut self, sample_idx: usize, splice_site_idx: usize, val: u32) {
        // expand if necessary
        let required_size = sample_idx * self.splice_junction_count + splice_site_idx;
        if self.matix.len() <= required_size {
            self.matix.resize(required_size + 1, 0); // initialize new elements to 0
        }

        self.matix[sample_idx * self.splice_junction_count + splice_site_idx] =
            self.matix[sample_idx * self.splice_junction_count + splice_site_idx] + val;
    }

    fn get(&self, sample_idx: usize, splice_site_idx: usize) -> u32 {
        let base_idx = sample_idx * self.splice_junction_count + splice_site_idx;
        if base_idx >= self.matix.len() {
            return 0;
        }

        self.matix[base_idx]
    }

    pub fn assemble(
        &mut self,
        queries: &Vec<(String, u64)>,
        all_res: &Vec<Vec<MergedIsoformOffsetPtr>>,
        archive: &mut ArchiveCache,
        cli: &AnnIsoCli,
        runtime: &mut Runtime,
    ) -> bool {
        let splice_junctions: Vec<u64> = queries.iter().map(|q| q.1).collect();
        // make Vec<(64,64)> for splice junciton
        let mut sj_pairs: Vec<(u64, u64)> = Vec::new();
        // let mut sj_pairs_set: std::collections::HashSet<(u64, u64), ahash::RandomState> = HashSet::default();
        for i in 0..splice_junctions.len() / 2 {
            sj_pairs.push((splice_junctions[2 * i], splice_junctions[2 * i + 1]));
            // sj_pairs_set.insert((splice_junctions[2*i], splice_junctions[2*i+1]));
        }

        // record the matrix size
        self.splice_junction_count = sj_pairs.len();

        // for each res, check if it has at least two splice junctions.

        let mut candidates_offsets = Vec::new();
        let mut candidates_offsets_dedup = HashSet::default();

        for offsets in all_res.iter() {
            for offset in offsets {
                if offset.n_splice_sites >= 2 && offset.n_splice_sites < queries.len() as u32 {
                    if !candidates_offsets_dedup.contains(&offset.offset) {
                        candidates_offsets_dedup.insert(offset.offset);
                        candidates_offsets.push(offset);
                    }
                }
            }
        }

        for cand_offset in candidates_offsets.iter() {
            // load the isoform from archive
            let record = archive.read_bytes(*cand_offset);

            let (is_all_matched, first_match, matched_count) =
                record.match_splice_junctions(&sj_pairs, cli.flank);

            if is_all_matched {
                if cli.info {
                    runtime
                        .ism_trigger_add_read_info(&record)
                        .unwrap_or_else(|_| panic!("can not add read info for ism reads"));
                }

                for idx_col in first_match as usize..(first_match as usize + matched_count as usize)
                {
                    for (idx_row, sample_evidence) in
                        record.get_sample_evidence_arr().iter().enumerate()
                    {
                        if *sample_evidence > 0 {
                            self.add_at(idx_row, idx_col, *sample_evidence);
                        }
                    }
                }
            }
        }

        let (is_hit, cov_vec) = self.estimate_sample_cov_vec();

        if is_hit {
            runtime
                .update_ism_record(&cov_vec)
                .unwrap_or_else(|_| panic!("can not update ism record"));
        }

        // self.print_matrix();

        is_hit
    }

    fn _estimate_cov(sample_vec: &Vec<u32>) -> u32 {
        // as of now, simply return the mean coverage
        let mut total_cov: u32 = 0;
        let mut count: u32 = 0;
        for sidx in 0..sample_vec.len() {
            let cov = sample_vec[sidx];
            if cov > 0 {
                total_cov += cov;
                count += 1;
            }
        }

        if count > 0 {
            if count == sample_vec.len() as u32 {
                total_cov / count
            } else {
                0
            }
        } else {
            0
        }
    }

    pub fn get_sample_vec(&self, sample_idx: usize) -> Option<Vec<u32>> {
        let base_idx = sample_idx * self.splice_junction_count;

        if base_idx >= self.matix.len() {
            return None;
        }
        if base_idx + self.splice_junction_count > self.matix.len() {
            return None;
        }

        Some(self.matix[base_idx..(base_idx + self.splice_junction_count)].to_vec())
    }

    pub fn estimate_sample_cov_vec(&self) -> (bool, Vec<u32>) {
        let mut is_hit = false;
        let mut sample_cov_vec: Vec<u32> = Vec::new();
        for sample_idx in 0..self.sample_count {
            match self.get_sample_vec(sample_idx) {
                Some(s) => {
                    let x_cov = Assembler::_estimate_cov(&s);
                    if x_cov > 0 {
                        is_hit = true;
                    }
                    sample_cov_vec.push(x_cov);
                }
                None => {
                    sample_cov_vec.push(0);
                }
            }
        }
        (is_hit, sample_cov_vec)
    }

    pub fn print_matrix(&self) {
        print!("SJ\t");
        for j in 0..self.sample_count {
            print!("Sample#{}\t", j);
        }
        println!();

        for i in 0..self.splice_junction_count {
            print!("SJ#{}\t", i);
            for j in 0..self.sample_count {
                print!("{}\t", self.get(j, i));
            }
            println!();
        }

        let mut count = 0;
        self.estimate_sample_cov_vec().1.iter().for_each(|cov| {
            count += 1;
            println!("Sample#{} estimated cov: {}", count, cov);
        });
    }
}
