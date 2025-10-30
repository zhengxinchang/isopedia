// maintain an function that inks fragmented reads into query
// a reusable matrix structure, that records the row and column max.
//

use crate::isoformarchive::read_record_from_mmap;
use crate::tmpidx::MergedIsoformOffsetPtr;
use ahash::HashSet;
use memmap2::Mmap;
pub struct Assembler {
    pub splice_sites: Vec<u32>,
    pub matix: Vec<Vec<u32>>,
    pub sample_count: usize,
    pub splice_site_count: usize,
    pub record_buf: Vec<u8>,
}

impl Assembler {
    pub fn init(index_size: usize) -> Self {
        //          sj1 sj2 sj3 ...
        // sample 1
        // sample 2

        let mut matrix = Vec::new();
        // for _ in 0..200 {
        //     // 100 sjs
        //     matrix.push(vec![0u32; index_size]);
        // }

        // outer is sample inner is splice site
        for _ in 0..index_size {
            matrix.push(vec![0u32; 200]); // 200 sjs
        }

        Assembler {
            splice_sites: Vec::new(),
            matix: matrix, //outer is splice site, inner is sample
            sample_count: index_size,
            splice_site_count: 0,
            record_buf: Vec::new(),
        }
    }

    /// reset the assembler state for reusable
    pub fn reset(&mut self) {
        // reset the matrix demensions
        // self.sample_count = 0;
        self.splice_site_count = 0;
        // clear the matrix content
        for row in self.matix.iter_mut() {
            for val in row.iter_mut() {
                *val = 0;
            }
        }
    }

    fn add_at(&mut self, sample_idx: usize, splice_site_idx: usize, val: u32) {
        self.matix[sample_idx][splice_site_idx] = self.matix[sample_idx][splice_site_idx] + val;
    }

    fn get(&self, sample_idx: usize, splice_site_idx: usize) -> u32 {
        self.matix[sample_idx][splice_site_idx]
    }

    pub fn assemble(
        &mut self,
        queries: &Vec<(String, u64)>,
        all_res: &Vec<Vec<MergedIsoformOffsetPtr>>,
        archive: &Mmap,
        flank: u64,
    ) -> (bool, Vec<u32>) {
        let splice_junctions: Vec<u64> = queries.iter().map(|q| q.1).collect();
        // make Vec<(64,64)> for splice junciton
        let mut sj_pairs: Vec<(u64, u64)> = Vec::new();
        // let mut sj_pairs_set: std::collections::HashSet<(u64, u64), ahash::RandomState> = HashSet::default();
        for i in 0..splice_junctions.len() / 2 {
            sj_pairs.push((splice_junctions[2 * i], splice_junctions[2 * i + 1]));
            // sj_pairs_set.insert((splice_junctions[2*i], splice_junctions[2*i+1]));
        }

        // record the matrix size
        self.splice_site_count = sj_pairs.len();

        // for each res, check if it has at least two splice junctions.

        let mut candidates = Vec::new();
        let mut candidates_sig_set = HashSet::default();

        // dbg!("Total isoform candidates: {}", all_res.iter().map(|r| r.len()).sum::<usize>());
        // dbg!(all_res.len());

        for rec in all_res.iter() {
            for r in rec {
                if r.n_splice_sites >= 2 && r.n_splice_sites < queries.len() as u32 {
                    // only consider isoforms with at least 2 splice junctions

                    // only consider the isoforms not covering all splice junctions
                    if !candidates_sig_set.contains(&r.offset) {
                        candidates_sig_set.insert(r.offset);
                        candidates.push(r);
                    }
                }
            }
        }

        // dbg!(candidates.len());

        for cand_offset in candidates.iter() {
            // load the isoform from archive
            let merged_rec = read_record_from_mmap(&archive, cand_offset, &mut self.record_buf);

            let (is_all_matched, first_match, matched_count) =
                merged_rec.match_splice_junctions(&sj_pairs, flank);

            // println!("query sj: {:?}\n", sj_pairs);
            // println!("isofo sj: {:?}\n", merged_rec.splice_junctions_vec);

            if is_all_matched {
                // mark the matrix
                // println!("Matched isoform: {}, sample_vec {:?}, matched_count: {}, first_match_pos: {}\n", merged_rec.signature, merged_rec.get_sample_evidence_arr(), matched_count, first_match);
                for idx_col in first_match as usize..(first_match as usize + matched_count as usize)
                {
                    for (idx_row, sample_evidence) in
                        merged_rec.get_sample_evidence_arr().iter().enumerate()
                    {
                        if *sample_evidence > 0 {
                            self.add_at(idx_row, idx_col, *sample_evidence);
                        }
                    }
                }
            }
        }

        self.estimate_sample_cov_vec()
    }

    fn _estimate_cov(sample_vec: &Vec<u32>, total_sample_size: &usize) -> u32 {
        // as of now, simply return the mean coverage
        let mut total_cov: u32 = 0;
        let mut count: u32 = 0;
        // for cov in sample_vec.iter() {
        //     if *cov > 0 {
        //         total_cov += *cov;
        //         count += 1;
        //     }
        // }
        for sidx in 0..*total_sample_size {
            let cov = sample_vec[sidx];
            if cov > 0 {
                total_cov += cov;
                count += 1;
            }
        }

        // dbg!(count);
        // dbg!(sample_vec.len());

        if count > 0 {
            if count == *total_sample_size as u32 {
                total_cov / count
            } else {
                0
            }
        } else {
            0
        }
    }

    pub fn estimate_sample_cov_vec(&self) -> (bool, Vec<u32>) {
        let mut is_hit = false;
        let mut sample_cov_vec: Vec<u32> = Vec::new();
        for sample_idx in 0..self.sample_count {
            let x_cov = Assembler::_estimate_cov(&self.matix[sample_idx], &self.splice_site_count);
            if x_cov > 0 {
                is_hit = true;
            }
            sample_cov_vec.push(x_cov);
        }
        (is_hit, sample_cov_vec)
    }

    pub fn print_matrix(&self) {
        // dbg!(&self.matix);

        // dbg!(self.sample_count);
        // dbg!(self.splice_site_count);

        print!("SJ\t");
        for j in 0..self.sample_count {
            print!("Sample#{}\t", j);
        }
        println!();

        for i in 0..self.splice_site_count {
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
