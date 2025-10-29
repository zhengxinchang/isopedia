// maintain an function that inks fragmented reads into query
// a reusable matrix structure, that records the row and column max.
//

use crate::tmpidx::MergedIsoformOffsetPtr;

pub struct Assembler {
    pub splice_sites: Vec<u32>,
    pub matix: Vec<Vec<u8>>,
    pub sample_count: usize,
    pub splice_site_count: usize,
}

impl Assembler {
    pub fn init(index_size: usize) -> Self {
        //          sj1 sj2 sj3 ...
        // sample 1
        // sample 2

        let mut matrix = Vec::new();
        for _ in 0..100 {
            matrix.push(vec![0u8; index_size]);
        }

        Assembler {
            splice_sites: Vec::new(),
            matix: matrix,
            sample_count: 0,
            splice_site_count: 0,
        }
    }

    /// reset the assembler state for reusable
    pub fn reset(&mut self) {}

    pub fn assemble(
        &mut self,
        queries: &Vec<(usize, u32)>,
        all_res: &Vec<Vec<MergedIsoformOffsetPtr>>,
    ) -> Vec<(usize, u32)> {
        todo!()
    }
}
