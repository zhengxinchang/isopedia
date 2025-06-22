use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
pub struct ChromMapping {
    pub chrom_map: IndexMap<String, u16>,
    pub idx_map: IndexMap<u16, String>,
}

// #[allow(unused)]
impl ChromMapping {
    pub fn new() -> ChromMapping {
        ChromMapping {
            chrom_map: IndexMap::new(),
            idx_map: IndexMap::new(),
        }
    }

    pub fn add_chrom(&mut self, chrom: String) -> u16 {
        if let Some(idx) = self.chrom_map.get(&chrom) {
            *idx
        } else {
            let idx = self.idx_map.len() as u16;
            self.chrom_map.insert(chrom.clone(), idx);
            self.idx_map.insert(idx, chrom);
            idx
        }
    }

    pub fn get_chrom_name(&self, idx: u16) -> &str {
        self.idx_map.get(&idx).unwrap()
    }

    pub fn get_chrom_idx(&self, name: &str) -> Option<u16> {
        self.chrom_map.get(name).map(|idx| *idx)
    }

    pub fn get_size(&self) -> usize {
        self.chrom_map.len()
    }

    pub fn get_chrom_names(&self) -> Vec<String> {
        self.idx_map.values().cloned().collect()
    }

    pub fn get_chrom_idxs(&self) -> Vec<u16> {
        self.chrom_map.values().cloned().collect()
    }

    pub fn sort_by_chrom_name(&mut self) {
        let mut sorted_chroms = self.get_chrom_names();
        sorted_chroms.sort();
        self.chrom_map = IndexMap::new();
        self.idx_map = IndexMap::new();
        for (idx, chrom) in sorted_chroms.into_iter().enumerate() {
            {
                self.chrom_map.insert(chrom.clone(), idx as u16);
                self.idx_map.insert(idx as u16, chrom);
            }
        }
    }

    pub fn udpate_from_string(&mut self, string: &str) {
        for chrom in string.trim().split('\t').filter(|x| !x.is_empty()) {
            self.add_chrom(chrom.to_string());
        }
    }

    pub fn encode(&self) -> Vec<u8> {
        bincode::serialize(&self).unwrap()
    }

    pub fn decode(bytes: &[u8]) -> ChromMapping {
        bincode::deserialize(bytes).unwrap()
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_sort_by_chrom_name() {
        // create a Chromosomes object
        let mut chroms = ChromMapping::new();
        // add some chromosomes
        chroms.add_chrom("chr1".to_string());
        chroms.add_chrom("chr2".to_string());
        chroms.add_chrom("chr4".to_string());
        chroms.add_chrom("chr5".to_string());
        chroms.add_chrom("chr6".to_string());
        chroms.add_chrom("chr3".to_string());
        chroms.add_chrom("chr8".to_string());
        chroms.add_chrom("chr9".to_string());
        chroms.add_chrom("chr10".to_string());
        chroms.add_chrom("chr7".to_string());
        // sort the chromosomes
        chroms.sort_by_chrom_name();
        // get the sorted chromosome names
        let chrom_names = chroms.get_chrom_names();
        dbg!(&chrom_names);
        dbg!(chroms.get_chrom_idx("chr10"));
    }

    #[test]
    fn test_encode_decode() {
        // create a Chromosomes object
        let mut chroms = ChromMapping::new();
        // add some chromosomes
        chroms.add_chrom("chr1".to_string());
        chroms.add_chrom("chr2".to_string());
        chroms.add_chrom("chr4".to_string());
        chroms.add_chrom("chr5".to_string());
        chroms.add_chrom("chr6".to_string());
        chroms.add_chrom("chr3".to_string());
        chroms.add_chrom("chr8".to_string());
        chroms.add_chrom("chr9".to_string());
        chroms.add_chrom("chr10".to_string());
        chroms.add_chrom("chr7".to_string());
        // encode the Chromosomes object
        let bytes = chroms.encode();
        // decode the bytes
        let chroms2 = ChromMapping::decode(&bytes);
        // get the chromosome names
        let chrom_names = chroms.get_chrom_names();
        let chrom_names2 = chroms2.get_chrom_names();
        // check if the chromosome names are the same
        assert_eq!(chrom_names, chrom_names2);
    }
}
