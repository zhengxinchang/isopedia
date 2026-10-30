use indexmap::IndexMap;
use log::info;
use serde::{Deserialize, Serialize};

use crate::utils::trim_chr_prefix_to_upper;

#[derive(Debug, Deserialize, Serialize)]
pub struct ChromMapping {
    pub name2id: IndexMap<String, u16>,
    pub id2name: IndexMap<u16, String>,
}

impl ChromMapping {
    pub fn new() -> ChromMapping {
        ChromMapping {
            name2id: IndexMap::new(),
            id2name: IndexMap::new(),
        }
    }

    pub fn add_chrom(&mut self, chrom: String) -> u16 {
        let normalized = trim_chr_prefix_to_upper(&chrom);

        if let Some(idx) = self.name2id.get(&normalized) {
            *idx
        } else {
            let idx = self.id2name.len() as u16;
            self.name2id.insert(normalized.clone(), idx);
            self.id2name.insert(idx, normalized);
            idx
        }
    }

    pub fn get_chrom_name(&self, idx: u16) -> &str {
        self.id2name.get(&idx).unwrap()
    }

    pub fn get_chrom_idx(&self, name: &str) -> Option<u16> {
        if let Some(idx) = self.name2id.get(name) {
            return Some(*idx);
        }

        let normalized = trim_chr_prefix_to_upper(name);
        if let Some(idx) = self.name2id.get(&normalized) {
            return Some(*idx);
        }

        self.name2id.iter().find_map(|(stored_name, idx)| {
            if trim_chr_prefix_to_upper(stored_name) == normalized {
                Some(*idx)
            } else {
                None
            }
        })
    }

    pub fn get_size(&self) -> usize {
        self.name2id.len()
    }

    pub fn get_chrom_names(&self) -> Vec<String> {
        self.id2name.values().cloned().collect()
    }

    pub fn get_chrom_idxs(&self) -> Vec<u16> {
        self.name2id.values().cloned().collect()
    }

    pub fn sort_by_chrom_name(&mut self) {
        let mut sorted_chroms = self.get_chrom_names();
        sorted_chroms.sort();
        self.name2id = IndexMap::new();
        self.id2name = IndexMap::new();
        for (idx, chrom) in sorted_chroms.into_iter().enumerate() {
            {
                self.name2id.insert(chrom.clone(), idx as u16);
                self.id2name.insert(idx as u16, chrom);
            }
        }
    }

    pub fn update_from_string(&mut self, string: &str) {
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
mod tests {
    use super::ChromMapping;

    #[test]
    fn chrom_mapping_normalizes_on_insert_and_lookup() {
        let mut mapping = ChromMapping::new();
        let idx = mapping.add_chrom("Chr1".to_string());

        assert_eq!(mapping.get_chrom_idx("Chr1"), Some(idx));
        assert_eq!(mapping.get_chrom_idx("chr1"), Some(idx));
        assert_eq!(mapping.get_chrom_idx("1"), Some(idx));
        assert_eq!(mapping.get_chrom_name(idx), "1");
    }
}

pub struct ChromMappingHelper {
    pub id2rec_count: IndexMap<u16, u64>,
}

impl ChromMappingHelper {
    pub fn new() -> ChromMappingHelper {
        let id2rec_count = IndexMap::new();
        ChromMappingHelper { id2rec_count }
    }

    pub fn add_record_to_chrom(&mut self, chrom_idx: u16) {
        if let Some(record) = self.id2rec_count.get_mut(&chrom_idx) {
            *record += 1;
        } else {
            self.id2rec_count.insert(chrom_idx, 1);
        }
    }

    pub fn get_record_count(&self, idx: u16) -> u64 {
        *self.id2rec_count.get(&idx).unwrap()
    }

    pub fn drop_zero(&mut self, chrmap: &mut ChromMapping) {
        let valid_chrom_id = self.id2rec_count.keys().cloned().collect::<Vec<u16>>();
        let mut new_name2id = IndexMap::new();
        let mut new_id2name = IndexMap::new();

        info!(
            "Drop {} chromosomes with zero records",
            chrmap.name2id.len() - valid_chrom_id.len()
        );

        let all_chrom_ids = chrmap.get_chrom_idxs();
        for chrom_id in all_chrom_ids {
            if !valid_chrom_id.contains(&chrom_id) {
                info!(
                    "Drop chromosome {} with zero records",
                    chrmap.get_chrom_name(chrom_id)
                );
            }
        }

        for chrom_id in self.id2rec_count.keys() {
            let chrom_name = chrmap.get_chrom_name(*chrom_id).to_string();
            new_name2id.insert(chrom_name.clone(), *chrom_id);
            new_id2name.insert(*chrom_id, chrom_name);
        }

        chrmap.name2id = new_name2id;
        chrmap.id2name = new_id2name;
    }
}
