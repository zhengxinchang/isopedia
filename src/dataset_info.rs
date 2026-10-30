use anyhow::Result;
use log::warn;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
    vec,
};

use crate::utils::split_line_fields;

#[derive(Serialize, Deserialize, Debug)]

pub struct DatasetInfo {
    name2idx: FxHashMap<String, u32>,
    name_vec: Vec<String>,
    path_vec: Vec<PathBuf>,
    pub sample_total_evidence_vec: Vec<u32>,
    pub sample_size: usize,
}

impl DatasetInfo {
    pub fn new() -> DatasetInfo {
        DatasetInfo {
            name2idx: FxHashMap::default(), //dont use it directly as its not ordered
            name_vec: Vec::new(),
            path_vec: Vec::new(),
            sample_total_evidence_vec: Vec::new(),
            sample_size: 0,
        }
    }

    /// parse sample from tab-separated file
    pub fn parse_manifest(path: &PathBuf) -> Result<DatasetInfo> {
        let content = std::fs::read_to_string(path)?;
        let mut dbinfo = DatasetInfo::new();
        for (i, line) in content.lines().skip(1).enumerate() {
            let fields = split_line_fields(line);

            if fields.is_empty() {
                warn!("Skipping empty line {} in manifest", i + 1);
                continue;
            }

            if fields.len() < 2 {
                return Err(anyhow::anyhow!("Invalid line in manifest: {}", line));
            }
            let name = fields[0].clone();
            let path = PathBuf::from(&fields[1]);
            dbinfo.add_sample(name, Some(path));
        }
        dbinfo.sample_total_evidence_vec = vec![0; dbinfo.sample_size];
        Ok(dbinfo)
    }

    pub fn save_to_file<P: AsRef<Path>>(self, path: P) -> std::io::Result<()> {
        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        let mut file = File::create(path)?;
        file.write_all(&encoded)?;
        Ok(())
    }

    pub fn load_from_file<P: AsRef<Path>>(path: P) -> std::io::Result<DatasetInfo> {
        let bytes: Vec<u8> = fs::read(path.as_ref())?;
        let decoded: DatasetInfo = bincode::deserialize(&bytes).unwrap();
        Ok(decoded)
    }

    pub fn add_sample(&mut self, name: String, file_path: Option<PathBuf>) -> u32 {
        if let Some(idx) = self.name2idx.get(&name) {
            *idx
        } else {
            let idx = self.name_vec.len() as u32;
            self.name2idx.insert(name.clone(), idx);
            self.name_vec.push(name);

            if let Some(f) = file_path {
                self.path_vec.push(f);
            }
            self.sample_size += 1;
            idx
        }
    }

    pub fn add_sample_evidence(&mut self, idx: usize, evidence: u32) {
        self.sample_total_evidence_vec[idx] += evidence;
    }

    pub fn get_path_list(&self) -> Vec<PathBuf> {
        self.path_vec.clone()
    }

    pub fn get_size(&self) -> usize {
        self.sample_size
    }

    pub fn get_sample_names(&self) -> Vec<String> {
        self.name_vec.clone()
    }

    pub fn get_sample_evidence_pair_vec(&self) -> Vec<(String, u32)> {
        let mut v = Vec::new();
        for name in self.name_vec.iter() {
            let idx = self.name2idx.get(name).unwrap();
            v.push((name.clone(), self.sample_total_evidence_vec[*idx as usize]));
        }
        v
    }

    pub fn get_sample_idx_by_name(&self, name: &str) -> Option<String> {
        if let Some(idx) = self.name2idx.get(name) {
            Some(idx.to_string())
        } else {
            None
        }
    }
}
