use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
    vec,
};
#[derive(Serialize, Deserialize, Debug)]

pub struct DatasetInfo {
    name2idx: FxHashMap<String, u32>,
    name_vec: Vec<String>,
    path_vec: Vec<PathBuf>,
    isoform_size: Vec<u32>,
    total_evidence_vec: Vec<u32>,
    pub sample_size: usize,
}

impl DatasetInfo {
    pub fn new() -> DatasetInfo {
        DatasetInfo {
            name2idx: FxHashMap::default(),
            name_vec: Vec::new(),
            path_vec: Vec::new(),
            isoform_size: Vec::new(),
            total_evidence_vec: Vec::new(),
            sample_size: 0,
        }
    }

    /// parse sample from tab-separated file
    pub fn parse_manifest(path: &PathBuf) -> DatasetInfo {
        let content = std::fs::read_to_string(path).unwrap();

        let mut meta = DatasetInfo::new();
        for line in content.lines().skip(1) {
            let fields = line
                .split('\t')
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let name = fields[1].clone();
            let path = PathBuf::from(&fields[0]);
            meta.add_sample(name, Some(path));
        }
        meta.isoform_size = vec![0; meta.sample_size];
        meta
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
        self.isoform_size[idx] += evidence;
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
}
