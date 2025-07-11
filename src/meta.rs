use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{fs::{self, File}, io::Write, path::{Path, PathBuf}, vec};
#[derive(Serialize, Deserialize, Debug)]

pub struct Meta {
    name2idx: FxHashMap<String, u32>,
    name_vec: Vec<String>,
    path_vec: Vec<PathBuf>,
    isoform_size: Vec<u32>,
    sample_size: usize,
}

impl Meta {
    pub fn new() -> Meta {
        Meta {
            name2idx: FxHashMap::default(),
            name_vec: Vec::new(),
            path_vec: Vec::new(),
            isoform_size: Vec::new(),
            sample_size: 0,
        }
    }

    /// parse sample from tab-separated file
    pub fn parse_manifest(path: &PathBuf) -> Meta {
        let content = std::fs::read_to_string(path).unwrap();
        // let header = content.lines().next().unwrap();
        // let meta_keys = header
        //     .split('\t')
        //     .map(|x| x.to_string())
        //     .collect::<Vec<String>>();
        // let meta_keys = meta_keys[2..].to_vec();
        // dbg!(&meta_keys);
        let mut meta = Meta::new();
        for line in content.lines().skip(1) {
            let fields = line
                .split('\t')
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let name = fields[1].clone();
            let path = PathBuf::from(&fields[0]);
            // let meta_vals = fields[2..].to_vec();
            meta.add_sample(name, Some(path));
        }
        meta.isoform_size = vec![0; meta.sample_size];
        meta
    }

    /// load sample meta from index
    pub fn load(path: &PathBuf) -> Meta {
        let content = std::fs::read_to_string(path)
            .expect(&format!("Cannot read sample meta file {}", path.display()));
        let mut samples = Meta::new();
        // let header = content.lines().next().unwrap();
        // let meta_keys = header
        //     .split('\t')
        //     .skip(1)
        //     .map(|x| x.to_string())
        //     .collect::<Vec<String>>();

        for line in content.lines().skip(1) {
            let fields = line
                .split('\t')
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let name = fields[0].clone();
            // let meta_vals = fields[1..].to_vec();
            samples.add_sample(name, None);
        }
        samples
    }

    /// dump sample meta to index
    // pub fn get_string(&self) -> String {
    //     let mut content = "sample_name".to_string();
    //     for (k, _) in self.meta_vec[0].iter() {
    //         content.push_str(&format!("\t{}", k));
    //     }
    //     content.push_str("\n");
    //     for (name, meta) in self.name_vec.iter().zip(self.meta_vec.iter()) {
    //         content.push_str(&format!("{}", name));
    //         for (_, v) in meta.iter() {
    //             content.push_str(&format!("\t{}", v));
    //         }
    //         content.push_str("\n");
    //     }
    //     content
    // }

    pub fn save_to_file<P: AsRef<Path>>(self, path: P) -> std::io::Result<()> {
        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        let mut file = File::create(path)?;
        file.write_all(&encoded)?;
        Ok(())
    }

    pub fn load_from_file(path: &str) -> std::io::Result<Meta> {
        let bytes = fs::read(path)?;
        let decoded: Meta = bincode::deserialize(&bytes).unwrap();
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

    // pub fn get_sample_name(&self, idx: u32) -> &str {
    //     &self.name_vec[idx as usize]
    // }

    // pub fn get_sample_idx(&self, name: &str) -> Option<u32> {
    //     self.name2idx.get(name).map(|idx| *idx)
    // }

    pub fn get_size(&self) -> usize {
        self.sample_size
    }

    pub fn get_sample_names(&self) -> Vec<String> {
        self.name_vec.clone()
    }
}
