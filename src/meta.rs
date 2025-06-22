use indexmap::IndexMap;
use rustc_hash::FxHashMap;
use std::path::PathBuf;

#[derive(Debug)]
pub struct Meta {
    name2idx: FxHashMap<String, u32>,
    name_vec: Vec<String>,
    path_vec: Vec<PathBuf>,
    meta_vec: Vec<IndexMap<String, String>>,
    size: usize,
}

impl Meta {
    pub fn new() -> Meta {
        Meta {
            name2idx: FxHashMap::default(),
            name_vec: Vec::new(),
            path_vec: Vec::new(),
            meta_vec: Vec::new(),
            size: 0,
        }
    }

    /// parse sample from tab-separated file
    pub fn parse(path: &PathBuf) -> Meta {
        let content = std::fs::read_to_string(path).unwrap();
        let header = content.lines().next().unwrap();
        let meta_keys = header
            .split('\t')
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        let meta_keys = meta_keys[2..].to_vec();
        // dbg!(&meta_keys);
        let mut samples = Meta::new();
        for line in content.lines().skip(1) {
            let fields = line
                .split('\t')
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let name = fields[1].clone();
            let path = PathBuf::from(&fields[0]);
            let meta_vals = fields[2..].to_vec();
            samples.add_sample(name, Some(path), &meta_keys, &meta_vals);
        }
        samples
    }

    /// load sample meta from index
    pub fn load(path: &PathBuf) -> Meta {
        let content = std::fs::read_to_string(path).unwrap();
        let mut samples = Meta::new();
        let header = content.lines().next().unwrap();
        let meta_keys = header
            .split('\t')
            .skip(1)
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        for line in content.lines().skip(1) {
            let fields = line
                .split('\t')
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            let name = fields[0].clone();
            let meta_vals = fields[1..].to_vec();
            samples.add_sample(name, None, &meta_keys, &meta_vals);
        }
        samples
    }

    /// dump sample meta to index
    pub fn get_string(&self) -> String {
        let mut content = "sample_name".to_string();
        for (k, _) in self.meta_vec[0].iter() {
            content.push_str(&format!("\t{}", k));
        }
        content.push_str("\n");
        for (name, meta) in self.name_vec.iter().zip(self.meta_vec.iter()) {
            content.push_str(&format!("{}", name));
            for (_, v) in meta.iter() {
                content.push_str(&format!("\t{}", v));
            }
            content.push_str("\n");
        }
        content
    }

    pub fn add_sample(
        &mut self,
        name: String,
        file_path: Option<PathBuf>,
        meta_key: &Vec<String>,
        meta_val: &Vec<String>,
    ) -> u32 {
        if let Some(idx) = self.name2idx.get(&name) {
            *idx
        } else {
            let idx = self.name_vec.len() as u32;
            self.name2idx.insert(name.clone(), idx);
            self.name_vec.push(name);

            if let Some(f) = file_path {
                self.path_vec.push(f);
            }

            let mut meta = IndexMap::new();
            for (k, v) in meta_key.iter().zip(meta_val.iter()) {
                meta.insert(k.clone(), v.clone());
            }
            self.meta_vec.push(meta);
            self.size += 1;
            idx
        }
    }

    pub fn get_path_list(&self) -> Vec<PathBuf> {
        self.path_vec.clone()
    }

    pub fn get_sample_name(&self, idx: u32) -> &str {
        &self.name_vec[idx as usize]
    }

    pub fn get_sample_idx(&self, name: &str) -> Option<u32> {
        self.name2idx.get(name).map(|idx| *idx)
    }

    pub fn get_size(&self) -> usize {
        self.size
    }

    pub fn get_sample_names(&self) -> Vec<String> {
        self.name_vec.clone()
    }
}
