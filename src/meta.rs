use anyhow::Result;
use indexmap::IndexMap;
use log::{error, info, warn};
use std::{
    fs::File,
    io::{BufRead, Write},
    path::Path,
};

use crate::{io::GeneralOutputIO, utils};

#[derive(Debug, Clone)]
pub struct SampleMetaEntry {
    // sample entry
    pub name: String,
    pub fields: IndexMap<String, String>, // attr -> value
}

impl SampleMetaEntry {
    fn new(name: String) -> Self {
        Self {
            name,
            fields: IndexMap::new(),
        }
    }

    fn add(&mut self, key: String, value: String) -> Result<()> {
        self.fields.insert(key, value);
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct Meta {
    pub samples: Vec<String>,                       // sample names
    pub header: Vec<String>,                        // atrribute names
    pub records: IndexMap<String, SampleMetaEntry>, // sample name -> MetaEntry
    is_empty: bool,
}

impl Meta {
    pub fn new_empty(sample_names: Vec<String>) -> Self {
        let mut records = IndexMap::new();
        let header = vec!["Sample".to_string(), "Path".to_string()];
        for sample in &sample_names {
            let mut entry = SampleMetaEntry::new(sample.clone());
            entry
                .add("Path".to_string(), "Fake_path".to_string())
                .expect("Failed to add Path field");
            records.insert(sample.clone(), entry);
        }
        Self {
            samples: sample_names,
            header,
            records,
            is_empty: true,
        }
    }
    pub fn parse<P: AsRef<Path>>(path: P, prefix: Option<&str>) -> Result<Meta> {
        info!("If tab is detected in the line, it will be used as the field separator,otherwise, space will be used as the field separator.");
        let meta = Meta::from_file(&path, prefix, None)?;
        info!(
            "Parsed meta file: {} with {} samples and {} attributes: {:?}",
            path.as_ref().display(),
            meta.samples.len(),
            meta.header.len(),
            meta.header
        );
        Ok(meta)
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut writer = File::create(path)?;
        write!(writer, "{}", self.to_table(None, None))?; // for internal use, no prefix
        Ok(())
    }

    pub fn get_record_by_name(&self, name: &str) -> Option<&SampleMetaEntry> {
        if self.is_empty {
            return None;
        }
        self.records.get(name)
    }

    pub fn get_samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    pub fn get_meta_table(&self, prefix: Option<&str>) -> String {
        let mut table = String::new();
        let mut header_clean = self.header.clone();
        header_clean.remove(1); // remove the path

        if let Some(p) = prefix {
            table.push_str(&format!("{}{}\n", p, header_clean.join("\t")));
            for sample in &self.samples {
                if let Some(entry) = self.records.get(sample) {
                    let fields: Vec<String> = header_clean
                        .iter()
                        .map(|h| entry.fields.get(h).cloned().unwrap_or_default())
                        .collect();
                    table.push_str(&format!("{}{}{}\n", p, entry.name, fields.join("\t")));
                }
            }
        } else {
            table.push_str(&format!("{}\n", self.header.join("\t")));
            for sample in &self.samples {
                if let Some(entry) = self.records.get(sample) {
                    let fields: Vec<String> = self
                        .header
                        .iter()
                        .map(|h| entry.fields.get(h).cloned().unwrap_or_default())
                        .collect();
                    table.push_str(&format!("{}{}\n", entry.name, fields.join("\t")));
                }
            }
        }

        table
    }

    /// remove the Path column from header and records
    pub fn remove_path_column(&mut self) {
        // remove the second column in header and records
        if self.header.len() < 2 {
            return;
        }
        let path_col = self.header[1].clone();
        self.header.remove(1);
        for (_sample, metaentry) in self.records.iter_mut() {
            metaentry.fields.shift_remove(&path_col);
        }
    }

    pub fn validate_path_column(&self) {
        //if the value in the scond column is a valid path

        let path_col = self.header[1].clone();
        for (sample, metaentry) in self.records.iter() {
            if let Some(path) = metaentry.fields.get(&path_col) {
                if !Path::new(path).exists() {
                    error!("The path for sample {} does not exist: {}", sample, path);
                    std::process::exit(1);
                }
            } else {
                error!(
                    "The path column {} is missing for sample {}",
                    path_col, sample
                );
                std::process::exit(1);
            }
        }
    }

    pub fn validate_sample_names(&self) {
        let mut name_set = std::collections::HashSet::new();
        for sample in &self.samples {
            if name_set.contains(sample) {
                error!("Duplicate sample name found: {}", sample);
                std::process::exit(1);
            } else {
                name_set.insert(sample);
            }
        }
    }

    pub fn merge(&mut self, other: &Meta) -> Result<()> {
        if self.header != other.header {
            warn!("Merging Meta with different headers: {:?} vs {:?}, missing attributes will be filled with 'NA'", self.header, other.header);
        }

        let mut all_headers = Vec::new();

        for h in &self.header {
            if !all_headers.contains(h) {
                all_headers.push(h.clone());
            }
        }

        for h in &other.header {
            if !all_headers.contains(h) {
                all_headers.push(h.clone());
            }
        }

        let mut all_samples = Vec::new();
        for s in &self.samples {
            if !all_samples.contains(s) {
                all_samples.push(s.clone());
            }
        }
        for s in &other.samples {
            if !all_samples.contains(s) {
                all_samples.push(s.clone());
            }
        }

        self.header = all_headers.clone();
        self.samples = all_samples.clone();

        let mut new_records = IndexMap::new();
        for sample in &all_samples {
            let mut entry = SampleMetaEntry::new(sample.clone());

            for attr in &all_headers {
                if attr == "Sample" {
                    continue;
                }
                if let Some(self_entry) = self.records.get(sample) {
                    if let Some(value) = self_entry.fields.get(attr) {
                        entry.add(attr.clone(), value.clone())?;
                        // continue;
                    } else {
                        entry.add(attr.clone(), "NA1".to_string())?;
                    }
                } else if let Some(other_entry) = other.records.get(sample) {
                    if let Some(value) = other_entry.fields.get(attr) {
                        entry.add(attr.clone(), value.clone())?;
                        // continue;
                    } else {
                        entry.add(attr.clone(), "NA2".to_string())?;
                    }
                } else {
                    // entry.add(attr.clone(), "NA".to_string())?;
                }
            }

            new_records.insert(sample.clone(), entry);
        }

        self.records = new_records;

        Ok(())
    }

    pub fn statement(&self) {}
}

impl GeneralOutputIO for Meta {
    /// the sep shouldnt be used
    fn to_table(&self, prefix: Option<&str>, _sep: Option<&str>) -> String {
        let mut table = String::new();
        let mut header_clean = self.header.clone();
        // standarize the fisrt column name with Sample
        header_clean[0] = "Sample".into();
        // header_clean.remove(1); // remove the path
        if let Some(p) = prefix {
            table.push_str(&format!("{}{}\n", p, header_clean.join("\t")));
            for sample in &self.samples {
                if let Some(entry) = self.records.get(sample) {
                    let fields: Vec<String> = header_clean
                        .iter()
                        .map(|h| entry.fields.get(h).cloned().unwrap_or_default())
                        .collect();
                    table.push_str(&format!("{}{}{}\n", p, entry.name, fields.join("\t")));
                }
            }
        } else {
            table.push_str(&format!("{}\n", self.header.join("\t")));
            for sample in &self.samples {
                if let Some(entry) = self.records.get(sample) {
                    let fields: Vec<String> = self
                        .header
                        .iter()
                        .map(|h| entry.fields.get(h).cloned().unwrap_or_default())
                        .collect();
                    table.push_str(&format!("{}{}\n", entry.name, fields.join("\t")));
                }
            }
        }

        table
    }

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        _sep: Option<&str>,
    ) -> Result<Self> {
        let mut header = String::new();
        let mut records = IndexMap::new();
        reader.read_line(&mut header)?;

        let header = if !prefix.is_none() {
            header.trim_start_matches(prefix.unwrap()).trim()
        } else {
            header.trim()
        };

        let header = utils::line2fields(&header);

        // info!("Parsed {} fields from header: {:?}", header.len(), header);

        let mut samples = Vec::new();
        let mut line_no = 0;
        for line in reader.lines() {
            let line = line?;

            let line = if !prefix.is_none() {
                line.trim_start_matches(prefix.unwrap()).trim().to_string()
            } else {
                line.trim().to_string()
            };

            line_no += 1;

            let parts = utils::line2fields(&line);
            if parts.len() == 0 {
                warn!("Skipping empty line at {}", line_no);
                continue;
            }

            if parts.len() != header.len() {
                return Err(anyhow::anyhow!(
                    "The record at line {} does not match header length: {} != {}\n\n> Header: {:?}\n> Record: {:?}\n\n> Please firstly check if the '\\t' and ' ' are both used as separators in the file. \n> Iosopedia does not support mixed separators in the same file, as some values may contain space within fields.\n> If '\\t' is detected in the line, it will be used as the field separator, otherwise, space will be used as the field separator.
                    ",
                    line_no,
                    parts.len(),
                    header.len(),
                    header,
                    parts
                ));
            }

            samples.push(parts[0].clone());

            let mut meta_entry = SampleMetaEntry::new(parts[0].clone());
            for (k, v) in header.iter().zip(&parts).skip(1) {
                meta_entry.add(k.clone(), v.clone())?;
            }
            records.insert(parts[0].clone(), meta_entry);
        }

        Ok(Meta {
            samples,
            header,
            records,
            is_empty: false,
        })
    }
}
