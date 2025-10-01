use std::{
    fs::File,
    io::{self, BufRead, Write},
    path::Path,
};

use anyhow::Result;
use indexmap::IndexMap;
use log::{info, warn};

use crate::{output::TableOutput, utils};

#[derive(Debug, Clone)]
pub struct MetaEntry {
    pub name: String,
    pub fields: IndexMap<String, String>,
}

impl MetaEntry {
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
    pub samples: Vec<String>,
    pub header: Vec<String>,
    pub records: IndexMap<String, MetaEntry>,
    is_empty: bool,
}

impl Meta {
    pub fn new_empty(sample_names: Vec<String>) -> Self {
        let mut records = IndexMap::new();
        let header = vec!["Sample".to_string(), "Path".to_string()];
        for sample in &sample_names {
            let mut entry = MetaEntry::new(sample.clone());
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
        Meta::from_table(&path, prefix)
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut writer = File::create(path)?;
        write!(writer, "{}", self.to_table(None))?; // for internal use, no prefix
        Ok(())
    }

    pub fn get_record_by_name(&self, name: &str) -> Option<&MetaEntry> {
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

    pub fn parse_meta_table<P: AsRef<Path>>(path: P) -> Result<Meta> {
        Err(anyhow::anyhow!("Not implemented yet"))
    }

    pub fn get_attr_by_name(&self, name: &str, attr: &str) -> Option<&String> {
        self.records
            .get(name)
            .and_then(|entry| entry.fields.get(attr))
    }
}

impl TableOutput for Meta {
    fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String {
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

    fn from_table<P: AsRef<Path>>(
        table: &P,
        prefix: Option<&str>,
        sep: Option<&str>,
    ) -> Result<Self>
    where
        Self: Sized,
    {
        let mut reader = io::BufReader::new(File::open(table)?);

        let mut header = String::new();
        let mut records = IndexMap::new();
        reader.read_line(&mut header)?;

        let header = if !prefix.is_none() {
            header.trim_start_matches(prefix.unwrap()).trim()
        } else {
            header.trim()
        };

        let header = utils::line2fields(&header);

        info!("Parsed {} fields from header: {:?}", header.len(), header);
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
                    "The record at line {} does not match header length: {} != {}\nHeader: {:?}\nRecord: {:?}",
                    line_no,
                    parts.len(),
                    header.len(),
                    header,
                    parts
                ));
            }

            samples.push(parts[0].clone());

            let mut meta_entry = MetaEntry::new(parts[0].clone());
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
