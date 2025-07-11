use std::{
    fs::File,
    io::{self, BufRead, Write},
    path::Path,
};

use anyhow::Result;
use indexmap::IndexMap;
use log::info;

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
}

impl Meta {
    pub fn parse<P: AsRef<Path>>(path: P) -> Result<Meta> {
        let mut reader = io::BufReader::new(File::open(path)?);

        let mut header = String::new();
        let mut records = IndexMap::new();
        reader.read_line(&mut header)?;

        let header = header
            .trim_end()
            .split('\t')
            .map(|s| s.to_string())
            .collect::<Vec<String>>();

        info!("Parsed {} fields from header: {:?}", header.len(), header);
        let mut samples = Vec::new();
        let mut line_no = 0;
        for line in reader.lines() {
            let line = line?;
            line_no += 1;
            let parts: Vec<String> = line
                .trim_end()
                .split('\t')
                // .into_iter()
                .map(|part| part.to_string())
                .collect();

            if parts.len() != header.len() {
                return Err(anyhow::anyhow!(
                    "The {} record does not match header length: {} != {}",
                    line_no,
                    parts.len(),
                    header.len()
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
        })
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut writer = File::create(path)?;
        writeln!(writer, "{}", self.header.join("\t"))?;

        for sample in &self.samples {
            if let Some(entry) = self.records.get(sample) {
                let fields: Vec<String> = self
                    .header
                    .iter()
                    .map(|h| entry.fields.get(h).cloned().unwrap_or_default())
                    .collect();
                writeln!(writer, "{}\t{}", entry.name, fields.join("\t"))?;
            }
        }

        Ok(())
    }

    pub fn get_record_by_name(&self, name: &str) -> Option<&MetaEntry> {
        self.records.get(name)
    }

    pub fn get_samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    pub fn get_attr_by_name(&self, name: &str, attr: &str) -> Option<&String> {
        self.records
            .get(name)
            .and_then(|entry| entry.fields.get(attr))
    }
}
// wrote tests for Meta
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::NamedTempFile;
    #[test]
    fn test_meta_parse() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "Sample\tField1\tField2").unwrap();
        writeln!(temp_file, "Sample1\tValue1\tValue2").unwrap();
        writeln!(temp_file, "Sample2\tValue3\tValue4").unwrap();
        temp_file.flush().unwrap(); // Ensure all data is written
        let meta = Meta::parse(temp_file.path()).unwrap();
        assert_eq!(meta.samples.len(), 2);
        assert_eq!(meta.samples[0], "Sample1");
        assert_eq!(meta.samples[1], "Sample2");
        assert_eq!(meta.header.len(), 3);
        assert_eq!(meta.header[0], "Sample");
        assert_eq!(meta.header[1], "Field1");
        assert_eq!(meta.header[2], "Field2");
        assert_eq!(meta.records.len(), 2);
        assert!(meta.records.contains_key("Sample1"));
        assert!(meta.records.contains_key("Sample2"));
        let sample1 = meta.records.get("Sample1").unwrap();
        assert_eq!(sample1.name, "Sample1");
        assert_eq!(sample1.fields.get("Field1").unwrap(), "Value1");
        assert_eq!(sample1.fields.get("Field2").unwrap(), "Value2");
        let sample2 = meta.records.get("Sample2").unwrap();
        assert_eq!(sample2.name, "Sample2");
        assert_eq!(sample2.fields.get("Field1").unwrap(), "Value3");
        assert_eq!(sample2.fields.get("Field2").unwrap(), "Value4");
        dbg!(&meta);
    }
}
