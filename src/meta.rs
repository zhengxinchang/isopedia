use std::{
    fs::File,
    io::{self, BufRead},
    path::Path,
};

use anyhow::Result;
use indexmap::IndexMap;
use rust_htslib::bam::header::HeaderRecord;

#[derive(Debug, Clone)]
struct MetaEntry {
    name: String,
    fields: IndexMap<String, String>,
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
struct Meta {
    samples: Vec<String>,
    header: Vec<String>,
    records: IndexMap<String, MetaEntry>,
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
        let mut samples = Vec::new();
        for line in reader.lines() {
            let line = line?;

            let parts: Vec<String> = line
                .trim_end()
                .split('\t')
                // .into_iter()
                .map(|part| part.to_string())
                .collect();
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