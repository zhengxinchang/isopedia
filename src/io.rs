use crate::{constants::FORMAT_STR_NAME, utils::add_prefix};
use anyhow::anyhow;
use anyhow::Result;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use indexmap::IndexMap;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct MyGzWriter {
    #[allow(dead_code)]
    file_name: String,
    inner: Option<GzEncoder<BufWriter<File>>>,
}

impl MyGzWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(&path)?;

        let fname = path.as_ref().to_string_lossy();
        let writer = GzEncoder::new(BufWriter::new(file), Compression::default());

        Ok(MyGzWriter {
            file_name: fname.to_string(),
            inner: Some(writer),
        })
    }

    pub fn write_all_bytes(&mut self, bytes: &[u8]) -> io::Result<()> {
        match &mut self.inner {
            Some(writer) => writer.write_all(bytes),
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been finished",
            )),
        }
    }

    pub fn path(&self) -> &str {
        &self.file_name
    }

    pub fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            Some(writer) => writer.flush(),
            None => Ok(()), // 已完成，无需 flush
        }
    }

    pub fn finish(&mut self) -> io::Result<()> {
        if let Some(writer) = self.inner.take() {
            // take() 移出所有权
            // let mut buf_writer =
            writer.finish()?;
            // buf_writer.flush()?;
        }
        Ok(())
    }
}

impl Drop for MyGzWriter {
    fn drop(&mut self) {
        self.finish().expect("Can not drop gz writer...");
    }
}

pub struct MyGzReader {
    #[allow(dead_code)]
    file_name: String,
    inner: Option<BufReader<GzDecoder<File>>>,
}

impl MyGzReader {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(&path)?;

        let fname = path.as_ref().to_string_lossy();
        let reader = BufReader::new(GzDecoder::new(file));

        Ok(MyGzReader {
            file_name: fname.to_string(),
            inner: Some(reader),
        })
    }

    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        match &mut self.inner {
            Some(reader) => reader.read_line(buf),
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Reader has been finished",
            )),
        }
    }
}

pub trait GeneralOutputIO: Sized {
    fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String;

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        sep: Option<&str>,
    ) -> Result<Self>;

    fn from_file<P: AsRef<Path>>(path: P, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
        let file = File::open(path)?;
        let mut buf = BufReader::new(file);
        Self::from_reader(&mut buf, prefix, sep)
    }

    fn from_str(record: &str, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
        Self::from_reader(&mut record.as_bytes(), prefix, sep)
    }
}

#[derive(Debug, Clone)]
pub struct Line {
    pub field_vec: Vec<String>,
    pub format_str: Option<String>,
    pub sample_vec: Vec<SampleChip>,
}

impl Line {
    pub fn new() -> Self {
        Line {
            field_vec: Vec::new(),
            format_str: None,
            sample_vec: Vec::new(),
        }
    }

    pub fn add_sample(&mut self, sample: SampleChip) {
        self.sample_vec.push(sample);
    }

    pub fn update_format_str(&mut self, format_str: &str) {
        self.format_str = Some(format_str.to_string());
    }

    pub fn add_field(&mut self, field: &str) {
        self.field_vec.push(field.to_string());
    }

    // pub fn merge(&mut self, other: &Line) -> Result<()> {
    //     if self.field_vec != other.field_vec {

    //         error!("Self field_vec: {:?}", self.field_vec);
    //         error!("Other field_vec: {:?}", other.field_vec);

    //         return Err(anyhow::anyhow!("Cannot merge lines with different field_vec"));
    //     }
    //     if self.format_str != other.format_str {
    //         return Err(anyhow::anyhow!("Cannot merge lines with different format_str"));
    //     }
    //     self.sample_vec.extend(other.sample_vec.clone());
    //     Ok(())
    // }
}

impl GeneralOutputIO for Line {
    /// sep does not work for Line
    fn to_table(&self, prefix: Option<&str>, _sep: Option<&str>) -> String {
        let mut parts = Vec::new();

        for value in &self.field_vec {
            parts.push(value.clone());
        }

        parts.push(self.format_str.clone().unwrap());

        for sample in &self.sample_vec {
            let sample_str = sample.to_table(None, Some(":"));
            parts.push(sample_str);
        }

        let line = parts.join("\t") + "\n";
        let line = add_prefix(&line, &prefix);
        line
    }

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        sep: Option<&str>,
    ) -> Result<Self> {
        // debug_assert!(sep.is_some());

        if prefix.is_some() {
            return Err(anyhow!("Line::from_reader: prefix must be None"));
        }

        if sep.is_none() {
            return Err(anyhow!("Line::from_reader: sep must be a single character"));
        }

        let sep = sep.unwrap();

        let mut bufline = String::new();
        reader.read_line(&mut bufline)?;
        let parts = bufline
            .trim_end_matches('\n')
            .split(sep)
            .collect::<Vec<&str>>();
        if parts.len() != 2 {
            return Err(anyhow::anyhow!(
                "Line does not contain the sample part separated by {}",
                sep
            ));
        }

        let mut field_vec = Vec::new();
        let format_str = sep.to_string();
        let mut sample_vec = Vec::new();
        // process field part
        let field_part = parts[0].trim_end_matches('\t').trim();

        let field_parts = field_part.split('\t').collect::<Vec<&str>>();
        for field in &field_parts {
            if field.is_empty() {
                return Err(anyhow::anyhow!(
                    "Field part contains empty field in line: {}",
                    bufline
                ));
            }
            field_vec.push(field.to_string());
        }

        // process sample part
        let sample_part = parts[1].trim_start_matches('\t').trim();

        for sample_str in sample_part.split('\t') {
            if sample_str.is_empty() {
                return Err(anyhow::anyhow!(
                    "Sample part contains empty sample in line: {}",
                    bufline
                ));
            }
            let sample = SampleChip::from_str(sample_str, None, Some(":"))?;
            sample_vec.push(sample);
        }

        Ok(Line {
            field_vec,
            format_str: Some(format_str),
            sample_vec,
        })
    }
}

#[derive(Debug, Clone)]
pub struct SampleChip {
    pub sample_name: Option<String>,
    pub fields: Vec<String>,
}

impl SampleChip {
    pub fn default(name: Option<String>) -> Self {
        let f = Vec::new();
        SampleChip {
            sample_name: name,
            fields: f,
        }
    }

    pub fn new(name: Option<String>, fields: Vec<String>) -> Self {
        SampleChip {
            sample_name: name,
            fields,
        }
    }

    pub fn add_item(&mut self, item: &str) {
        self.fields.push(item.to_string());
    }
}

impl GeneralOutputIO for SampleChip {
    fn to_table(&self, _prefix: Option<&str>, _sep: Option<&str>) -> String {
        let mut parts = Vec::new();
        for value in &self.fields {
            parts.push(value.clone());
        }
        let line = parts.join(":");
        line
    }

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        sep: Option<&str>,
    ) -> Result<Self> {
        if sep.is_none() {
            return Err(anyhow!(
                "SampleChip::from_reader: sep must be a single character"
            ));
        }

        if prefix.is_some() {
            return Err(anyhow!("SampleChip::from_reader: prefix must be None"));
        }

        let content = {
            let mut buf = String::new();
            reader.read_to_string(&mut buf)?;
            buf
        };

        let sep = sep.unwrap();

        let content = content.trim();
        let parts = content.split(sep).collect::<Vec<&str>>();

        let mut sample_fields = Vec::new();
        for field in &parts {
            if field.is_empty() {
                return Err(anyhow::anyhow!(
                    "SampleChip contains empty field in line: {}",
                    content
                ));
            }
            sample_fields.push(field.to_string());
        }

        Ok(SampleChip {
            sample_name: None,
            fields: sample_fields,
        })
    }
}

#[derive(Debug, Clone)]
pub struct DBInfos {
    pub sample_total_evidence_map: IndexMap<String, u32>,
}

impl DBInfos {
    pub fn new() -> Self {
        DBInfos {
            sample_total_evidence_map: IndexMap::new(),
        }
    }

    pub fn add_sample_evidence(&mut self, sample: &str, total: u32) {
        self.sample_total_evidence_map
            .insert(sample.to_string(), total);
    }

    pub fn merge(&mut self, other: &DBInfos) -> Result<()> {
        for (sample, total) in &other.sample_total_evidence_map {
            if self.sample_total_evidence_map.contains_key(sample) {
                return Err(anyhow::anyhow!(
                    "Duplicate sample name found when merging DBInfos: {}",
                    sample
                ));
            }
            self.sample_total_evidence_map
                .insert(sample.clone(), *total);
        }
        Ok(())
    }

    pub fn get_total_evidence_vec(&self) -> (Vec<String>, Vec<u32>) {
        let mut vec = Vec::new();
        let mut vec_u32 = Vec::new();
        for (sample, total) in &self.sample_total_evidence_map {
            vec.push(sample.clone());
            vec_u32.push(*total);
        }
        (vec, vec_u32)
    }
}

impl GeneralOutputIO for DBInfos {
    fn to_table(&self, prefix: Option<&str>, _sep: Option<&str>) -> String {
        let mut parts = Vec::new();
        for (sample, total) in &self.sample_total_evidence_map {
            parts.push(format!("{}:{}", sample, total));
        }
        let table = parts.join(",") + "\n";
        let table = add_prefix(&table, &prefix);
        table
    }

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        _sep: Option<&str>,
    ) -> Result<Self> {
        let mut content = String::new();
        reader.read_to_string(&mut content)?; // 10 MB limit
        let content = if !prefix.is_none() {
            content.trim_start_matches(prefix.unwrap()).trim()
        } else {
            content.trim()
        };

        let mut sample_total_evidence_map = IndexMap::new();
        for part in content.split(',') {
            let kv: Vec<&str> = part.split(':').collect();
            if kv.len() != 2 {
                return Err(anyhow::anyhow!(
                    "DbInfos line is not in the correct format: {}",
                    content
                ));
            }
            sample_total_evidence_map.insert(kv[0].to_string(), kv[1].parse::<u32>()?);
        }

        Ok(DBInfos {
            sample_total_evidence_map,
        })
    }
}

#[derive(Debug, Clone)]
pub struct Header {
    pub columns: Vec<String>,
    pub sample_names: Vec<String>,
}

impl Header {
    pub fn new() -> Self {
        Header {
            columns: Vec::new(),
            sample_names: Vec::new(),
        }
    }

    pub fn add_column(&mut self, column: &str) -> Result<()> {
        self.columns.push(column.to_string());
        Ok(())
    }

    pub fn add_sample_name(&mut self, sample_name: &str) -> Result<()> {
        self.sample_names.push(sample_name.to_string());
        Ok(())
    }

    pub fn merge(&mut self, other: &Header) -> Result<()> {
        for sample in &other.sample_names {
            if !self.sample_names.contains(sample) {
                self.sample_names.push(sample.clone());
            } else {
                return Err(anyhow::anyhow!(
                    "Duplicate sample name found when merging headers: {}",
                    sample
                ));
            }
        }
        Ok(())
    }
}

impl GeneralOutputIO for Header {
    fn to_table(&self, prefix: Option<&str>, _sep: Option<&str>) -> String {
        let mut header_line = self.columns.join("\t");
        let sample_line = self.sample_names.join("\t");
        header_line.push_str("\t");
        header_line.push_str(FORMAT_STR_NAME);
        header_line.push_str("\t");
        header_line.push_str(&sample_line);
        header_line.push_str("\n");
        let header_line = add_prefix(&header_line, &prefix);
        header_line
    }

    fn from_reader<R: BufRead>(
        reader: &mut R,
        prefix: Option<&str>,
        _sep: Option<&str>,
    ) -> Result<Self> {
        let mut header_line = String::new();
        reader.read_line(&mut header_line)?;

        let header_line = if !prefix.is_none() {
            header_line.trim_start_matches(prefix.unwrap()).trim()
        } else {
            header_line.trim()
        };
        // dbg!(&header_line);
        let parts = header_line.split(FORMAT_STR_NAME).collect::<Vec<&str>>();
        if parts.len() != 2 {
            return Err(anyhow::anyhow!("Header line does not contain FORMAT field"));
        }

        let columns_part = parts[0].trim_end_matches('\t').trim();
        let columns: Vec<String> = columns_part.split('\t').map(|s| s.to_string()).collect();
        let samples_part = parts[1].trim_start_matches('\t').trim();
        let samples: Vec<String> = samples_part.split('\t').map(|s| s.to_string()).collect();

        Ok(Header {
            columns,
            sample_names: samples,
        })
    }
}

// pub struct OutputVersion{
//     type_: String,
//     version: String,
// }

// impl OutputVersion {
//     pub fn new(type_: &str, version: &str) -> Self {
//         OutputVersion {
//             type_: type_.to_string(),
//             version: version.to_string(),
//         }
//     }
// }

// impl GeneralOutputIO for OutputVersion {
//     fn to_table(&self, prefix: Option<&str>, _sep: Option<&str>) -> String {
//         let mut parts = Vec::new();
//         parts.push(format!("type={}",self.type_));
//         parts.push(format!("version={}",self.version));
//         let line = parts.join(";") + "\n";
//         let line = add_prefix(&line, &prefix);
//         line
//     }

//     fn from_reader<R: BufRead>(
//         reader: &mut R,
//         prefix: Option<&str>,
//         _sep: Option<&str>,
//     ) -> Result<Self> {
//         let mut bufline = String::new();
//         reader.read_line(&mut bufline)?;
//         let line = if !prefix.is_none() {
//             bufline.trim_start_matches(prefix.unwrap()).trim()
//         } else {
//             bufline.trim()
//         };
//         let mut type_ = String::new();
//         let mut version = String::new();
//         for part in line.split(';') {
//             let mut kv = part.splitn(2, '=');
//             let key = kv.next().unwrap();
//             let value = kv.next().unwrap();
//             match key {
//                 "type" => type_ = value.to_string(),
//                 "version" => version = value.to_string(),
//                 _ => {}
//             }
//         }
//         Ok(OutputVersion { type_, version })
//     }
// }
