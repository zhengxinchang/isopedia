use crate::io::{DBInfos, GeneralOutputIO, Header, Line};
use crate::meta::Meta;
use crate::writer::MyGzWriter;
use anyhow::Result;
use std::path::Path;

// pub trait GeneralOutputIO: Sized {
//     fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String;

//     fn from_reader<R: BufRead>(reader: &mut R, prefix: Option<&str>, sep: Option<&str>) -> Result<Self>;

//     fn from_file<P: AsRef<Path>>(path: P, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
//         let file = File::open(path)?;
//         let mut buf = BufReader::new(file);
//         Self::from_reader(&mut buf, prefix, sep)
//     }

//     fn from_str(record: &str, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
//         Self::from_reader(&mut record.as_bytes(), prefix, sep)
//     }
// }

// #[derive(Debug, Clone)]
// pub struct Line {
//     pub field_vec: Vec<String>,
//     pub format_str: Option<String>,
//     pub sample_vec: Vec<SampleChip>,
// }

// impl Line {

//     pub fn new() -> Self {
//         Line {
//             field_vec: Vec::new(),
//             format_str: None,
//             sample_vec: Vec::new(),
//         }
//     }

//     pub fn add_sample(&mut self, sample: SampleChip) {
//         self.sample_vec.push(sample);
//     }

//     pub fn update_format_str(&mut self, format_str: &str) {
//         self.format_str = Some(format_str.to_string());
//     }

//     pub fn add_field(&mut self, field: &str) {
//         self.field_vec.push(field.to_string());
//     }
// }

// impl GeneralOutputIO for Line {
//     fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String {
//         let mut parts = Vec::new();

//         for value in &self.field_vec {
//             parts.push(value.clone());
//         }

//         // if let Some(sep) = sep {
//         //     parts.push(sep.to_string());
//         // }

//         parts.push(self.format_str.clone().unwrap());

//         for  sample in &self.sample_vec {
//             let sample_str = sample.to_table(None, Some(":"));
//             parts.push(sample_str);
//         }

//         let line = parts.join("\t") + "\n";
//         let line = add_prefix(&line, &prefix);
//         line
//     }

//     fn from_reader<R: BufRead>(reader: &mut R, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {

//         // debug_assert!(sep.is_some());

//         if prefix.is_some() {
//             return Err(anyhow!("Line::from_reader: prefix must be None"));
//         }

//         if sep.is_none()  {
//             return Err(anyhow!("Line::from_reader: sep must be a single character"));
//         }

//         let sep = sep.unwrap();

//         let mut bufline = String::new();
//         reader.read_line(&mut bufline)?;
//         let parts = bufline.trim_end_matches('\n').split(sep).collect::<Vec<&str>>();
//         if parts.len() != 2 {
//             return Err(anyhow::anyhow!(
//                 "Line does not contain the sample part separated by {}",
//                 sep
//             ));
//         }

//         let mut field_vec = Vec::new();
//         let format_str = sep.to_string();
//         let mut sample_vec = Vec::new();
//         // process field part
//         let field_part = parts[0].trim_end_matches('\t').trim();

//         let field_parts = field_part.split('\t').collect::<Vec<&str>>();
//         for field in &field_parts {
//             if field.is_empty() {
//                 return Err(anyhow::anyhow!(
//                     "Field part contains empty field in line: {}",
//                     bufline
//                 ));
//             }
//             field_vec.push(field.to_string());
//         }

//         // process sample part
//         let sample_part = parts[1].trim_start_matches('\t').trim();

//         for sample_str in sample_part.split('\t') {
//             if sample_str.is_empty() {
//                 return Err(anyhow::anyhow!(
//                     "Sample part contains empty sample in line: {}",
//                     bufline
//                 ));
//             }
//             let sample = SampleChip::from_str(sample_str, None, Some(":"))?;
//             sample_vec.push(sample);
//         }

//         Ok(Line {
//             field_vec,
//             format_str: Some(format_str),
//             sample_vec,
//         })
//     }
// }

// #[derive(Debug, Clone)]
// pub struct SampleChip {
//     pub sample_name: Option<String>,
//     pub fields: Vec<String>,
// }

// impl SampleChip {
//     pub fn new(name: Option<String>, fields: Vec<String>) -> Self {
//         SampleChip {
//             sample_name: name,
//             fields,
//         }
//     }
// }

// impl GeneralOutputIO for SampleChip {
//     fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String {
//         let mut parts = Vec::new();
//         for value in &self.fields {
//             parts.push(value.clone());
//         }
//         let line = parts.join(":");
//         line
//     }

//     fn from_reader<R: BufRead>(reader: &mut R, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {

//         if sep.is_none() {
//             return Err(anyhow!("SampleChip::from_reader: sep must be a single character"));
//         }

//         if prefix.is_some() {
//             return Err(anyhow!("SampleChip::from_reader: prefix must be None"));
//         }

//         let content =  {
//             let mut buf = String::new();
//             reader.read_to_string(&mut buf)?;
//             buf
//         };

//         let sep = sep.unwrap();

//         let content = content.trim();
//         let parts = content.split(sep).collect::<Vec<&str>>();

//         let mut sample_fields = Vec::new();
//         for field in &parts {

//             if field.is_empty() {
//                 return Err(anyhow::anyhow!(
//                     "SampleChip contains empty field in line: {}",
//                     content
//                 ));
//             }
//             sample_fields.push(field.to_string());
//         }

//         Ok(SampleChip {
//             sample_name: None,
//             fields: sample_fields,
//         })
//     }
// }

// #[derive(Debug, Clone)]
// pub struct DBInfos {
//     pub sample_total_evidence_map: IndexMap<String, u32>,
// }

// impl  DBInfos {
//     pub fn new() -> Self {
//         DBInfos {
//             sample_total_evidence_map: IndexMap::new(),
//         }
//     }

//     pub fn add_sample_evidence(&mut self, sample: &str, total: u32) {
//         self.sample_total_evidence_map.insert(sample.to_string(), total);
//     }
// }

// impl GeneralOutputIO for DBInfos {
//     fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String {
//         let mut parts = Vec::new();
//         for (sample, total) in &self.sample_total_evidence_map {
//             parts.push(format!("{}:{}", sample, total));
//         }
//         let table = parts.join(",") + "\n";
//         let table = add_prefix(&table, &prefix);
//         table
//     }

//     fn from_reader<R: BufRead>(reader: &mut R, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
//         let mut content = String::new();
//         reader.read_to_string(&mut content)?; // 10 MB limit
//         let content = if !prefix.is_none() {
//             content.trim_start_matches(prefix.unwrap()).trim()
//         } else {
//             content.trim()
//         };

//         let mut sample_total_evidence_map = IndexMap::new();
//         for part in content.split(',') {
//             let kv: Vec<&str> = part.split(':').collect();
//             if kv.len() != 2 {
//                 return Err(anyhow::anyhow!(
//                     "DbInfos line is not in the correct format: {}",
//                     content
//                 ));
//             }
//             sample_total_evidence_map.insert(kv[0].to_string(), kv[1].parse::<u32>()?);
//         }

//         Ok(DBInfos {
//             sample_total_evidence_map,
//         })
//     }
// }

// #[derive(Debug, Clone)]
// pub struct Header {
//     pub columns: Vec<String>,
//     pub sample_names: Vec<String>,
// }

// impl Header {
//     pub fn new() -> Self {

//         Header {
//             columns: Vec::new(),
//             sample_names: Vec::new(),
//         }
//     }

//     pub fn add_column(&mut self, column: &str) -> Result<()> {
//         self.columns.push(column.to_string());
//         Ok(())
//     }

//     pub fn add_sample_name(&mut self, sample_name: &str) -> Result<()> {
//         self.sample_names.push(sample_name.to_string());
//         Ok(())
//     }
// }

// impl GeneralOutputIO for Header {
//     fn to_table(&self, prefix: Option<&str>, sep: Option<&str>) -> String {
//         let mut header_line = self.columns.join("\t");
//         let sample_line = self.sample_names.join("\t");

//         dbg!(&header_line);
//         dbg!(&sample_line);
//         header_line.push_str("\t");
//         header_line.push_str(FORMAT_STR);
//         header_line.push_str("\t");
//         header_line.push_str(&sample_line);
//         header_line.push_str("\n");
//         let header_line = add_prefix(&header_line, &prefix);
//         header_line
//     }

//     fn from_reader<R: BufRead>(reader: &mut R, prefix: Option<&str>, sep: Option<&str>) -> Result<Self> {
//         let mut header_line = String::new();
//         reader.read_line(&mut header_line)?;

//         let header_line = if !prefix.is_none() {
//             header_line.trim_start_matches(prefix.unwrap()).trim()
//         } else {
//             header_line.trim()
//         };

//         let parts = header_line.split(FORMAT_STR).collect::<Vec<&str>>();
//         if parts.len() != 2 {
//             return Err(anyhow::anyhow!("Header line does not contain FORMAT field"));
//         }

//         let columns_part = parts[0].trim_end_matches('\t').trim();
//         let columns: Vec<String> = columns_part.split('\t').map(|s| s.to_string()).collect();
//         let samples_part = parts[1].trim_start_matches('\t').trim();
//         let samples: Vec<String> = samples_part.split('\t').map(|s| s.to_string()).collect();

//         Ok(Header {
//             columns,
//             sample_names: samples,
//         })
//     }
// }

/*

GeneralOutput trait defines the interface for output handling

*/

pub trait GeneralTableOutput {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self
    where
        Self: Sized;
    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized;
    fn save_to_file<P: AsRef<Path>>(&self, out: P) -> Result<()>;
    fn add_line(&mut self, line: &Line) -> Result<()>;

    fn merge(&mut self, other: &Self) -> Result<()>;
}

pub struct IsoformOutput {
    header: Header,
    db_infos: DBInfos,
    meta: Meta,
    lines: Vec<Line>,
    format_str: String,
}

impl GeneralTableOutput for IsoformOutput {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self {
        IsoformOutput {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn save_to_file<P: AsRef<Path>>(&self, out: P) -> Result<()> {
        // let mut file = File::create(out)?;
        // let mut writer = BufWriter::new(&mut file);

        let mut mywriter = MyGzWriter::new(out.as_ref())?;

        // write meta
        let meta_str = self.meta.to_table(Some("##[SAMPLE]"), None);
        mywriter.write_all_bytes(meta_str.as_bytes())?;
        // write db infos
        let dbinfo_str = self.db_infos.to_table(Some("##[DBINFO]"), None);
        mywriter.write_all_bytes(dbinfo_str.as_bytes())?;
        // write header
        let header_str = self.header.to_table(Some("#"), None);
        mywriter.write_all_bytes(header_str.as_bytes())?;

        // write lines
        for line in &self.lines {
            let line_str = line.to_table(None, None);
            mywriter.write_all_bytes(line_str.as_bytes())?;
        }
        Ok(())
    }

    fn add_line(&mut self, line: &Line) -> Result<()> {
        self.lines.push(line.clone());
        Ok(())
    }

    fn merge(&mut self, other: &Self) -> Result<()> {
        todo!()
    }
}
