use crate::constants::FORMAT_STR_NAME;
use crate::io::MyGzWriter;
use crate::io::{DBInfos, GeneralOutputIO, Header, Line, MyGzReader};
use crate::meta::Meta;
use crate::utils;
use anyhow::Result;
use log::info;
use std::path::Path;

/*
GeneralOutput trait defines the interface for output handling
*/

pub trait GeneralTableOutput {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self
    where
        Self: Sized;
    fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized;
    fn load_elements<P: AsRef<Path>>(path: P) -> Result<(Header, DBInfos, Meta, Vec<Line>, String)>
    where
        Self: Sized,
    {
        let mut mygzreader = MyGzReader::new(path.as_ref())?;

        let mut header_strings = String::new();
        let mut db_infos_strings = String::new();
        let mut meta_strings = String::new();
        let mut lines = Vec::new();

        let mut line_buf = String::new();
        while let Ok(bytes) = mygzreader.read_line(&mut line_buf) {
            if bytes == 0 {
                break;
            }
            if line_buf.starts_with("##[DBINFO]") {
                db_infos_strings.push_str(&line_buf);
            } else if line_buf.starts_with("##[SAMPLE]") {
                meta_strings.push_str(&line_buf);
            } else if line_buf.starts_with('#') {
                // must behind all other lines that start with ##
                // dbg!(&line_buf);
                header_strings.push_str(&line_buf);
            } else {
                lines.push(line_buf.clone());
            }

            line_buf.clear();
        }

        let header = Header::from_str(&header_strings, Some("#"), Some(FORMAT_STR_NAME))?;
        let db_infos = DBInfos::from_str(&db_infos_strings, Some("##[DBINFO]"), None)?;
        let meta = Meta::from_str(&meta_strings, Some("##[SAMPLE]"), None)?;

        let first_line = lines[0].clone();
        let first_line_fields = first_line.split('\t').collect::<Vec<&str>>();
        let format_str = first_line_fields.get(header.columns.len()).unwrap();

        let mut line_objs = Vec::new();
        for line in lines {
            let line_obj = Line::from_str(&line, None, Some(format_str))?;
            // add to lines
            line_objs.push(line_obj);
        }

        // dbg!(&line_objs[0]);

        Ok((header, db_infos, meta, line_objs, format_str.to_string()))
    }
    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()>;
    fn add_line(&mut self, line: &Line) -> Result<()>;

    fn merge(&mut self, other: &Self) -> Result<()>;
}

pub struct IsoformTableOut {
    header: Header,
    db_infos: DBInfos,
    pub meta: Meta,
    lines: Vec<Line>,
    #[allow(dead_code)]
    format_str: String,
}

impl IsoformTableOut {
    pub fn get_mem_size(&self) -> usize {
        let mut total = 0usize;

        // header/db_infos/meta 不太大，可以略过
        for line in &self.lines {
            total += std::mem::size_of_val(line);

            // field_vec: Vec<String>
            total += line.field_vec.capacity() * std::mem::size_of::<String>();
            for s in &line.field_vec {
                total += s.capacity();
            }

            // sample_vec: Vec<SampleChip>
            total += line.sample_vec.capacity() * std::mem::size_of::<crate::io::SampleChip>();
            for sample in &line.sample_vec {
                total += sample.fields.capacity() * std::mem::size_of::<String>();
                for s in &sample.fields {
                    total += s.capacity();
                }
            }
        }

        total
    }
}

impl GeneralTableOutput for IsoformTableOut {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self {
        IsoformTableOut {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized,
    {
        let (header, db_infos, meta, line_objs, format_str) = Self::load_elements(path)?;

        Ok(IsoformTableOut {
            header,
            db_infos,
            meta,
            lines: line_objs,
            format_str: format_str.to_string(),
        })
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
        // let outname = utils::add_gz_suffix_if_needed(out.as_ref());
        let outname = utils::add_gz_suffix_if_needed(&out);

        info!("Saving output to {:?}", outname);

        let mut mywriter = MyGzWriter::new(&outname)?;

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
        for line in &mut self.lines {
            if line.format_str.is_none() {
                line.update_format_str(&self.format_str);
            }

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
        if self.header.columns != other.header.columns {
            return Err(anyhow::anyhow!(
                "Cannot merge tables with different headers"
            ));
        }

        if self.lines.len() != other.lines.len() {
            return Err(anyhow::anyhow!(
                "Cannot merge tables with different number of lines"
            ));
        }

        self.meta.merge(&other.meta)?;
        self.header.merge(&other.header)?;
        self.db_infos.merge(&other.db_infos)?;

        // dbg!(&self.db_infos);

        let (sample_name_vec, sample_total_evidence_vec) = self.db_infos.get_total_evidence_vec();
        info!(
            "Merging output tables with {} samples: {:?}",
            sample_name_vec.len(),
            sample_name_vec
        );
        // dbg!(&sample_total_evidence_vec);
        for (i, line) in self.lines.iter_mut().enumerate() {
            let other_line = &other.lines[i];

            //custom merge function for line to handle isoform specific fields

            let mut self_fixed_fields = Vec::new();
            let mut other_fixed_fields = Vec::new();
            for j in vec![0, 1, 2, 3, 4, 5, 6, 11] {
                self_fixed_fields.push(line.field_vec[j].clone());
                other_fixed_fields.push(other_line.field_vec[j].clone());
            }

            if self_fixed_fields != other_fixed_fields {
                return Err(anyhow::anyhow!(
                    "Cannot merge lines with different fixed fields at line {}",
                    i + 1
                ));
            }

            if line.field_vec[9] != other_line.field_vec[9] {
                return Err(anyhow::anyhow!(
                    "Cannot merge lines with different 'min_read' field at line {}",
                    i + 1
                ));
            }

            let self_pos_sample_parts = line.field_vec[10]
                .split('/')
                .map(|s| s.parse::<u64>().unwrap())
                .collect::<Vec<u64>>();
            let other_pos_sample_parts = other_line.field_vec[10]
                .split('/')
                .map(|s| s.parse::<u64>().unwrap())
                .collect::<Vec<u64>>();

            let merged_parts = self_pos_sample_parts
                .iter()
                .zip(other_pos_sample_parts.iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<u64>>();
            // update the is_detected field
            if merged_parts[0] > 0 {
                line.field_vec[8] = "Yes".to_string();
            } else {
                line.field_vec[8] = "No".to_string();
            }
            // update the pos/sample field
            line.field_vec[10] = format!("{}/{}", merged_parts[0], merged_parts[1]);

            // merge the sample fields
            line.sample_vec.extend(other_line.sample_vec.clone());

            let mut acc_sample_evidence_arr = Vec::new();
            for sample in &line.sample_vec {
                // dbg!(&sample);
                acc_sample_evidence_arr.push(sample.fields[1].parse::<u32>().unwrap());
                // evidence is at index 1
            }

            // dbg!(&acc_sample_evidence_arr);
            // dbg!(&sample_total_evidence_vec);
            // dbg!(&sample_name_vec);

            let confidence = crate::isoform::MergedIsoform::get_confidence_value(
                &acc_sample_evidence_arr,
                sample_total_evidence_vec.len(),
                &sample_total_evidence_vec,
            );
            line.field_vec[7] = format!("{}", confidence); // confidence is at index 7
        }
        info!("Finished");
        Ok(())
    }
}

pub struct SpliceTableOut {
    header: Header,
    db_infos: DBInfos,
    meta: Meta,
    lines: Vec<Line>,
    #[allow(dead_code)]
    format_str: String,
}

impl GeneralTableOutput for SpliceTableOut {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self {
        SpliceTableOut {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized,
    {
        let (header, db_infos, meta, line_objs, format_str) = Self::load_elements(path)?;

        Ok(SpliceTableOut {
            header,
            db_infos,
            meta,
            lines: line_objs,
            format_str: format_str.to_string(),
        })
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
        let outname = utils::add_gz_suffix_if_needed(&out);

        let mut mywriter = MyGzWriter::new(&outname)?;

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
        for line in &mut self.lines {
            if line.format_str.is_none() {
                line.update_format_str(&self.format_str);
            }
            let line_str = line.to_table(None, None);
            mywriter.write_all_bytes(line_str.as_bytes())?;
        }
        info!("Saved output to {:?}", outname);
        Ok(())
    }

    fn add_line(&mut self, line: &Line) -> Result<()> {
        self.lines.push(line.clone());
        Ok(())
    }
    #[allow(unused_variables)]
    fn merge(&mut self, other: &Self) -> Result<()> {
        todo!()
    }
}

pub struct FusionBrkPtTableOut {
    header: Header,
    db_infos: DBInfos,
    meta: Meta,
    lines: Vec<Line>,
    #[allow(dead_code)]
    format_str: String,
}

impl GeneralTableOutput for FusionBrkPtTableOut {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self {
        FusionBrkPtTableOut {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized,
    {
        let (header, db_infos, meta, line_objs, format_str) = Self::load_elements(path)?;
        Ok(FusionBrkPtTableOut {
            header,
            db_infos,
            meta,
            lines: line_objs,
            format_str: format_str.to_string(),
        })
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
        let outname = utils::add_gz_suffix_if_needed(&out);
        let mut mywriter = MyGzWriter::new(&outname)?;

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
        for line in &mut self.lines {
            if line.format_str.is_none() {
                line.update_format_str(&self.format_str);
            }
            let line_str = line.to_table(None, None);
            mywriter.write_all_bytes(line_str.as_bytes())?;
        }
        info!("Saved output to {:?}", outname);
        Ok(())
    }

    fn add_line(&mut self, line: &Line) -> Result<()> {
        self.lines.push(line.clone());
        Ok(())
    }
    #[allow(unused_variables)]
    fn merge(&mut self, other: &Self) -> Result<()> {
        todo!()
    }
}

pub struct FusionDiscoveryTableOut {
    header: Header,
    db_infos: DBInfos,
    meta: Meta,
    lines: Vec<Line>,
    #[allow(dead_code)]
    format_str: String,
}

impl GeneralTableOutput for FusionDiscoveryTableOut {
    fn new(header: Header, db_infos: DBInfos, meta: Meta, format_str: String) -> Self {
        FusionDiscoveryTableOut {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized,
    {
        let (header, db_infos, meta, line_objs, format_str) = Self::load_elements(path)?;
        Ok(FusionDiscoveryTableOut {
            header,
            db_infos,
            meta,
            lines: line_objs,
            format_str: format_str.to_string(),
        })
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
        let outname = utils::add_gz_suffix_if_needed(&out);

        let mut mywriter = MyGzWriter::new(&outname)?;

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
        for line in &mut self.lines {
            if line.format_str.is_none() {
                line.update_format_str(&self.format_str);
            }
            let line_str = line.to_table(None, None);
            mywriter.write_all_bytes(line_str.as_bytes())?;
        }
        info!("Saved output to {:?}", outname);
        Ok(())
    }

    fn add_line(&mut self, line: &Line) -> Result<()> {
        self.lines.push(line.clone());
        Ok(())
    }

    #[allow(unused_variables)]
    fn merge(&mut self, other: &Self) -> Result<()> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_isoform_table_load() {
        let path = "/ssd1/stix-iso-devspace/isopedia-dev/test/test.output.isoform.gz";
        let isoform_table = IsoformTableOut::load(path).unwrap();
        dbg!("header: {:?}", isoform_table.header);
        dbg!("db_infos: {:?}", isoform_table.db_infos);
        dbg!("meta: {:?}", isoform_table.meta);
        dbg!("lines: {:?}", isoform_table.lines.len());
    }
}
