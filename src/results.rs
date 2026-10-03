use crate::cmd::isoform::ISOFORM_FORMAT;
use crate::constants::FORMAT_STR_NAME;
use crate::meta::Meta;
use crate::myio::MyGzWriter;
use crate::myio::{DBInfos, GeneralOutputIO, Header, Line, MyGzReader};
// use crate::output_traits::GeneralTableOutputTrait;
use crate::utils;
use anyhow::Result;
use log::info;
use std::path::Path;

pub struct TableOutput {
    header: Header,
    db_infos: DBInfos,
    pub meta: Meta,
    lines: Vec<Line>,
    pub format_str: String,
    writer: Option<MyGzWriter>,
}

impl TableOutput {
    pub fn get_mem_size(&self) -> usize {
        let mut total = 0usize;

        for line in &self.lines {
            total += std::mem::size_of_val(line);

            total += line.field_vec.capacity() * std::mem::size_of::<String>();
            for s in &line.field_vec {
                total += s.capacity();
            }

            total += line.sample_vec.capacity() * std::mem::size_of::<crate::myio::SampleChip>();
            for sample in &line.sample_vec {
                total += sample.init_string.capacity() * std::mem::size_of::<String>();
                for s in &sample.init_string {
                    total += s.capacity();
                }
            }
        }

        total
    }

    pub fn new<P: AsRef<Path>>(
        path: P,
        header: Header,
        db_infos: DBInfos,
        meta: Meta,
        format_str: String,
    ) -> Self {
        let outname = utils::add_gz_suffix_if_needed(&path);
        let mut writer = MyGzWriter::new(&outname).ok();

        // dump header dbinfos meta to file first
        let header_t = header.to_table(Some("#"), None);
        let dbinfo_t = db_infos.to_table(Some("##[DBINFO]"), None);
        let meta_t = meta.to_table(Some("##[SAMPLE]"), None);

        if let Some(w) = &mut writer {
            w.write_all_bytes(meta_t.as_bytes()).unwrap();
            w.write_all_bytes(dbinfo_t.as_bytes()).unwrap();
            w.write_all_bytes(header_t.as_bytes()).unwrap();
        }

        TableOutput {
            header,
            db_infos,
            meta,
            lines: Vec::new(),
            format_str,
            writer,
        }
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self>
    where
        Self: Sized,
    {
        let (header, db_infos, meta, line_objs, format_str) = Self::load_elements(path)?;

        Ok(TableOutput {
            header,
            db_infos,
            meta,
            lines: line_objs,
            format_str: format_str.to_string(),
            writer: None,
        })
    }

    pub fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
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

    pub fn add_line(&mut self, line: &mut Line) -> Result<()> {
        // if format string is not set, set it
        if line.format_str.is_none() {
            line.update_format_str(&self.format_str);
        }

        self.lines.push(line.clone());

        if self.lines.len() > 0 && self.lines.len() % 1_000 == 0 {
            // info!("Flushing {} lines to file...", self.lines.len());
            if let Some(w) = &mut self.writer {
                // write lines
                for each_line in &mut self.lines {
                    if each_line.format_str.is_none() {
                        each_line.update_format_str(&self.format_str);
                    }
                    let line_str = each_line.to_table(None, None);

                    w.write_all_bytes(line_str.as_bytes())?;
                }

                w.flush()?;
            }

            self.lines.clear();
            self.lines.shrink_to_fit();
        }

        Ok(())
    }

    pub fn write_bytes(&mut self, bytes: &[u8]) -> Result<()> {
        if let Some(w) = &mut self.writer {
            w.write_all_bytes(bytes)?;
        }

        Ok(())
    }

    pub fn write_format_str(&mut self) -> Result<()> {
        if let Some(w) = &mut self.writer {
            w.write_all_bytes(self.format_str.as_bytes())?;
        }

        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        // write remaining lines
        if let Some(w) = &mut self.writer {
            if self.lines.len() > 0 {
                for each_line in &mut self.lines {
                    if each_line.format_str.is_none() {
                        each_line.update_format_str(&self.format_str);
                    }
                    let line_str = each_line.to_table(None, None);

                    w.write_all_bytes(line_str.as_bytes())?;
                }
            }

            w.flush()?;
        }

        self.lines.clear();
        self.lines.shrink_to_fit();

        // info!(
        //     "Saved output to file {:?}",
        //     self.writer.as_ref().unwrap().path()
        // );

        Ok(())
    }

    pub fn get_out_path(&self) -> Option<String> {
        self.writer.as_ref().map(|w| w.path().to_string())
    }

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

            if line_buf.starts_with("##[SAMPLE]") {
                meta_strings.push_str(&line_buf);
            } else if line_buf.starts_with("##[DBINFO]") {
                db_infos_strings.push_str(&line_buf);
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

        Ok((header, db_infos, meta, line_objs, format_str.to_string()))
    }
}

#[cfg(test)]
mod isoform_merge_tests {
    use super::*;
    use crate::myio::SampleChip;

    fn table(sample_name: &str, direct: u64, em: f32) -> TableOutput {
        let mut header = Header::new();
        for column in [
            "chrom",
            "start",
            "end",
            "length",
            "exon_count",
            "trans_id",
            "gene_id",
            "ranking_score",
            "detected",
            "min_read",
            "n_pos_samples/sample_size",
            "attributes",
        ] {
            header.add_column(column).unwrap();
        }
        header.add_sample_name(sample_name).unwrap();
        let mut db_infos = DBInfos::new();
        db_infos.add_sample_evidence(sample_name, 100);
        let line = Line {
            field_vec: [
                "chr1", "100", "500", "400", "2", "tx", "gene", "0", "no", "1", "0/1", "attrs",
            ]
            .iter()
            .map(|x| x.to_string())
            .collect(),
            format_str: Some(ISOFORM_FORMAT.to_string()),
            sample_vec: vec![SampleChip {
                sample_name: None,
                init_string: vec![
                    direct.to_string(),
                    "0".into(),
                    "0".into(),
                    "0".into(),
                    direct.to_string(),
                    "0".into(),
                    em.to_string(),
                    "0".into(),
                    "0".into(),
                    "0".into(),
                    "0".into(),
                    "0".into(),
                    "0".into(),
                ],
            }],
        };
        TableOutput {
            header,
            db_infos,
            meta: Meta::new_empty(vec![sample_name.to_string()]),
            lines: vec![line],
            format_str: ISOFORM_FORMAT.to_string(),
            writer: None,
        }
    }

    #[test]
    fn merge_isoform_uses_direct_support_not_em() {
        let mut left = table("s1", 1, 0.0);
        let right = table("s2", 0, 5.0);
        left.merge_isoform(&right).unwrap();
        assert_eq!(left.lines[0].field_vec[8], "yes");
        assert_eq!(left.lines[0].field_vec[10], "1/2");
        assert_eq!(left.lines[0].sample_vec.len(), 2);
        assert_eq!(left.lines[0].sample_vec[0].init_string[4], "1");
        assert_eq!(left.lines[0].sample_vec[1].init_string[5], "0");
        assert_eq!(left.lines[0].field_vec[7], "NA");
    }

    #[test]
    fn merge_isoform_rejects_old_order_or_short_format_before_mutation() {
        let mut left = table("s1", 1, 0.0);
        let mut old = table("s2", 1, 0.0);
        old.format_str = "RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES".to_string();
        old.lines[0].format_str = Some(old.format_str.clone());
        old.lines[0].sample_vec[0].init_string.truncate(11);
        assert!(left.merge_isoform(&old).is_err());
        assert_eq!(left.lines[0].sample_vec.len(), 1);

        let mut previous_13 = table("s2", 1, 0.0);
        previous_13.format_str = "RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES:RC_FSM_JC_EXACT:RC_FSM_JC_WOBBLE_ONLY".to_string();
        previous_13.lines[0].format_str = Some(previous_13.format_str.clone());
        assert!(left.merge_isoform(&previous_13).is_err());
        assert_eq!(left.lines[0].sample_vec.len(), 1);

        let mut short = table("s2", 1, 0.0);
        short.lines[0].sample_vec[0].init_string.pop();
        assert!(left.merge_isoform(&short).is_err());
        assert_eq!(left.lines[0].sample_vec.len(), 1);
    }

    #[test]
    fn merge_isoform_rejects_duplicate_samples_without_mutation() {
        let mut left = table("s1", 1, 0.0);
        let right = table("s1", 2, 0.0);
        assert!(left.merge_isoform(&right).is_err());
        assert_eq!(left.header.sample_names, vec!["s1"]);
        assert_eq!(left.db_infos.sample_total_evidence_map.len(), 1);
        assert_eq!(left.lines[0].sample_vec.len(), 1);
    }

    #[test]
    fn new_isoform_format_round_trips_through_table_loader() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("isoform.tsv.gz");
        table("s1", 2, 1.5).save_to_file(&path).unwrap();
        let loaded = TableOutput::load(&path).unwrap();
        assert_eq!(loaded.format_str, ISOFORM_FORMAT);
        assert_eq!(loaded.format_str.split(':').count(), 13);
        assert_eq!(loaded.lines[0].sample_vec[0].init_string.len(), 13);
        assert_eq!(loaded.lines[0].sample_vec[0].init_string[0], "2");
        assert_eq!(loaded.lines[0].sample_vec[0].init_string[4], "2");
        assert_eq!(loaded.lines[0].sample_vec[0].init_string[5], "0");
        assert_eq!(loaded.lines[0].sample_vec[0].init_string[6], "1.5");
    }
}

impl TableOutput {
    pub fn merge_isoform(&mut self, other: &Self) -> Result<()> {
        if self.format_str != ISOFORM_FORMAT || other.format_str != ISOFORM_FORMAT {
            return Err(anyhow::anyhow!(
                "Cannot merge isoform tables with old or different FORMAT fields"
            ));
        }
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

        let mut merged_counts = Vec::with_capacity(self.lines.len());
        for (i, (line, other_line)) in self.lines.iter().zip(&other.lines).enumerate() {
            if line.field_vec.len() != 12
                || other_line.field_vec.len() != 12
                || line.sample_vec.len() != self.db_infos.sample_total_evidence_map.len()
                || other_line.sample_vec.len() != other.db_infos.sample_total_evidence_map.len()
            {
                return Err(anyhow::anyhow!(
                    "Invalid isoform fields or sample count at line {}",
                    i + 1
                ));
            }
            if (0..=6)
                .chain(std::iter::once(11))
                .any(|j| line.field_vec[j] != other_line.field_vec[j])
            {
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
            let min_read = line.field_vec[9].parse::<u64>()?;
            if min_read == 0 {
                return Err(anyhow::anyhow!("Invalid min_read=0 at line {}", i + 1));
            }
            let mut counts =
                Vec::with_capacity(line.sample_vec.len() + other_line.sample_vec.len());
            for sample in line.sample_vec.iter().chain(&other_line.sample_vec) {
                if sample.init_string.len() != 13 {
                    return Err(anyhow::anyhow!(
                        "Expected 13 sample fields at line {}",
                        i + 1
                    ));
                }
                counts.push(sample.init_string[0].parse::<u64>()?);
            }
            merged_counts.push((counts, min_read));
        }

        let mut merged_meta = self.meta.clone();
        let mut merged_header = self.header.clone();
        let mut merged_db_infos = self.db_infos.clone();
        merged_meta.merge(&other.meta)?;
        merged_header.merge(&other.header)?;
        merged_db_infos.merge(&other.db_infos)?;
        self.meta = merged_meta;
        self.header = merged_header;
        self.db_infos = merged_db_infos;
        for ((line, other_line), (counts, min_read)) in
            self.lines.iter_mut().zip(&other.lines).zip(merged_counts)
        {
            let positive = counts.iter().filter(|&&count| count >= min_read).count();
            line.field_vec[7] = "NA".to_string();
            line.field_vec[8] = if positive > 0 { "yes" } else { "no" }.to_string();
            line.field_vec[10] = format!("{}/{}", positive, counts.len());
            line.sample_vec.extend(other_line.sample_vec.clone());
        }
        info!("Finished");
        Ok(())
    }
}
