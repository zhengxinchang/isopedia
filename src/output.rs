use crate::io::{DBInfos, GeneralOutputIO, Header, Line};
use crate::meta::Meta;
use crate::writer::MyGzWriter;
use anyhow::Result;
use std::path::Path;

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
    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()>;
    fn add_line(&mut self, line: &Line) -> Result<()>;

    fn merge(&mut self, other: &Self) -> Result<()>;
}

pub struct IsoformTableOut {
    header: Header,
    db_infos: DBInfos,
    meta: Meta,
    lines: Vec<Line>,
    #[allow(dead_code)]
    format_str: String,
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

    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
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
        todo!()
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

    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
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

    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
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

    fn load<P: AsRef<Path>>(path: P) -> Self
    where
        Self: Sized,
    {
        todo!()
    }

    fn save_to_file<P: AsRef<Path>>(&mut self, out: P) -> Result<()> {
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
        todo!()
    }
}
