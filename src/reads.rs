use crate::config::PROFILE_SCHEMA_VERSION;
use crate::constants::BUF_SIZE_4M;
use crate::utils::{self, isoform_signature};
use anyhow::{bail, ensure, Context, Result};
use flate2::read::GzDecoder;
use indexmap::IndexSet;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::strand::Strand;

const VERSION_HEADER: &str = "##ISOPEDIA <VERSION> ";
const CHROM_HEADER: &str = "##ISOPEDIA <CHROM> ";
pub const PROFILE_COLUMNS: &str = "#signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n";

/// Append an encoded leaf field; structural separators are written by the caller.
pub fn append_escaped(out: &mut String, value: &str) {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut start = 0;
    for (i, byte) in value.bytes().enumerate() {
        if matches!(
            byte,
            b'\\' | b'"' | b'\t' | b'\r' | b'\n' | b'|' | b':' | b',' | b'%' | b';' | b'='
        ) {
            out.push_str(&value[start..i]);
            out.push_str("\\x");
            out.push(HEX[(byte >> 4) as usize] as char);
            out.push(HEX[(byte & 15) as usize] as char);
            start = i + 1;
        }
    }
    out.push_str(&value[start..]);
}

/// Decode one escaped leaf field, after all enclosing delimiters have been split.
pub fn decode_field(value: &str) -> Result<String> {
    let Some(first) = value.bytes().position(|b| b == b'\\') else {
        return Ok(value.to_owned());
    };
    let bytes = value.as_bytes();
    let mut decoded = Vec::with_capacity(bytes.len());
    decoded.extend_from_slice(&bytes[..first]);
    let mut i = first;
    while i < bytes.len() {
        if bytes[i] == b'\\' {
            ensure!(
                bytes.get(i + 1) == Some(&b'x') && i + 3 < bytes.len(),
                "incomplete or invalid escape at byte {i}"
            );
            let hi = (bytes[i + 2] as char)
                .to_digit(16)
                .context("invalid escape hex digit")?;
            let lo = (bytes[i + 3] as char)
                .to_digit(16)
                .context("invalid escape hex digit")?;
            decoded.push((hi * 16 + lo) as u8);
            i += 4;
        } else {
            decoded.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(decoded).context("decoded field is not UTF-8")
}

pub fn profile_header<'a>(chroms: impl IntoIterator<Item = &'a str>) -> String {
    let mut header = format!(
        "{VERSION_HEADER}version=\"{}\"; program=\"isopedia\"; schema=\"{}\";\n",
        env!("CARGO_PKG_VERSION"),
        PROFILE_SCHEMA_VERSION
    );
    let mut chroms: Vec<_> = chroms.into_iter().collect();
    chroms.sort_unstable();
    chroms.dedup();
    for chrom in chroms {
        header.push_str(CHROM_HEADER);
        header.push_str("ID=\"");
        append_escaped(&mut header, chrom);
        header.push_str("\";\n");
    }
    header.push_str(PROFILE_COLUMNS);
    header
}

fn header_attributes(text: &str) -> Result<FxHashMap<&str, String>> {
    let text = text
        .trim()
        .strip_suffix(';')
        .context("missing metadata terminator ';'")?;
    let mut attributes = FxHashMap::default();
    for attribute in text.split(';') {
        let (key, value) = attribute
            .trim()
            .split_once('=')
            .context("invalid metadata attribute")?;
        let key = key.trim();
        ensure!(!key.is_empty(), "empty metadata key");
        let value = value
            .trim()
            .strip_prefix('"')
            .and_then(|s| s.strip_suffix('"'))
            .context("metadata values must be quoted")?;
        ensure!(!value.contains('"'), "unescaped quote in metadata value");
        let value =
            decode_field(value).with_context(|| format!("invalid metadata value for {key}"))?;
        ensure!(
            attributes.insert(key, value).is_none(),
            "duplicate metadata attribute {key}"
        );
    }
    Ok(attributes)
}

fn parse_profile_strand(value: &str) -> Result<Strand> {
    match value {
        "+" => Ok(Strand::Plus),
        "-" => Ok(Strand::Minus),
        _ => bail!("invalid strand {value:?}"),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Segment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand, // 0 for +, 1 for -
    pub is_supp: bool,
} // using segment instead of single SJ is more informative and save memory(only one chrom and strand for all SJ)

impl Segment {
    pub fn new(chrom: String, start: u64, end: u64, strand: Strand, is_supp: bool) -> Self {
        let chrom2 = utils::trim_chr_prefix_to_upper(&chrom);

        Self {
            chrom: chrom2,
            start,
            end,
            strand,
            is_supp,
        }
    }

    pub fn cmp_end(&self, other: &Self) -> std::cmp::Ordering {
        if self.end < other.end {
            return std::cmp::Ordering::Less;
        } else if self.end > other.end {
            return std::cmp::Ordering::Greater;
        } else {
            return std::cmp::Ordering::Equal;
        }
    }
}

impl PartialEq for Segment {
    fn ne(&self, other: &Self) -> bool {
        if self.chrom == other.chrom
            && self.start == other.start
            && self.end == other.end
            && self.strand == other.strand
            && self.is_supp == other.is_supp
        {
            return false;
        } else {
            return true;
        }
    }

    fn eq(&self, other: &Self) -> bool {
        !self.ne(other)
    }
}

impl PartialOrd for Segment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if !self.is_supp && other.is_supp {
            // make sure the primary alignment is always before the supplementary alignment
            return Some(std::cmp::Ordering::Less);
        } else if self.is_supp && !other.is_supp {
            return Some(std::cmp::Ordering::Greater);
        } else {
            if self.chrom < other.chrom {
                return Some(std::cmp::Ordering::Less);
            } else if self.chrom > other.chrom {
                return Some(std::cmp::Ordering::Greater);
            } else {
                if self.start < other.start {
                    return Some(std::cmp::Ordering::Less);
                } else if self.start > other.start {
                    return Some(std::cmp::Ordering::Greater);
                } else {
                    Some(self.cmp_end(other))
                }
            }
        }
    }
}

#[derive(Debug)]
pub struct SingleRead {
    pub chrom: String,
    pub seg_size: u64,
    pub read_len: u64,
    pub read_start: u64,
    pub right: u64, // not sure if it is necessary
    pub mapq: u8,
    pub evidence: u32,
    pub segment_list: Vec<Segment>,
    pub pri_seg_size: u64,
    pub signature: u128,
    pub strand: Strand,
    pub info: Vec<(String, String)>,
}

impl SingleRead {
    pub fn new(
        chrom: String,
        mapq: u8,
        supp_raw_read: u32,
        readlen: u64,
        start: u64,
        strand: Strand,
    ) -> Self {
        Self {
            chrom: chrom,
            seg_size: 0,
            segment_list: Vec::new(),
            mapq: mapq,
            evidence: supp_raw_read,
            read_len: readlen,
            read_start: start,
            right: 0,
            pri_seg_size: 0,
            signature: 0,
            strand: strand,
            info: Vec::new(),
        }
    }

    pub fn add_info(&mut self, key: String, value: String) {
        self.info.push((key, value));
    }

    pub fn add_segment(
        &mut self,
        chrom: String,
        start: u64,
        end: u64,
        strand: &Strand,
        is_supp: bool,
    ) {
        self.segment_list
            .push(Segment::new(chrom, start, end, *strand, is_supp));
        self.seg_size += 1;
        self.pri_seg_size += 1;
    }

    pub fn add_supp_segment(
        &mut self,
        chrom: String,
        start: u64,
        end: u64,
        strand: &str,
        is_supp: bool,
    ) {
        let strand = Strand::from_string(strand);

        self.segment_list
            .push(Segment::new(chrom, start, end, strand, is_supp));
        self.seg_size += 1
    }

    pub fn update_right(&mut self, right: u64) {
        self.right = right;
    }

    pub fn sort_segments(&mut self) {
        self.segment_list.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    pub fn get_pri_splice_junction(&self) -> Vec<(u64, u64)> {
        let mut res = Vec::new();
        for seg in &self.segment_list {
            if !seg.is_supp {
                res.push(seg.start);
                res.push(seg.end);
            }
        }
        // remove the first and last postion
        if res.len() > 2 {
            res.pop();
            res.remove(0);
        }
        // chunk the vector into pairs for splice junction
        res.chunks(2).into_iter().map(|x| (x[0], x[1])).collect()
    }

    pub fn get_supp_segments(&self) -> Vec<Segment> {
        self.segment_list
            .clone()
            .into_iter()
            .filter(|x| x.is_supp)
            .map(|x| x.clone())
            .collect()
    }

    /// Process the isoform
    /// the signature is the hash value of the splice junctions
    pub fn process(&mut self) {
        self.signature = isoform_signature(
            &self.chrom,
            self.pri_seg_size - 1,
            &self.get_pri_splice_junction(),
        );
        self.sort_segments();
    }
}

#[derive(Debug, Clone)]
pub struct ReadDiff {
    pub left: u64,
    pub right: u64,
    pub strand: Strand,
    pub supp_segs: Vec<Segment>,
    pub info: Vec<(String, String)>,
}

impl ReadDiff {
    pub fn from_string(rec: &str) -> Result<ReadDiff> {
        let mut fields = rec.split(':');
        let (left, right) = fields
            .next()
            .context("missing read span")?
            .split_once('-')
            .context("invalid read span")?;
        let strand = parse_profile_strand(fields.next().context("missing read strand")?)?;
        let mut supp_segs = Vec::new();
        for segment in fields
            .next()
            .context("missing supplementary segments")?
            .split(',')
            .filter(|s| !s.is_empty())
        {
            let mut parts = segment.split('%');
            let chrom = decode_field(parts.next().context("missing SA chromosome")?)
                .context("invalid SA chromosome")?;
            let start = parts.next().context("missing SA start")?.parse()?;
            let end = parts.next().context("missing SA end")?.parse()?;
            let strand = parse_profile_strand(parts.next().context("missing SA strand")?)?;
            ensure!(parts.next().is_none(), "extra supplementary segment fields");
            supp_segs.push(Segment {
                chrom,
                start,
                end,
                strand,
                is_supp: true,
            });
        }
        let mut info = Vec::new();
        for entry in fields
            .next()
            .unwrap_or("")
            .split(';')
            .filter(|s| !s.is_empty())
        {
            let (key, value) = entry.split_once('=').context("invalid info entry")?;
            info.push((
                decode_field(key).context("invalid info key")?,
                decode_field(value).with_context(|| format!("invalid info value for {key}"))?,
            ));
        }
        ensure!(
            fields.next().is_none(),
            "extra read fields (unescaped colon)"
        );
        Ok(ReadDiff {
            left: left.parse()?,
            right: right.parse()?,
            strand,
            supp_segs,
            info,
        })
    }

    fn append_to(&self, out: &mut String) {
        let _ = write!(out, "{}-{}:{}:", self.left, self.right, self.strand);
        for (i, seg) in self.supp_segs.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            append_escaped(out, &seg.chrom);
            let _ = write!(out, "%{}%{}%{}", seg.start, seg.end, seg.strand);
        }
        out.push(':');
        for (i, (key, value)) in self.info.iter().enumerate() {
            if i > 0 {
                out.push(';');
            }
            append_escaped(out, key);
            out.push('=');
            append_escaped(out, value);
        }
    }
}

/// Aggregated isoform from single sample in extractor
/// also used in aggregator
#[derive(Debug, Clone)]
pub struct AggrRead {
    pub signature: u128,
    pub chrom: String,
    pub splice_junctions: Vec<(u64, u64)>,
    pub evidence: u32,
    pub read_diffs: Vec<ReadDiff>,
}

impl AggrRead {
    pub fn new(sr: SingleRead) -> AggrRead {
        AggrRead {
            signature: sr.signature,
            chrom: sr.chrom.clone(),
            splice_junctions: sr.get_pri_splice_junction(),
            evidence: 1,
            read_diffs: vec![ReadDiff {
                left: sr.read_start,
                right: sr.right,
                strand: sr.strand,
                supp_segs: sr.get_supp_segments(),
                info: sr.info.clone(),
            }],
        }
    }

    pub fn add_read(&mut self, isoform: SingleRead) {
        self.evidence += 1;
        self.read_diffs.push(ReadDiff {
            left: isoform.read_start,
            right: isoform.right,
            strand: isoform.strand,
            supp_segs: isoform.get_supp_segments(),
            info: isoform.info.clone(),
        });
    }

    /// Serialize an escaped record. Hashes and in-memory fields remain unescaped.
    pub fn to_record(&self) -> String {
        let mut record = String::new();
        let _ = write!(record, "{:032x}\t{}\t", self.signature, self.evidence);
        append_escaped(&mut record, &self.chrom);
        record.push('\t');
        for (i, (left, right)) in self.splice_junctions.iter().enumerate() {
            if i > 0 {
                record.push(',');
            }
            let _ = write!(record, "{left}-{right}");
        }
        record.push('\t');
        for (i, diff) in self.read_diffs.iter().enumerate() {
            if i > 0 {
                record.push('|');
            }
            diff.append_to(&mut record);
        }
        record.push('\n');
        record
    }

    pub fn from_string(rec: &str) -> Result<AggrRead> {
        let mut fields = rec.trim_end_matches(['\r', '\n']).split('\t');
        let field = fields.next().context("missing signature")?;
        ensure!(
            field.len() == 32 && field.bytes().all(|b| b.is_ascii_hexdigit()),
            "signature must contain exactly 32 hexadecimal characters"
        );
        let signature = u128::from_str_radix(field, 16)?;
        let evidence = fields.next().context("missing evidence")?.parse()?;
        let chrom = decode_field(fields.next().context("missing chromosome")?)
            .context("invalid chromosome")?;
        let mut splice_junctions = Vec::new();
        for junction in fields
            .next()
            .context("missing junctions")?
            .split(',')
            .filter(|s| !s.is_empty())
        {
            let (left, right) = junction.split_once('-').context("invalid junction")?;
            splice_junctions.push((left.parse()?, right.parse()?));
        }
        let read_diffs = fields
            .next()
            .context("missing read differences")?
            .split('|')
            .filter(|s| !s.is_empty())
            .map(|s| ReadDiff::from_string(s))
            .collect::<Result<Vec<_>>>()?;
        ensure!(fields.next().is_none(), "extra profile columns");
        Ok(AggrRead {
            signature,
            evidence,
            chrom,
            splice_junctions,
            read_diffs,
        })
    }

    pub fn merge(&mut self, other: AggrRead) {
        self.evidence += other.evidence;
        self.read_diffs.extend(other.read_diffs);
    }
}

pub struct SingleSampleReader {
    pub file_name: String,
    pub chroms: IndexSet<String>,
    reader: BufReader<GzDecoder<File>>,
    curr_rec_string: String,
    line_number: u64,
}

impl SingleSampleReader {
    /// Read the annotated header and check the configured schema once per file.
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("open profile {}", path.display()))?;
        let mut reader = BufReader::with_capacity(BUF_SIZE_4M, GzDecoder::new(file));
        let mut line = String::new();
        let mut line_number = 1;
        let header = (|| -> Result<IndexSet<String>> {
            ensure!(reader.read_line(&mut line)? > 0, "missing VERSION header");
            let attributes = header_attributes(
                line.trim_end_matches(['\r', '\n'])
                    .strip_prefix(VERSION_HEADER)
                    .context("missing valid VERSION header")?,
            )?;
            for key in ["version", "program", "schema"] {
                ensure!(
                    attributes.get(key).is_some_and(|v| !v.is_empty()),
                    "missing {key} attribute"
                );
            }
            let schema: u32 = attributes["schema"]
                .parse()
                .context("schema must be a u32 integer")?;
            ensure!(
                schema == PROFILE_SCHEMA_VERSION,
                "unsupported profile schema {schema}; expected {PROFILE_SCHEMA_VERSION}"
            );
            let mut chroms = IndexSet::new();
            loop {
                line.clear();
                line_number += 1;
                ensure!(reader.read_line(&mut line)? > 0, "missing column header");
                let text = line.trim_end_matches(['\r', '\n']);
                if text == PROFILE_COLUMNS.trim_end() {
                    return Ok(chroms);
                }
                ensure!(
                    !text.starts_with("##ISOPEDIA <VERSION>"),
                    "duplicate VERSION header"
                );
                if text.starts_with("##ISOPEDIA <CHROM>") {
                    let mut attributes = header_attributes(
                        text.strip_prefix(CHROM_HEADER)
                            .context("invalid CHROM header")?,
                    )?;
                    let chrom = attributes.remove("ID").context("missing CHROM ID")?;
                    ensure!(!chrom.is_empty(), "empty CHROM ID");
                    ensure!(chroms.insert(chrom.clone()), "duplicate CHROM ID {chrom:?}");
                } else {
                    ensure!(text.starts_with("##"), "invalid or missing column header");
                }
            }
        })()
        .with_context(|| format!("parse {} line {}", path.display(), line_number))?;
        Ok(Self {
            file_name: path.display().to_string(),
            chroms: header,
            reader,
            curr_rec_string: String::new(),
            line_number,
        })
    }

    pub fn next_rec(&mut self) -> Result<Option<AggrRead>> {
        self.curr_rec_string.clear();
        let count = self
            .reader
            .read_line(&mut self.curr_rec_string)
            .with_context(|| format!("read {} line {}", self.file_name, self.line_number + 1))?;
        if count == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        let record = (|| -> Result<AggrRead> {
            ensure!(
                !self.curr_rec_string.starts_with('#'),
                "header line in record section"
            );
            let record = AggrRead::from_string(&self.curr_rec_string)?;
            ensure!(
                self.chroms.contains(&record.chrom),
                "undeclared chromosome {:?}",
                record.chrom
            );
            Ok(record)
        })()
        .with_context(|| format!("parse {} line {}", self.file_name, self.line_number))?;
        Ok(Some(record))
    }
}

#[cfg(test)]
#[path = "tests/reads.rs"]
mod tests;
