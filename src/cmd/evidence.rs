use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::PathBuf,
};

use std::io::Write;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::{bufread::GzDecoder, write::GzEncoder, Compression};
use log::error;
use serde::Serialize;

use crate::{
    constants::MERGED_FILE_NAME,
    ptir::PTIR,
    tools::ToolCmdValidate,
    utils::{greetings2, log_wall_time, start_wall_timer},
};

/// Written by isopedia >= 1.8 merge; its PTIR layout is incompatible with 1.x.
const ARCHIVE_INFO_FILE: &str = "archive_info.json";

#[derive(Parser, Debug, Serialize, Clone)]
#[command(name = "isopedia call")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[Call] Call transcripts from indexed catalog.
", long_about = None)]
#[clap(after_long_help = "")]

pub struct EvidenceCli {
    /// Path to the index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// Minimal sample required for output evidence
    #[arg(short, long, default_value_t = 2)]
    pub min_sample: usize,

    /// Minimum number of reads required to define a positive sample
    #[arg(short = 'r', long, default_value_t = 2)]
    pub min_read: u32,

    /// Output prefix; writes <output>.evidence.gtf.gz
    #[arg(short, long)]
    pub output: PathBuf,

    /// Region to output
    #[arg(short='R', long)]
    pub region: String,
}

impl ToolCmdValidate for EvidenceCli {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.idxdir.exists() {
            error!(
                "--idxdir: index directory {} does not exist",
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join(MERGED_FILE_NAME).exists() {
            error!(
                "--idxdir: merged isoform data {} does not exist in {}, please run `isopedia merge` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if self.min_read == 0 {
            error!("--min-read: must be greater than 0");
            is_ok = false;
        }

        if self.idxdir.join(ARCHIVE_INFO_FILE).exists() {
            error!(
                "--idxdir: {} was built by isopedia >= 1.8 and cannot be read by {}; use isopedia >= 1.8 or rebuild the index with this version",
                self.idxdir.display(),
                env!("CARGO_PKG_VERSION")
            );
            is_ok = false;
        }

        is_ok
    }
}

pub fn run_evidence(cli: &EvidenceCli) -> Result<()> {
    let started = start_wall_timer();
    greetings2(&cli);
    if !cli.validate() {
        bail!("invalid arguments, see errors above");
    }

    let ptir_reader = SequentialArchiveReader::load_from_disk(&cli.idxdir.join(MERGED_FILE_NAME))?;

    let mut output_path = cli.output.as_os_str().to_os_string();
    output_path.push(".evidence.gz");
    let mut writer = GzEncoder::new(
        BufWriter::new(File::create(&output_path)?),
        Compression::default(),
    );

    for ptir in ptir_reader {
        let ptir = ptir?;

        if filter_ptir(&ptir, cli.min_read, cli.min_sample) {
            ptir_to_gtf_records(&ptir, &mut writer)?
        }
    }

    writer.finish()?.flush()?;
    log_wall_time("call", started);

    Ok(())
}

struct SequentialArchiveReader {
    decoder: GzDecoder<Box<dyn BufRead + Send>>,
    decoded: Vec<u8>,
    needs_reset: bool,
    finished: bool,
}

impl SequentialArchiveReader {
    pub fn load_from_disk(path: &PathBuf) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mut reader = BufReader::new(file);
        let finished = reader.fill_buf()?.is_empty();
        Ok(Self {
            decoder: GzDecoder::new(Box::new(reader)),
            decoded: Vec::new(),
            needs_reset: false,
            finished,
        })
    }
}

impl Iterator for SequentialArchiveReader {
    type Item = Result<PTIR>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let result = (|| -> Result<Option<PTIR>> {
            if self.needs_reset {
                if self.decoder.get_mut().fill_buf()?.is_empty() {
                    return Ok(None);
                }
                // Preserve buffered input while resetting the completed gzip member.
                let reader = std::mem::replace(self.decoder.get_mut(), Box::new(std::io::empty()));
                self.decoder.reset(reader);
                self.needs_reset = false;
            }
            self.decoded.clear();
            self.decoder.read_to_end(&mut self.decoded)?;
            let ptir: PTIR =
                bincode::deserialize(&self.decoded).context("deserialize sequential PTIR")?;
            self.needs_reset = true;
            Ok(Some(ptir))
        })();
        self.finished = !matches!(&result, Ok(Some(_)));
        result.transpose()
    }
}

/// Keep a PTIR when at least `min_sample` samples each have at least `min_read` reads.
fn filter_ptir(ptir: &PTIR, min_read: u32, min_sample: usize) -> bool {
    ptir.get_positive_count(&min_read) as usize >= min_sample
}

#[test]
fn filter_counts_samples_meeting_min_read() {
    let ptir = PTIR {
        signature: 0,
        total_evidence: 6,
        sample_size: 3,
        sample_evidence_arr: vec![2, 1, 3],
        sample_offset_arr: vec![0, 2, 3],
        chrom: "1".into(),
        chrom_id: 0,
        rec_type: crate::ptir::RecordType::Spliced,
        splice_junctions_vec: vec![],
        isoform_reads_slim_vec: vec![],
        supp_segs_vec: vec![],
    };
    assert!(filter_ptir(&ptir, 1, 3));
    assert!(filter_ptir(&ptir, 2, 2));
    assert!(!filter_ptir(&ptir, 2, 3));
    assert!(filter_ptir(&ptir, 3, 1));
    assert!(!filter_ptir(&ptir, 4, 1));
}

fn ptir_to_gtf_records(ptir: &PTIR, writer: &mut dyn Write) -> Result<()> {
    let first = ptir
        .isoform_reads_slim_vec
        .first()
        .with_context(|| format!("PTIR on {} has no reads", ptir.chrom))?;
    let (mut left, mut right) = (first.left, first.right);
    let mut strand = match first.strand {
        crate::strand::Strand::Plus => "+",
        crate::strand::Strand::Minus => "-",
    };
    for read in &ptir.isoform_reads_slim_vec {
        if read.left >= read.right {
            bail!(
                "invalid PTIR read span on {}: {}-{}",
                ptir.chrom,
                read.left,
                read.right
            );
        }
        left = left.min(read.left);
        right = right.max(read.right);
        // ponytail: mixed read strands remain unknown until strand inference is implemented.
        if read.strand != first.strand {
            strand = ".";
        }
    }

    let mut exons = Vec::new();
    let mut start = left;
    if !ptir.is_mono_exonic() {
        for &(donor, acceptor) in &ptir.splice_junctions_vec {
            if start > donor || donor >= acceptor || acceptor > right {
                bail!(
                    "invalid PTIR splice junction on {}: {donor}-{acceptor}",
                    ptir.chrom
                );
            }
            if start == donor {
                log::warn!(
                    "PTIR on {}:{}-{}: skipping zero-length exon at {}; continuing GTF output",
                    ptir.chrom,
                    left,
                    right,
                    start
                );
            } else {
                exons.push((start, donor));
            }
            start = acceptor;
        }
    }
    if start == right {
        log::warn!(
            "PTIR on {}:{}-{}: skipping zero-length terminal exon at {}; continuing GTF output",
            ptir.chrom,
            left,
            right,
            start
        );
    } else {
        exons.push((start, right));
    }

    // PTIR coordinates are 0-based, half-open; GTF coordinates are 1-based, inclusive.
    // ponytail: use the PTIR signature as a placeholder gene ID until gene grouping exists.
    let id = ptir.signature;
    writeln!(writer,
        "{}\tisopedia\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"isopedia.{id}\"; transcript_id \"isopedia.{id}\";",
        ptir.chrom, left + 1, right, strand
    )?;
    for (index, &(start, end)) in exons.iter().enumerate() {
        let exon_number = if strand == "-" {
            exons.len() - index
        } else {
            index + 1
        };
        writeln!(writer,
            "{}\tisopedia\texon\t{}\t{}\t.\t{}\t.\tgene_id \"isopedia.{id}\"; transcript_id \"isopedia.{id}\"; exon_number \"{exon_number}\";",
            ptir.chrom, start + 1, end, strand
        )?;
    }
    Ok(())
}

#[test]
fn gtf_skips_zero_length_exon_and_preserves_strand_numbering() {
    let mut ptir = PTIR {
        signature: 0,
        total_evidence: 1,
        sample_size: 1,
        sample_evidence_arr: vec![1],
        sample_offset_arr: vec![0],
        chrom: "1".into(),
        chrom_id: 0,
        rec_type: crate::ptir::RecordType::Spliced,
        splice_junctions_vec: vec![(9554665, 9554954), (9554954, 9557217)],
        isoform_reads_slim_vec: vec![crate::ptir::ReadDiffSlim {
            left: 9554545,
            right: 9557902,
            strand: crate::strand::Strand::Plus,
            supp_seg_vec_offset: 0,
            supp_seg_vec_length: 0,
            info: vec![],
        }],
        supp_segs_vec: vec![],
    };
    for strand in [crate::strand::Strand::Plus, crate::strand::Strand::Minus] {
        ptir.isoform_reads_slim_vec[0].strand = strand;
        let mut output = Vec::new();
        ptir_to_gtf_records(&ptir, &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();
        let rows: Vec<Vec<&str>> = text
            .lines()
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(rows.len(), 3);
        assert_eq!(&rows[0][2..5], &["transcript", "9554546", "9557902"]);
        assert_eq!(&rows[1][2..5], &["exon", "9554546", "9554665"]);
        assert_eq!(&rows[2][2..5], &["exon", "9557218", "9557902"]);
        for (index, row) in rows.iter().enumerate() {
            assert_eq!(row[6], strand.to_string());
            if index > 0 {
                let number = if strand == crate::strand::Strand::Plus {
                    index
                } else {
                    3 - index
                };
                assert!(row[8].ends_with(&format!("exon_number \"{number}\";")));
            }
        }
    }
    for junction in [(9554545, 9557217), (9554665, 9557902)] {
        ptir.splice_junctions_vec = vec![junction];
        let mut output = Vec::new();
        ptir_to_gtf_records(&ptir, &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();
        let rows: Vec<Vec<&str>> = text
            .lines()
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1][2], "exon");
        assert!(rows[1][3].parse::<u64>().unwrap() <= rows[1][4].parse::<u64>().unwrap());
        assert!(rows[1][8].ends_with("exon_number \"1\";"));
    }
    ptir.splice_junctions_vec = vec![(9554665, 9554954), (9554953, 9557217)];
    let mut output = Vec::new();
    assert!(ptir_to_gtf_records(&ptir, &mut output).is_err());
    assert!(output.is_empty());
}

#[cfg(test)]
mod sequential_reader_tests {
    use super::*;

    #[test]
    fn reuses_decoder_and_stops_at_eof_or_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(MERGED_FILE_NAME);
        let mut record = PTIR {
            signature: 0,
            total_evidence: 0,
            sample_size: 0,
            sample_evidence_arr: vec![],
            sample_offset_arr: vec![],
            chrom: "chr1".to_string(),
            chrom_id: 0,
            rec_type: crate::ptir::RecordType::MonoExon,
            splice_junctions_vec: vec![],
            isoform_reads_slim_vec: vec![],
            supp_segs_vec: vec![],
        };
        let mut first = Vec::new();
        record.gz_encode(&mut first);
        record.chrom = "chr2".to_string();
        let mut second = Vec::new();
        record.gz_encode(&mut second);

        std::fs::write(&path, [&first[..], &second[..]].concat()).unwrap();
        let mut reader = SequentialArchiveReader::load_from_disk(&path).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().chrom, "chr1");
        assert_eq!(reader.next().unwrap().unwrap().chrom, "chr2");
        assert!(reader.next().is_none());
        assert!(reader.next().is_none());

        std::fs::write(&path, []).unwrap();
        let mut reader = SequentialArchiveReader::load_from_disk(&path).unwrap();
        assert!(reader.next().is_none());

        let mut invalid_ptir = Vec::new();
        flate2::bufread::GzEncoder::new(&b"x"[..], flate2::Compression::default())
            .read_to_end(&mut invalid_ptir)
            .unwrap();
        let mut bad_crc = second.clone();
        let crc_start = bad_crc.len() - 8;
        bad_crc[crc_start] ^= 1;
        for (bad_member, codec_error) in [
            (second[..second.len() - 1].to_vec(), false),
            (second[..10].to_vec(), false),
            (b"invalid gzip".to_vec(), false),
            (bad_crc, false),
            (invalid_ptir, true),
        ] {
            for prefix in [&[][..], first.as_slice()] {
                std::fs::write(&path, [prefix, &bad_member].concat()).unwrap();
                let mut reader = SequentialArchiveReader::load_from_disk(&path).unwrap();
                if !prefix.is_empty() {
                    assert_eq!(reader.next().unwrap().unwrap().chrom, "chr1");
                }
                let error = reader.next().unwrap().unwrap_err();
                if codec_error {
                    assert!(error.downcast_ref::<bincode::Error>().is_some());
                } else {
                    assert!(error.downcast_ref::<std::io::Error>().is_some());
                }
                assert!(reader.next().is_none());
                assert!(reader.next().is_none());
            }
        }
    }
}
