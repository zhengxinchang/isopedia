use std::{
    cmp::{Ordering, Reverse},
    collections::{BTreeMap, BinaryHeap},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::PathBuf,
};

use std::io::Write;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::{bufread::GzDecoder, write::GzEncoder, Compression};
use log::{error, warn};
use rustc_hash::FxHashSet;
use serde::Serialize;

use crate::{
    bptree::BPTree,
    chromosome::ChromMapping,
    constants::{CHROM_FILE_NAME, DATASET_INFO_FILE_NAME, MERGED_FILE_NAME},
    dataset_info::DatasetInfo,
    ptir::PTIR,
    ptir_archive::PTIRArchiveCache,
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

    /// Only output PTIRs with a splice site or read end in this region: chr or chr:start-end (1-based, inclusive).
    /// PTIRs spanning the whole region without such a position inside it are not reported.
    /// Requires `isopedia index`; without --region the whole merged archive is scanned.
    #[arg(short = 'R', long)]
    pub region: Option<String>,

    /// BED file of regions (chrom, start, end; 0-based, half-open), each queried like --region.
    /// Overlapping regions are merged so every PTIR is reported once, ordered by chromosome and position.
    #[arg(long, conflicts_with = "region")]
    pub region_bed: Option<PathBuf>,

    /// Also write every supporting read to <output>.evidence.src.tsv.gz
    /// (transcript_id, sample_id, sample_name, start, end, strand; 1-based, inclusive).
    /// One row per read: can be very large without --region/--region-bed.
    #[arg(long, default_value_t = false)]
    pub src: bool,
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

        for file in [MERGED_FILE_NAME, DATASET_INFO_FILE_NAME] {
            if !self.idxdir.join(file).exists() {
                error!(
                    "--idxdir: {} does not exist in {}, please run `isopedia merge` first",
                    file,
                    self.idxdir.display()
                );
                is_ok = false;
            }
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

        if let Some(region) = &self.region {
            if let Err(e) = parse_region(region) {
                error!("--region: {e}");
                is_ok = false;
            }
        }

        if let Some(bed) = &self.region_bed {
            if !bed.exists() {
                error!("--region-bed: {} does not exist", bed.display());
                is_ok = false;
            }
        }

        if self.region.is_some() || self.region_bed.is_some() {
            for file in [CHROM_FILE_NAME, "bptree_0.idx"] {
                if !self.idxdir.join(file).exists() {
                    error!(
                        "--region/--region-bed: {} does not exist in {}, please run `isopedia index` first",
                        file,
                        self.idxdir.display()
                    );
                    is_ok = false;
                }
            }
        }

        is_ok
    }
}

/// Parse `chr` or `chr:start-end` (1-based, inclusive) into a chromosome and B+ tree key range.
/// Keys are boundary coordinates (0-based starts, half-open ends), so [start, end] maps to [start - 1, end].
fn parse_region(region: &str) -> Result<(&str, u64, u64)> {
    let Some((chrom, range)) = region.rsplit_once(':') else {
        if region.is_empty() {
            bail!("empty region");
        }
        return Ok((region, 0, u64::MAX));
    };
    let (start, end) = range
        .split_once('-')
        .with_context(|| format!("{region}: expected chr:start-end"))?;
    let (start, end): (u64, u64) = (start.parse()?, end.parse()?);
    if chrom.is_empty() || start == 0 || start > end {
        bail!("{region}: expected chr:start-end with 1 <= start <= end");
    }
    Ok((chrom, start - 1, end))
}

/// Parse BED lines (chrom, start, end; 0-based, half-open) into B+ tree key ranges [start, end].
fn parse_bed(path: &PathBuf) -> Result<Vec<(String, u64, u64)>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut regions = Vec::new();
    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let mut fields = line.split_whitespace();
        let Some(chrom) = fields.next() else {
            continue;
        };
        if chrom.starts_with('#') || chrom == "track" || chrom == "browser" {
            continue;
        }
        let (start, end) = (|| -> Result<(u64, u64)> {
            let start: u64 = fields.next().context("missing start")?.parse()?;
            let end: u64 = fields.next().context("missing end")?.parse()?;
            if start >= end {
                bail!("start must be less than end");
            }
            Ok((start, end))
        })()
        .with_context(|| format!("{}:{}: {line}", path.display(), i + 1))?;
        regions.push((chrom.to_string(), start, end));
    }
    Ok(regions)
}

pub fn run_evidence(cli: &EvidenceCli) -> Result<()> {
    let started = start_wall_timer();
    greetings2(&cli);
    if !cli.validate() {
        bail!("invalid arguments, see errors above");
    }

    // Parse before creating the output so a bad region file leaves no empty result behind.
    let regions = match (&cli.region, &cli.region_bed) {
        (Some(region), _) => {
            let (chrom, start, end) = parse_region(region)?;
            Some(vec![(chrom.to_string(), start, end)])
        }
        (None, Some(bed)) => Some(parse_bed(bed)?),
        (None, None) => None,
    };

    let sample_names =
        DatasetInfo::load_from_file(cli.idxdir.join(DATASET_INFO_FILE_NAME))?.get_sample_names();
    let output_path = |suffix: &str| {
        let mut path = cli.output.as_os_str().to_os_string();
        path.push(suffix);
        PathBuf::from(path)
    };
    let gz_writer = |path: &PathBuf| -> Result<GzEncoder<BufWriter<File>>> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(GzEncoder::new(BufWriter::new(file), Compression::default()))
    };

    let src_path = cli.src.then(|| output_path(".evidence.src.tsv.gz"));
    let mut writer = gz_writer(&output_path(".evidence.gtf.gz"))?;
    write_header(&mut writer, &sample_names, src_path.as_ref())?;
    let mut src_writer = match &src_path {
        Some(path) => {
            let mut w = gz_writer(path)?;
            writeln!(
                w,
                "transcript_id\tsample_id\tsample_name\tstart\tend\tstrand"
            )?;
            Some(w)
        }
        None => None,
    };

    // Filter, then write the GTF records and optionally the per-read source rows.
    let mut emit = |ptir: &PTIR| -> Result<()> {
        if let Some(n_samples) = filter_ptir(ptir, cli.min_read, cli.min_sample) {
            ptir_to_gtf_records(ptir, n_samples, &mut writer)?;
            if let Some(src) = src_writer.as_mut() {
                write_src(ptir, &sample_names, src)?;
            }
        }
        Ok(())
    };

    match regions {
        Some(regions) => write_regions(cli, regions, &mut emit)?,
        None => {
            for ptir in SequentialArchiveReader::load_from_disk(&cli.idxdir.join(MERGED_FILE_NAME))?
            {
                emit(&ptir?)?;
            }
        }
    }

    writer.finish()?.flush()?;
    if let Some(src) = src_writer {
        src.finish()?.flush()?;
    }
    log_wall_time("call", started);

    Ok(())
}

/// Stream the regions' index hits from the B+ tree, loading each PTIR once.
/// Regions are grouped by chromosome (one tree load each) and merged, so keys stay ascending per chromosome.
fn write_regions(
    cli: &EvidenceCli,
    regions: Vec<(String, u64, u64)>,
    emit: &mut dyn FnMut(&PTIR) -> Result<()>,
) -> Result<()> {
    let chroms = ChromMapping::decode(&std::fs::read(cli.idxdir.join(CHROM_FILE_NAME))?);
    let mut by_chrom: BTreeMap<u16, Vec<(u64, u64)>> = BTreeMap::new();
    let mut missing = FxHashSet::default();
    for (chrom, start, end) in regions {
        match chroms.get_chrom_idx(&chrom) {
            Some(id) => by_chrom.entry(id).or_default().push((start, end)),
            None if missing.insert(chrom.clone()) => {
                warn!("chromosome {chrom} is not in the index; skipping its regions")
            }
            None => {}
        }
    }

    // Same cache defaults as `isopedia splice`: 1000 tree nodes, 4 x 128 MB archive chunks.
    let mut archive = PTIRArchiveCache::new(cli.idxdir.join(MERGED_FILE_NAME), 128 << 20, 4);
    for (chrom_id, ranges) in by_chrom {
        let mut tree = BPTree::from_disk(&cli.idxdir, chrom_id, 1000);
        let mut window = SeenWindow::default();
        for (start, end) in merge_ranges(ranges) {
            tree.scan_range(start, end, |pos, ptr| {
                if window.contains(pos, ptr.offset) {
                    return Ok(());
                }
                let ptir = archive.load_from_disk(ptr);
                let right = ptir.isoform_reads_slim_vec.iter().map(|r| r.right).max();
                window.insert(ptr.offset, right.unwrap_or(0));
                emit(&ptir)
            })?;
        }
    }
    Ok(())
}

/// Sort and merge key ranges sharing any key. SeenWindow needs ascending keys, and overlapping
/// scans would revisit keys and report PTIRs twice.
fn merge_ranges(mut ranges: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    ranges.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ranges.len());
    for (start, end) in ranges {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

#[test]
fn merge_ranges_joins_ranges_sharing_a_key() {
    let ranges = vec![(20, 30), (0, 5), (5, 9), (10, 12), (25, 28), (11, 15)];
    assert_eq!(merge_ranges(ranges), vec![(0, 9), (10, 15), (20, 30)]);
}

#[test]
fn parse_bed_skips_headers_and_rejects_bad_lines() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("r.bed");
    std::fs::write(&bed, "track name=x\n#c\n\nchr1\t10\t20\tgeneA\n22 0 5\n").unwrap();
    assert_eq!(
        parse_bed(&bed).unwrap(),
        vec![("chr1".to_string(), 10, 20), ("22".to_string(), 0, 5)]
    );
    for bad in ["chr1\t5\t5\n", "chr1\t5\n", "chr1\tx\t9\n"] {
        std::fs::write(&bed, bad).unwrap();
        assert!(parse_bed(&bed).is_err(), "{bad:?}");
    }
}

/// Dedup for keys arriving in ascending order. Every index key of a PTIR (splice sites, read ends)
/// is <= its rightmost read end, so it can be forgotten once the scan passes that position.
/// ponytail: memory grows with the PTIRs overlapping one position (coverage depth), not region size.
#[derive(Default)]
struct SeenWindow {
    offsets: FxHashSet<u64>,
    rights: BinaryHeap<Reverse<(u64, u64)>>, // (right, offset), smallest right first
}

impl SeenWindow {
    /// Forget PTIRs ending before `pos`, then report whether `offset` was already seen.
    fn contains(&mut self, pos: u64, offset: u64) -> bool {
        while let Some(&Reverse((right, old))) = self.rights.peek() {
            if right >= pos {
                break;
            }
            self.rights.pop();
            self.offsets.remove(&old);
        }
        self.offsets.contains(&offset)
    }

    fn insert(&mut self, offset: u64, right: u64) {
        self.offsets.insert(offset);
        self.rights.push(Reverse((right, offset)));
    }
}

#[test]
fn seen_window_forgets_only_after_right_end() {
    let mut window = SeenWindow::default();
    assert!(!window.contains(10, 1));
    window.insert(1, 20);
    assert!(window.contains(15, 1));
    assert!(window.contains(20, 1)); // a key equal to the right end may still arrive
    assert!(!window.contains(21, 1));
    assert!(window.offsets.is_empty() && window.rights.is_empty());
}

#[test]
fn parse_region_maps_to_key_range() {
    assert_eq!(parse_region("chr22").unwrap(), ("chr22", 0, u64::MAX));
    assert_eq!(parse_region("chr22:101-200").unwrap(), ("chr22", 100, 200));
    for bad in ["", "chr22:0-5", "chr22:9-5", "chr22:5", ":1-2", "chr22:a-5"] {
        assert!(parse_region(bad).is_err(), "{bad}");
    }
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

/// Keep a PTIR when at least `min_sample` samples each have at least `min_read` reads;
/// returns that number of positive samples.
fn filter_ptir(ptir: &PTIR, min_read: u32, min_sample: usize) -> Option<u32> {
    let n_samples = ptir.get_positive_count(&min_read);
    (n_samples as usize >= min_sample).then_some(n_samples)
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
    assert_eq!(filter_ptir(&ptir, 1, 3), Some(3));
    assert_eq!(filter_ptir(&ptir, 2, 2), Some(2));
    assert_eq!(filter_ptir(&ptir, 2, 3), None);
    assert_eq!(filter_ptir(&ptir, 3, 1), Some(1));
    assert_eq!(filter_ptir(&ptir, 4, 1), None);
}

/// Read ends jitter by a few bp, so transcript ends count read ends within ±END_WINDOW bp.
const END_WINDOW: u64 = 20;

/// The position with the most values within ±END_WINDOW (`positions` must be non-empty);
/// ties go to the outermost position: the smallest if `prefer_min`, else the largest.
fn window_mode(mut positions: Vec<u64>, prefer_min: bool) -> u64 {
    positions.sort_unstable();
    let (mut lo, mut hi) = (0, 0);
    let (mut best, mut best_count) = (positions[0], 0);
    for &p in &positions {
        while positions[lo] + END_WINDOW < p {
            lo += 1;
        }
        while hi < positions.len() && positions[hi] <= p + END_WINDOW {
            hi += 1;
        }
        let count = hi - lo;
        if count > best_count || (count == best_count && !prefer_min) {
            (best, best_count) = (p, count);
        }
    }
    best
}

#[test]
fn window_mode_picks_densest_end_and_breaks_ties_outward() {
    // One outlier far left: the cluster at 1000..1010 wins for both directions.
    let xs = vec![10, 1000, 1005, 1010, 1100];
    assert_eq!(window_mode(xs.clone(), true), 1000);
    assert_eq!(window_mode(xs, false), 1010);
    // Exactly END_WINDOW apart still counts as one window.
    assert_eq!(window_mode(vec![500, 520, 900], true), 500);
    assert_eq!(window_mode(vec![500, 521, 900], true), 500); // three singletons tie -> smallest
    assert_eq!(window_mode(vec![500, 521, 900], false), 900); // tie -> largest
    assert_eq!(window_mode(vec![7, 7, 7], true), 7);
}

/// Write one row per supporting read: transcript_id, sample_id, sample_name, start, end, strand (1-based, inclusive).
fn write_src(ptir: &PTIR, sample_names: &[String], writer: &mut dyn Write) -> Result<()> {
    let id = format!("{:016x}", ptir.signature);
    let per_sample = ptir.sample_offset_arr.iter().zip(&ptir.sample_evidence_arr);
    for (sample, (&offset, &count)) in per_sample.enumerate() {
        if count == 0 {
            continue;
        }
        let name = sample_names
            .get(sample)
            .with_context(|| format!("PTIR isopedia.{id} has reads for unknown sample {sample}"))?;
        let reads = ptir
            .isoform_reads_slim_vec
            .get(offset as usize..(offset + count) as usize)
            .with_context(|| {
                format!("PTIR isopedia.{id}: reads of sample {sample} out of range")
            })?;
        for read in reads {
            writeln!(
                writer,
                "isopedia.{id}\tS{}\t{name}\t{}\t{}\t{}",
                sample + 1,
                read.left + 1,
                read.right,
                read.strand
            )?;
        }
    }
    Ok(())
}

/// GTF header in the style of isomatch `##ISOM` lines.
fn write_header(
    writer: &mut dyn Write,
    sample_names: &[String],
    src: Option<&PathBuf>,
) -> Result<()> {
    let escape = |value: &str| value.replace('\\', "\\\\").replace('"', "\\\"");
    writeln!(
        writer,
        "##ISOP <VERSION> version=\"{}\"; program=\"isopedia\"; schema=\"1\"",
        env!("CARGO_PKG_VERSION")
    )?;
    for (i, name) in sample_names.iter().enumerate() {
        writeln!(
            writer,
            "##ISOP <SAMPLE> id=\"S{}\"; name=\"{}\";",
            i + 1,
            escape(name)
        )?;
    }
    for (id, description) in [
        ("ISOP_READ_CNT", "supporting reads of this transcript as total,plus,minus".to_string()),
        ("ISOP_SAMPLE_CNT", "number of samples with at least --min-read supporting reads".to_string()),
        ("ISOP_SPAN", format!("start:end covered by all supporting reads (1-based, inclusive); transcript start/end are the read ends with the most read ends within +/-{END_WINDOW} bp")),
    ] {
        writeln!(writer, "##ISOP <FORMAT> ID=\"{id}\"; Description=\"{description}\";")?;
    }
    if let Some(name) = src.and_then(|path| path.file_name()) {
        writeln!(
            writer,
            "##ISOP <SRC> file=\"{}\";",
            escape(&name.to_string_lossy())
        )?;
    }
    let command = std::env::args_os()
        .map(|arg| escape(&arg.to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(writer, "##ISOP <COMMAND> cmd={command}")?;
    Ok(())
}

#[test]
fn write_src_maps_reads_to_samples() {
    let read = |left, strand| crate::ptir::ReadDiffSlim {
        left,
        right: left + 100,
        strand,
        supp_seg_vec_offset: 0,
        supp_seg_vec_length: 0,
        info: vec![],
    };
    use crate::strand::Strand::{Minus, Plus};
    let ptir = PTIR {
        signature: 9,
        total_evidence: 3,
        sample_size: 3,
        sample_evidence_arr: vec![1, 0, 2],
        sample_offset_arr: vec![0, 1, 1],
        chrom: "1".into(),
        chrom_id: 0,
        rec_type: crate::ptir::RecordType::MonoExon,
        splice_junctions_vec: vec![],
        isoform_reads_slim_vec: vec![read(10, Plus), read(20, Minus), read(30, Plus)],
        supp_segs_vec: vec![],
    };
    let names = ["a", "b", "c"].map(String::from);
    let mut output = Vec::new();
    write_src(&ptir, &names, &mut output).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "isopedia.0000000000000009\tS1\ta\t11\t110\t+\nisopedia.0000000000000009\tS3\tc\t21\t120\t-\nisopedia.0000000000000009\tS3\tc\t31\t130\t+\n"
    );
    assert!(write_src(&ptir, &names[..2], &mut Vec::new()).is_err());
}

/// `n_samples` is the positive sample count written as ISOP_SAMPLE_CNT on the transcript line.
fn ptir_to_gtf_records(ptir: &PTIR, n_samples: u32, writer: &mut dyn Write) -> Result<()> {
    let first = ptir
        .isoform_reads_slim_vec
        .first()
        .with_context(|| format!("PTIR on {} has no reads", ptir.chrom))?;
    let (mut left, mut right) = (first.left, first.right);
    let (mut plus, mut minus) = (0u32, 0u32);
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
        match read.strand {
            crate::strand::Strand::Plus => plus += 1,
            crate::strand::Strand::Minus => minus += 1,
        }
    }
    // Majority strand of the supporting reads; a tie is unknown.
    let strand = match plus.cmp(&minus) {
        Ordering::Greater => "+",
        Ordering::Less => "-",
        Ordering::Equal => ".",
    };
    // Transcript ends are the densest read ends; the full read envelope [left, right) is kept as ISOP_SPAN.
    // Every read carries the whole junction chain, so these ends always enclose it.
    let reads = &ptir.isoform_reads_slim_vec;
    let tx_start = window_mode(reads.iter().map(|r| r.left).collect(), true);
    let tx_end = window_mode(reads.iter().map(|r| r.right).collect(), false);

    let mut exons = Vec::new();
    let mut start = tx_start;
    if !ptir.is_mono_exonic() {
        for &(donor, acceptor) in &ptir.splice_junctions_vec {
            if start > donor || donor >= acceptor || acceptor > tx_end {
                bail!(
                    "invalid PTIR splice junction on {}: {donor}-{acceptor}",
                    ptir.chrom
                );
            }
            if start == donor {
                log::warn!(
                    "PTIR on {}:{}-{}: skipping zero-length exon at {}; continuing GTF output",
                    ptir.chrom,
                    tx_start,
                    tx_end,
                    start
                );
            } else {
                exons.push((start, donor));
            }
            start = acceptor;
        }
    }
    if start == tx_end {
        log::warn!(
            "PTIR on {}:{}-{}: skipping zero-length terminal exon at {}; continuing GTF output",
            ptir.chrom,
            tx_start,
            tx_end,
            start
        );
    } else {
        exons.push((start, tx_end));
    }

    // PTIR coordinates are 0-based, half-open; GTF coordinates are 1-based, inclusive.
    // ponytail: use the PTIR signature as a placeholder gene ID until gene grouping exists.
    let id = format!("{:016x}", ptir.signature);
    writeln!(writer,
        "{}\tisopedia\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"isopedia.{id}\"; transcript_id \"isopedia.{id}\"; ISOP_READ_CNT \"{},{plus},{minus}\"; ISOP_SAMPLE_CNT \"{n_samples}\"; ISOP_SPAN \"{}:{right}\";",
        ptir.chrom, tx_start + 1, tx_end, strand, plus + minus, left + 1
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
        ptir_to_gtf_records(&ptir, 1, &mut output).unwrap();
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
        ptir_to_gtf_records(&ptir, 1, &mut output).unwrap();
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
    assert!(ptir_to_gtf_records(&ptir, 1, &mut output).is_err());
    assert!(output.is_empty());

    // Strand is the majority of read strands (a tie is "."); counts go on the transcript line only.
    use crate::strand::Strand::{Minus, Plus};
    ptir.splice_junctions_vec = vec![(9554665, 9557217)];
    let read = ptir.isoform_reads_slim_vec[0].clone();
    for (strands, expected, counts) in [
        (&[Plus, Minus, Plus][..], "+", "3,2,1"),
        (&[Minus, Plus, Minus][..], "-", "3,1,2"),
        (&[Plus, Minus][..], ".", "2,1,1"),
    ] {
        ptir.isoform_reads_slim_vec = strands
            .iter()
            .map(|&strand| crate::ptir::ReadDiffSlim {
                strand,
                ..read.clone()
            })
            .collect();
        let mut output = Vec::new();
        ptir_to_gtf_records(&ptir, 7, &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.lines().all(|l| l.split('\t').nth(6) == Some(expected)));
        let attrs = format!(
            "ISOP_READ_CNT \"{counts}\"; ISOP_SAMPLE_CNT \"7\"; ISOP_SPAN \"9554546:9557902\";"
        );
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].ends_with(&attrs));
        assert!(lines[1..].iter().all(|l| !l.contains("ISOP_")));
    }
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
