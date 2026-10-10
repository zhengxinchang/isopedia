use clap::Parser;
use indexmap::IndexSet;
use log::info;
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
    Read, Record,
};
use rustc_hash::FxHashMap;
use std::fs;
use std::path::PathBuf;

use crate::error::ProfileError;
use crate::strand::Strand;
use crate::{
    myio::MyGzWriter,
    reads::{AggrRead, SingleRead},
    utils::{greetings2, log_wall_time, start_wall_timer},
};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia profile")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[build index, step1] Profile raw isoform signals from each single BAM/CRAM file.
", long_about = None)]
#[clap(after_long_help = "

Example: isopedia profile -i /path/to/bam/file.bam -o /path/to/output/isoform.gz

Note that if you are using CRAM file, you must provide the reference fasta file with --reference option.

")]
pub struct ProfileCli {
    /// Input file in BAM/CRAM format
    #[arg(short = 'i', long = "bam")]
    pub bam: PathBuf,

    /// Reference file for CRAMs. Must provide for CRAM format input
    #[arg(short, long)]
    pub reference: Option<PathBuf>,

    /// Name of the output signal file
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimal mapping quality of reads
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Include secondary alignments
    #[arg(long = "use-secondary", default_value_t = false)]
    pub use_secondary: bool,

    /// Include read names in the output
    #[arg(long = "rname", default_value_t = false)]
    pub rname: bool,

    /// Output alignment stats with the name of output.stats.txt
    #[arg(long = "stats", default_value_t = false)]
    pub stats: bool,

    /// Verbose mode
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
}

pub struct ProfileStats {
    mapq: [u32; 62],
}

impl ProfileStats {
    pub fn new() -> Self {
        ProfileStats { mapq: [0u32; 62] }
    }

    pub fn add_mapq(&mut self, mapq: u8) {
        let mapq_idx = usize::min(mapq as usize, 61);
        self.mapq[mapq_idx] += 1;
    }

    pub fn get_stats(&self) -> String {
        let mut report = "MAPQ Hist:\n".to_string();
        let total_count: u32 = self.mapq.iter().sum();

        report.push_str(&format!(
            "Total alignments counted: {}\n",
            total_count.to_formatted_string(&Locale::en)
        ));

        if total_count == 0 {
            report.push_str("No alignments available for MAPQ statistics.\n");
            return report;
        }

        let max_count = *self.mapq.iter().max().unwrap_or(&0);
        let max_bar_width = 50usize;

        for (idx, count) in self.mapq.iter().enumerate() {
            if *count == 0 {
                continue;
            }

            let pct = (*count as f64) * 100.0 / (total_count as f64);
            let bar_len = if max_count > 0 {
                ((*count as f64) / (max_count as f64) * (max_bar_width as f64)).round() as usize
            } else {
                0
            };
            let bar = "*".repeat(bar_len.max(1));
            let mapq_label = if idx == 61 {
                ">60".to_string()
            } else {
                idx.to_string()
            };

            report.push_str(&format!(
                "{:>4}: {:>12} ({:>6.2}%) {}\n",
                mapq_label,
                count.to_formatted_string(&Locale::en),
                pct,
                bar
            ));
        }

        report
    }
}

impl ProfileCli {
    fn is_cram_input(&self) -> bool {
        self.bam
            .extension()
            .and_then(|ext| ext.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("cram"))
            .unwrap_or(false)
    }

    pub fn validate(&self) -> Result<(), ProfileError> {
        if !self.bam.exists() {
            return Err(ProfileError::InvalidArgument(format!(
                "--bam: input file {} does not exist",
                self.bam.display()
            )));
        }
        if self.is_cram_input() && self.reference.is_none() {
            return Err(ProfileError::InvalidArgument(format!(
                "--reference must be provided for CRAM input: {}",
                self.bam.display()
            )));
        }
        if self.output.is_dir() || self.output.file_name().is_none() {
            return Err(ProfileError::InvalidArgument(format!(
                "--output: {} must name a file",
                self.output.display()
            )));
        }
        if self.rname {
            info!("--rname option is set, read names will be included in the output.");
        }
        Ok(())
    }
}

pub fn run_profile(cli: &ProfileCli) -> Result<(), ProfileError> {
    let started = start_wall_timer();

    greetings2(&cli);
    cli.validate()?;

    let mut chrom_set: IndexSet<String> = IndexSet::new();
    let mut agg_isoform_map: FxHashMap<u64, AggrRead> = FxHashMap::default();
    let mut profile_stats = ProfileStats::new();

    let mut n_record = 0u32;
    // let mut batches = 0;
    let mut total_count = 0;
    let batch_size = 1000000; // 1 million records per batch
    let mut skipped_unmapped = 0;
    let mut skipped_secondary = 0;
    let mut skipped_low_mapq = 0;

    let mut mywriter = MyGzWriter::new(&cli.output).map_err(|source| ProfileError::Io {
        context: format!("create output {}", cli.output.display()),
        source,
    })?;

    {
        info!("Using BAM/CRAM file for isoform profiling.");
        let bam_path = &cli.bam;
        let mut bam_reader =
            bam::Reader::from_path(bam_path).map_err(|source| ProfileError::Hts {
                context: format!("open input {}", bam_path.display()),
                source,
            })?;

        if cli.is_cram_input() {
            let reference = cli.reference.as_ref().ok_or_else(|| {
                ProfileError::InvalidArgument("--reference is required for CRAM input".into())
            })?;
            bam_reader
                .set_reference(reference)
                .map_err(|source| ProfileError::Hts {
                    context: format!(
                        "set reference {} for {}",
                        reference.display(),
                        bam_path.display()
                    ),
                    source,
                })?;
        }

        let mut record = Record::new();
        let header = bam_reader.header().to_owned();

        while let Some(result) = bam_reader.read(&mut record) {
            total_count += 1;

            result.map_err(|source| ProfileError::Hts {
                context: format!("read record {} from {}", total_count, bam_path.display()),
                source,
            })?;

            if total_count % batch_size == 0 {
                info!(
                    "Processed {} records",
                    total_count.to_formatted_string(&Locale::en)
                );
                // batches += 1;
            }

            if record.is_unmapped() {
                skipped_unmapped += 1;
                continue;
            }

            if record.is_secondary() {
                if !cli.use_secondary {
                    skipped_secondary += 1;
                    continue;
                }
            }

            if cli.stats {
                profile_stats.add_mapq(record.mapq());
            }

            if record.mapq() < cli.mapq {
                skipped_low_mapq += 1;
                continue;
            }

            n_record += 1;

            let invalid_record = |reason: String| ProfileError::InvalidRecord {
                read_name: String::from_utf8_lossy(record.qname()).into_owned(),
                reason,
            };
            if record.tid() < 0 || record.tid() as u32 >= header.target_count() {
                return Err(invalid_record(format!(
                    "invalid reference ID {}",
                    record.tid()
                )));
            }
            let chrom = String::from_utf8(header.tid2name(record.tid() as u32).to_vec())
                .map_err(|e| invalid_record(format!("invalid reference name: {e}")))?;
            let chrom = crate::utils::trim_chr_prefix_to_upper(&chrom);
            chrom_set.insert(chrom.clone());

            let strand = record.strand().to_owned(); // MUST move out of match, Because of mutable borrow by strand().
            let my_strand = Strand::from_bam_strand(strand);
            let mut left = u64::try_from(record.pos())
                .map_err(|e| invalid_record(format!("invalid alignment position: {e}")))?;
            let mut right = left;
            let mut single_read = SingleRead::new(
                chrom.clone(),
                record.mapq(),
                1,
                record.seq().len() as u64,
                left,
                my_strand,
            );
            // process the main CIGAR

            for cigar in record.cigar().iter() {
                // Match(u32),    // M
                // Ins(u32),      // I
                // Del(u32),      // D
                // RefSkip(u32),  // N
                // SoftClip(u32), // S
                // HardClip(u32), // H
                // Pad(u32),      // P
                // Equal(u32),    // =
                // Diff(u32),     // X
                match cigar {
                    Cigar::Match(n) | Cigar::Del(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                        right = right
                            .checked_add(u64::from(*n))
                            .ok_or_else(|| invalid_record("CIGAR coordinate overflow".into()))?;
                    }
                    Cigar::RefSkip(n) => {
                        single_read.add_segment(chrom.clone(), left, right, &my_strand, false);
                        left = right
                            .checked_add(u64::from(*n))
                            .ok_or_else(|| invalid_record("CIGAR coordinate overflow".into()))?;
                        right = left;
                    }
                    _ => {}
                }
            }
            // handle the last segment in main CIGAR
            // record.tname()
            single_read.add_segment(chrom.clone(), left, right, &my_strand, false);
            // update ref span
            single_read.update_right(right);
            // process the supplementary alignment
            match record.aux(b"SA") {
                Err(rust_htslib::errors::Error::BamAuxTagNotFound) => {}
                Err(source) => {
                    return Err(ProfileError::Hts {
                        context: format!(
                            "read SA tag for {:?} in {}",
                            String::from_utf8_lossy(record.qname()),
                            bam_path.display()
                        ),
                        source,
                    })
                }
                Ok(Aux::String(sa)) => {
                    for entry in sa.split(';').filter(|entry| !entry.is_empty()) {
                        let fields: Vec<_> = entry.split(',').collect();
                        if fields.len() != 6
                            || fields[0].is_empty()
                            || !matches!(fields[2], "+" | "-")
                        {
                            return Err(invalid_record(format!("invalid SA entry: {entry}")));
                        }
                        let mut left = fields[1]
                            .parse::<u64>()
                            .map_err(|e| invalid_record(format!("invalid SA position: {e}")))?
                            .checked_sub(1)
                            .ok_or_else(|| invalid_record("SA position must be positive".into()))?;
                        fields[4]
                            .parse::<u8>()
                            .map_err(|e| invalid_record(format!("invalid SA MAPQ: {e}")))?;
                        fields[5]
                            .parse::<u64>()
                            .map_err(|e| invalid_record(format!("invalid SA NM: {e}")))?;
                        // rust-htslib 0.48 panics on a trailing CIGAR length; validate before parsing.
                        if fields[3].is_empty()
                            || fields[3].as_bytes().last().is_some_and(u8::is_ascii_digit)
                        {
                            return Err(invalid_record(format!("invalid SA CIGAR: {}", fields[3])));
                        }
                        let cigar = CigarString::try_from(fields[3])
                            .map_err(|e| invalid_record(format!("invalid SA CIGAR: {e}")))?;
                        let mut right = left;
                        for op in &cigar {
                            if op.len() == 0 {
                                return Err(invalid_record(
                                    "SA CIGAR operation length must be positive".into(),
                                ));
                            }
                            match op {
                                Cigar::Match(n)
                                | Cigar::Del(n)
                                | Cigar::Equal(n)
                                | Cigar::Diff(n) => {
                                    right = right.checked_add(u64::from(*n)).ok_or_else(|| {
                                        invalid_record("SA coordinate overflow".into())
                                    })?;
                                }
                                Cigar::RefSkip(n) => {
                                    single_read.add_supp_segment(
                                        fields[0].into(),
                                        left,
                                        right,
                                        fields[2],
                                        true,
                                    );
                                    left = right.checked_add(u64::from(*n)).ok_or_else(|| {
                                        invalid_record("SA coordinate overflow".into())
                                    })?;
                                    right = left;
                                }
                                _ => {}
                            }
                        }
                        single_read.add_supp_segment(
                            fields[0].into(),
                            left,
                            right,
                            fields[2],
                            true,
                        );
                    }
                }
                Ok(_) => return Err(invalid_record("SA tag must be a string".into())),
            }

            // add read name to single read
            if cli.rname {
                single_read.add_info(
                    "RN".to_string(),
                    String::from_utf8(record.qname().to_vec())
                        .map_err(|e| invalid_record(format!("invalid read name: {e}")))?,
                );
            }

            single_read.process();

            match agg_isoform_map.get_mut(&single_read.signature) {
                Some(agg_isoform) => {
                    agg_isoform.add_read(single_read);
                }
                None => {
                    let agg_isoform = AggrRead::new(single_read);
                    agg_isoform_map.insert(agg_isoform.signature, agg_isoform);
                }
            }
        }
    }

    // write the chrom_set
    info!("Processing chromosomes");
    let mut chrom_str = Vec::new();
    chrom_set.into_iter().for_each(|x| {
        chrom_str.extend_from_slice(x.as_bytes());
        chrom_str.push(b'\t');
    });
    chrom_str.push(b'\n');

    info!("Write to output file: {}", cli.output.display());

    chrom_str.extend_from_slice(b"signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n");
    mywriter
        .write_all_bytes(&chrom_str)
        .map_err(|source| ProfileError::Io {
            context: format!("write header to {}", cli.output.display()),
            source,
        })?;

    if cli.verbose {
        // aggr isoform size
        info!("Number of aggr isoforms: {}", agg_isoform_map.len());
    }

    let mut sorted_agg_isoform_map = agg_isoform_map
        .into_iter()
        .collect::<Vec<(u64, AggrRead)>>();
    sorted_agg_isoform_map.sort_by_key(|x| x.0);
    for (_, agg_isoform) in sorted_agg_isoform_map.iter() {
        // writer.write(agg_isoform.to_record().as_bytes()).unwrap();
        mywriter
            .write_all_bytes(agg_isoform.to_record().as_bytes())
            .map_err(|source| ProfileError::Io {
                context: format!("write record to {}", cli.output.display()),
                source,
            })?;
    }
    mywriter.finish().map_err(|source| ProfileError::Io {
        context: format!("finish output {}", cli.output.display()),
        source,
    })?;
    info!(
        "Total records processed: {}",
        total_count.to_formatted_string(&Locale::en)
    );

    info!(
        "Total valid records: {}",
        n_record.to_formatted_string(&Locale::en)
    );
    // info!("Total records : {}", n_record);
    let skipped_records = skipped_unmapped + skipped_secondary + skipped_low_mapq;
    info!(
        "Filtered records summary: unmapped={}, low_mapq(< {})={}, secondary(disabled by --use-secondary)={}",
        skipped_unmapped.to_formatted_string(&Locale::en),
        cli.mapq,
        skipped_low_mapq.to_formatted_string(&Locale::en),
        skipped_secondary.to_formatted_string(&Locale::en)
    );
    info!(
        "Total records filtered out: {}",
        skipped_records.to_formatted_string(&Locale::en)
    );

    if cli.stats {
        let stats_path = cli.output.with_file_name(format!(
            "{}.stats.txt",
            cli.output
                .file_name()
                .ok_or_else(|| ProfileError::InvalidArgument("output file name is missing".into()))?
                .to_string_lossy()
        ));
        let mut stats_report = profile_stats.get_stats();
        stats_report.push('\n');
        stats_report.push_str(&format!(
            "Filter summary:\nunmapped: {}\nlow_mapq(< {}): {}\nsecondary(disabled by --use-secondary): {}\ntotal_filtered: {}\n",
            skipped_unmapped.to_formatted_string(&Locale::en),
            cli.mapq,
            skipped_low_mapq.to_formatted_string(&Locale::en),
            skipped_secondary.to_formatted_string(&Locale::en),
            skipped_records.to_formatted_string(&Locale::en)
        ));
        fs::write(&stats_path, stats_report).map_err(|source| ProfileError::Io {
            context: format!("write statistics {}", stats_path.display()),
            source,
        })?;
        info!("Alignment stats written to {}", stats_path.display());
    }

    info!("Finished");
    log_wall_time("profile", started);
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{
        cmd::profile::{run_profile, ProfileCli},
        error::ProfileError,
    };
    use clap::Parser;
    use flate2::read::GzDecoder;
    use rust_htslib::bam::{self, Read};
    use std::{fs, io::Read as _, path::Path};

    #[test]
    fn profile_errors_and_output() {
        let dir = tempfile::tempdir().unwrap();
        let sam = dir.path().join("input.sam");
        let bam = dir.path().join("input.bam");
        let output = dir.path().join("output.gz");
        let mut cli = ProfileCli::parse_from([
            "profile",
            "-i",
            bam.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--rname",
            "--stats",
        ]);
        assert!(matches!(
            run_profile(&cli),
            Err(ProfileError::InvalidArgument(_))
        ));

        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
        let record = "read1\t0\tchr1\t1\t60\t10M20N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t*";
        let write_bam = |tag: &str| {
            fs::write(&sam, format!("{header}{record}{tag}\n")).unwrap();
            let mut reader = bam::Reader::from_path(&sam).unwrap();
            let header = bam::Header::from_template(reader.header());
            let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
            for record in reader.records() {
                writer.write(&record.unwrap()).unwrap();
            }
        };
        write_bam("");
        run_profile(&cli).unwrap();
        let mut text = String::new();
        GzDecoder::new(fs::File::open(&output).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert!(text.contains("\t1\t1\t10-30\t0-40:+::RN=read1\n"));
        assert!(dir.path().join("output.gz.stats.txt").is_file());

        write_bam("\tSA:Z:chr1,101,+,5M10N5M,60,0;");
        run_profile(&cli).unwrap();
        text.clear();
        GzDecoder::new(fs::File::open(&output).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert!(text.contains("1%100%105%+,1%115%120%+"));

        for tag in [
            "\tSA:i:1",
            "\tSA:Z:chr1,101,+;",
            "\tSA:Z:chr1,0,+,10M,60,0;",
            "\tSA:Z:chr1,101,?,10M,60,0;",
            "\tSA:Z:chr1,101,+,10M5,60,0;",
            "\tSA:Z:chr1,101,+,0M,60,0;",
            "\tSA:Z:chr1,101,+,10Q,60,0;",
            "\tSA:Z:chr1,101,+,10M,bad,0;",
            "\tSA:Z:chr1,101,+,10M,60,bad;",
            "\tSA:Z:chr1,18446744073709551615,+,10M,60,0;",
        ] {
            write_bam(tag);
            assert!(
                matches!(run_profile(&cli), Err(ProfileError::InvalidRecord { read_name, .. }) if read_name == "read1"),
                "{tag}"
            );
        }

        write_bam("");
        let bytes = fs::read(&bam).unwrap();
        fs::write(&bam, &bytes[..bytes.len() - 40]).unwrap();
        assert!(matches!(run_profile(&cli), Err(ProfileError::Hts { .. })));
        fs::write(&bam, b"not a BAM").unwrap();
        assert!(matches!(run_profile(&cli), Err(ProfileError::Hts { .. })));

        write_bam("");
        cli.output = dir.path().join("missing/output.gz");
        assert!(matches!(run_profile(&cli), Err(ProfileError::Io { .. })));
        cli.output = dir.path().to_owned();
        assert!(matches!(
            run_profile(&cli),
            Err(ProfileError::InvalidArgument(_))
        ));
        cli.output = output.clone();
        fs::remove_file(dir.path().join("output.gz.stats.txt")).unwrap();
        fs::create_dir(dir.path().join("output.gz.stats.txt")).unwrap();
        assert!(
            matches!(run_profile(&cli), Err(ProfileError::Io { context, .. }) if context.contains("statistics"))
        );

        if Path::new("/dev/full").exists() {
            cli.output = "/dev/full".into();
            cli.stats = false;
            assert!(matches!(run_profile(&cli), Err(ProfileError::Io { .. })));
            write_bam("\tSA:i:1");
            assert!(matches!(
                run_profile(&cli),
                Err(ProfileError::InvalidRecord { .. })
            ));
        }

        cli.bam = dir.path().join("input.cram");
        fs::write(&cli.bam, b"").unwrap();
        assert!(matches!(
            run_profile(&cli),
            Err(ProfileError::InvalidArgument(_))
        ));
    }
}
