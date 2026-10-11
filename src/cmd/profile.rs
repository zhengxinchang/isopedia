use clap::Parser;
use indexmap::IndexSet;
use log::info;
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
    Read, Record,
};
use rustc_hash::FxHashMap;
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};

use crate::error::ProfileError;
use crate::strand::Strand;
use crate::{
    myio::MyGzWriter,
    reads::{profile_header, AggrRead, SingleRead},
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

    /// Verbose mode
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
}

pub struct ProfileStats {
    total: u64,
    skipped_unmapped: u64,
    skipped_secondary: u64,
    skipped_low_mapq: u64,
    mapq: [u64; 256],
    // Summary indices: 0 = all, 1 = mono exon, 2 = multiple exon.
    alignment_count: [u64; 3],
    length_observed: [u64; 3],
    total_bases: [u64; 3],
    length_min: [Option<u64>; 3],
    length_max: [Option<u64>; 3],
    group_count: [u64; 3],
    group_observed: [u64; 3],
    group_mean_sum: [f64; 3],
    group_mean_min: [Option<f64>; 3],
    group_mean_max: [Option<f64>; 3],
}

impl ProfileStats {
    pub fn new() -> Self {
        Self {
            total: 0,
            skipped_unmapped: 0,
            skipped_secondary: 0,
            skipped_low_mapq: 0,
            mapq: [0; 256],
            alignment_count: [0; 3],
            length_observed: [0; 3],
            total_bases: [0; 3],
            length_min: [None; 3],
            length_max: [None; 3],
            group_count: [0; 3],
            group_observed: [0; 3],
            group_mean_sum: [0.0; 3],
            group_mean_min: [None; 3],
            group_mean_max: [None; 3],
        }
    }

    pub fn add_mapq(&mut self, mapq: u8) {
        self.mapq[mapq as usize] += 1;
    }

    fn add_alignment(&mut self, length: u64, mono: bool) {
        for i in [0, if mono { 1 } else { 2 }] {
            self.alignment_count[i] += 1;
            if length > 0 {
                self.length_observed[i] += 1;
                self.total_bases[i] += length;
                self.length_min[i] = Some(self.length_min[i].map_or(length, |v| v.min(length)));
                self.length_max[i] = Some(self.length_max[i].map_or(length, |v| v.max(length)));
            }
        }
    }

    fn add_group(&mut self, observed: u64, total_bases: u64, mono: bool) {
        for i in [0, if mono { 1 } else { 2 }] {
            self.group_count[i] += 1;
            if observed > 0 {
                let mean = total_bases as f64 / observed as f64;
                self.group_observed[i] += 1;
                self.group_mean_sum[i] += mean;
                self.group_mean_min[i] = Some(self.group_mean_min[i].map_or(mean, |v| v.min(mean)));
                self.group_mean_max[i] = Some(self.group_mean_max[i].map_or(mean, |v| v.max(mean)));
            }
        }
    }

    pub fn get_stats(&self) -> Value {
        let ratio = |sum: f64, count: u64| (count > 0).then(|| sum / count as f64);
        let mut alignments = json!({});
        let mut groups = json!({});
        for (i, name) in ["all", "mono_exon", "multiple_exon"].iter().enumerate() {
            alignments[*name] = json!({
                "count": self.alignment_count[i],
                "length": {
                    "observed": self.length_observed[i],
                    "missing": self.alignment_count[i] - self.length_observed[i],
                    "total_bases": self.total_bases[i],
                    "min": self.length_min[i], "max": self.length_max[i],
                    "mean": ratio(self.total_bases[i] as f64, self.length_observed[i]),
                }
            });
            groups[*name] = json!({
                "count": self.group_count[i],
                "mean_alignment_length": {
                    "observed": self.group_observed[i],
                    "missing": self.group_count[i] - self.group_observed[i],
                    "min": self.group_mean_min[i], "max": self.group_mean_max[i],
                    "mean": ratio(self.group_mean_sum[i], self.group_observed[i]),
                }
            });
            if i != 0 {
                alignments[*name]["fraction"] = json!(ratio(
                    self.alignment_count[i] as f64,
                    self.alignment_count[0]
                ));
                groups[*name]["fraction"] =
                    json!(ratio(self.group_count[i] as f64, self.group_count[0]));
            }
        }
        let counted: u64 = self.mapq.iter().sum();
        let bins: Vec<_> = self
            .mapq
            .iter()
            .enumerate()
            .map(|(i, &count)| {
                json!({
                    "label": i.to_string(),
                    "count": count, "fraction": ratio(count as f64, counted),
                })
            })
            .collect();
        let compact_bins: Vec<_> = self
            .mapq
            .chunks(20)
            .enumerate()
            .map(|(i, bucket)| {
                let count: u64 = bucket.iter().sum();
                let start = i * 20;
                json!({
                    "label": format!("{}-{}", start, start + bucket.len() - 1),
                    "count": count, "fraction": ratio(count as f64, counted),
                })
            })
            .collect();
        json!({
            "program_version": env!("CARGO_PKG_VERSION"),
            "alignment_filter_summary": {
                "total": self.total,
                "retained": self.alignment_count[0],
                "skipped": self.skipped_unmapped + self.skipped_secondary + self.skipped_low_mapq,
                "skipped_by_reason": {
                    "unmapped": self.skipped_unmapped,
                    "secondary": self.skipped_secondary,
                    "low_mapq": self.skipped_low_mapq,
                }
            },
            "mapq_compact": { "counted": counted, "bins": compact_bins },
            "mapq": { "counted": counted, "bins": bins },
            "profile_summary": {
                "alignment_records": alignments,
                "aggregated_isoform_records": groups,
            }
        })
    }

    fn format_stats(&self) -> serde_json::Result<Vec<u8>> {
        let stats = self.get_stats();
        let mut output = String::from("{\n");
        for (i, key) in [
            "program_version",
            "alignment_filter_summary",
            "profile_summary",
            "mapq_compact",
            "mapq",
        ]
        .iter()
        .enumerate()
        {
            if i > 0 {
                output.push_str(",\n");
            }
            output.push_str("  ");
            output.push_str(&serde_json::to_string(key)?);
            output.push_str(": ");
            if *key == "mapq" {
                output.push_str(&serde_json::to_string(&stats[*key])?);
            } else {
                let value = serde_json::to_string_pretty(&stats[*key])?;
                output.push_str(&value.replace('\n', "\n  "));
            }
        }
        output.push_str("\n}\n");
        Ok(output.into_bytes())
    }
}

fn write_profile_json(
    path: &Path,
    contents: serde_json::Result<Vec<u8>>,
) -> Result<(), ProfileError> {
    let temporary = path.with_extension("json.tmp");
    let result = (|| -> std::io::Result<()> {
        let bytes = contents.map_err(std::io::Error::other)?;
        fs::write(&temporary, bytes)?;
        fs::rename(&temporary, path)
    })();
    result.map_err(|source| ProfileError::Io {
        context: format!("write profile JSON {}", path.display()),
        source,
    })
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

    cli.validate()?;
    let mut info_name = cli.output.as_os_str().to_owned();
    info_name.push(".profile_info.json");
    let info_path = PathBuf::from(info_name);
    let mut params_name = cli.output.as_os_str().to_owned();
    params_name.push(".profile_params.json");
    let params_path = PathBuf::from(params_name);
    match fs::remove_file(&info_path) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(source) => {
            return Err(ProfileError::Io {
                context: format!("remove previous statistics {}", info_path.display()),
                source,
            })
        }
    }
    write_profile_json(&params_path, serde_json::to_vec_pretty(cli))?;
    greetings2(&cli);

    let mut chrom_set: IndexSet<String> = IndexSet::new();
    let mut agg_isoform_map: FxHashMap<u128, (AggrRead, u64, u64, bool)> = FxHashMap::default();
    let mut profile_stats = ProfileStats::new();

    let batch_size = 1_000_000;

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
            result.map_err(|source| ProfileError::Hts {
                context: format!(
                    "read record {} from {}",
                    profile_stats.total + 1,
                    bam_path.display()
                ),
                source,
            })?;

            profile_stats.total += 1;

            if profile_stats.total % batch_size == 0 {
                info!(
                    "Processed {} records",
                    profile_stats.total.to_formatted_string(&Locale::en)
                );
            }

            if record.is_unmapped() {
                profile_stats.skipped_unmapped += 1;
                continue;
            }

            if record.is_secondary() {
                if !cli.use_secondary {
                    profile_stats.skipped_secondary += 1;
                    continue;
                }
            }

            profile_stats.add_mapq(record.mapq());

            if record.mapq() < cli.mapq {
                profile_stats.skipped_low_mapq += 1;
                continue;
            }

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

            let mono = single_read.pri_seg_size == 1;
            let length = single_read.read_len;
            let observed = u64::from(length > 0);
            profile_stats.add_alignment(length, mono);
            match agg_isoform_map.entry(single_read.signature) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    let (agg_isoform, length_count, total_bases, _) = entry.get_mut();
                    *length_count += observed;
                    *total_bases += length;
                    agg_isoform.add_read(single_read);
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert((AggrRead::new(single_read), observed, length, mono));
                }
            }
        }
    }

    // write the chrom_set
    info!("Processing chromosomes");
    let chrom_str = profile_header(chrom_set.iter().map(String::as_str));
    info!("Write to output file: {}", cli.output.display());
    mywriter
        .write_all_bytes(chrom_str.as_bytes())
        .map_err(|source| ProfileError::Io {
            context: format!("write header to {}", cli.output.display()),
            source,
        })?;

    if cli.verbose {
        // aggr isoform size
        info!("Number of aggr isoforms: {}", agg_isoform_map.len());
    }

    let mut sorted_agg_isoform_map = agg_isoform_map.into_iter().collect::<Vec<_>>();
    sorted_agg_isoform_map.sort_by_key(|x| x.0);
    for (_, (agg_isoform, observed, total_bases, mono)) in &sorted_agg_isoform_map {
        profile_stats.add_group(*observed, *total_bases, *mono);
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
        profile_stats.total.to_formatted_string(&Locale::en)
    );

    info!(
        "Total valid records: {}",
        profile_stats.alignment_count[0].to_formatted_string(&Locale::en)
    );
    let skipped_records = profile_stats.skipped_unmapped
        + profile_stats.skipped_secondary
        + profile_stats.skipped_low_mapq;
    info!(
        "Filtered records summary: unmapped={}, low_mapq(< {})={}, secondary(disabled by --use-secondary)={}",
        profile_stats.skipped_unmapped.to_formatted_string(&Locale::en),
        cli.mapq,
        profile_stats.skipped_low_mapq.to_formatted_string(&Locale::en),
        profile_stats.skipped_secondary.to_formatted_string(&Locale::en)
    );
    info!(
        "Total records filtered out: {}",
        skipped_records.to_formatted_string(&Locale::en)
    );

    write_profile_json(&info_path, profile_stats.format_stats())?;

    info!("Finished");
    log_wall_time("profile", started);
    Ok(())
}

#[cfg(test)]
#[path = "../tests/cmd_profile.rs"]
mod tests;
