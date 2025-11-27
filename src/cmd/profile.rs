use clap::Parser;
use core::panic;
use indexmap::IndexSet;
use log::{error, info, warn};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar},
    Read, Record,
};
use rustc_hash::FxHashMap;
use std::path::PathBuf;

use crate::{
    gtf::{open_gtf_reader, TranscriptChunker},
    io::MyGzWriter,
    reads::{AggrRead, SingleRead},
};
use anyhow::Result;
use num_format::{Locale, ToFormattedString};
use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

impl Strand {
    pub fn from_string(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    #[allow(unused)]
    pub fn from_num_string(s: &str) -> Self {
        match s {
            "0" => Strand::Plus,
            "1" => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    #[allow(unused)]
    pub fn from_u8(n: u8) -> Self {
        match n {
            0 => Strand::Plus,
            1 => Strand::Minus,
            other => panic!("Invalid strand {:?}", other),
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            Strand::Plus => "+".to_string(),
            Strand::Minus => "-".to_string(),
        }
    }
}

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia profile")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[build index, step1] Profile raw isoform signals from each single BAM/CRAM or GTF file.
", long_about = None)]
#[clap(after_long_help = "

Example: isopedia profile -i /path/to/bam/file.bam -o /path/to/output/isoform.gz

Note that if you are using CRAM file, you must provide the reference fasta file with --reference option.

")]
pub struct ProfileCli {
    /// Input file in BAM/CRAM format
    #[arg(short = 'i', long = "bam")]
    pub bam: Option<PathBuf>,

    /// Input file in GTF format
    #[arg(short = 'g', long = "gtf")]
    pub gtf: Option<PathBuf>,

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

    /// Include read names in the output, only for BAM/CRAM input
    #[arg(long = "rname", default_value_t = false)]
    pub rname: bool,

    /// Include transcript IDs in the output, only for GTF input
    #[arg(long = "tid", default_value_t = false)]
    pub tid: bool,

    /// Include gene IDs in the output, only for GTF input
    #[arg(long = "gid", default_value_t = false)]
    pub gid: bool,

    /// Verbose mode
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
}

impl ProfileCli {
    pub fn validate(&self) {
        let mut is_ok = true;

        if self.gtf.is_none() && self.bam.is_none() {
            error!("Either --gtf or --bam must be provided as input.");
            is_ok = false;
        }

        if self.gtf.is_some() && self.bam.is_some() {
            error!("Only one of --gtf or --bam should be provided as input.");
            is_ok = false;
        }

        if self.gtf.is_some() {
            let gtf_path = self.gtf.as_ref().unwrap();
            if !gtf_path.exists() {
                error!("--gtf: input file {} does not exist", gtf_path.display());
                is_ok = false;
            }

            if self.rname {
                warn!("--rname can only be used with BAM/CRAM input.");
            }

            if self.tid {
                info!("--tid option is set, transcript IDs will be included in the output.");
            }
            if self.gid {
                info!("--gid option is set, gene IDs will be included in the output.");
            }
        }

        if self.bam.is_some() {
            let bam_path = self.bam.as_ref().unwrap();
            if !bam_path.exists() {
                error!("--bam: input file {} does not exist", bam_path.display());
                is_ok = false;
            }

            if self.tid {
                warn!("--tid can only be used with GTF input.");
            }

            if self.gid {
                warn!("--gid can only be used with GTF input.");
            }

            if self.rname {
                info!("--rname option is set, read names will be included in the output.");
            }
        }

        // check output is not a directory
        if self.output.is_dir() {
            error!(
                "--output: {} is a directory, not a file",
                self.output.display()
            );
            is_ok = false;
        }

        if is_ok != true {
            panic!("Invalid arguments, please check the error messages above.");
        }
    }
}

fn greetings(args: &ProfileCli) {
    // eprintln!("\nIsopedia: [Profile isoform signals from BAM/CRAM or GTF]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

pub fn run_profile(cli: &ProfileCli) -> Result<()> {
    // env::set_var("RUST_LOG", "info");
    // env_logger::init();

    // let cli = ExtrCli::parse();

    greetings(&cli);
    cli.validate();

    let mut chrom_set: IndexSet<String> = IndexSet::new();
    let mut agg_isoform_map: FxHashMap<u64, AggrRead> = FxHashMap::default();

    let mut n_record = 0u32;
    // let mut batches = 0;
    let mut total_count = 0;
    let batch_size = 1000000; // 1 million records per batch
    let mut skipped_records = 0;

    let mut mywriter = MyGzWriter::new(&cli.output).expect(&format!(
        "Can not create output file {} .",
        &cli.output.display()
    ));

    if cli.gtf.is_none() && cli.bam.is_some() {
        info!("Using BAM/CRAM file for isoform profiling.");
        let bam_path = &cli.bam.as_ref().unwrap();
        let mut bam_reader = bam::Reader::from_path(bam_path).expect("Failed to open BAM file");

        if (&bam_path.ends_with(".cram")) & (&cli.reference.is_some()) {
            bam_reader
                .set_reference(
                    cli.reference
                        .as_ref()
                        .expect("Reference must be provided for CRAM files"),
                )
                .expect("Failed to set reference for CRAM file");
        }

        let mut record = Record::new();
        let header = bam_reader.header().to_owned();

        while let Some(result) = bam_reader.read(&mut record) {
            total_count += 1;

            match result {
                Err(_) => {
                    eprintln!("Can not read record");
                    break;
                } //exit if the last one was processed.
                Ok(()) => {}
            }

            if total_count % batch_size == 0 {
                info!(
                    "Processed {} records",
                    total_count.to_formatted_string(&Locale::en)
                );
                // batches += 1;
            }

            if record.is_unmapped() {
                skipped_records += 1;
                continue;
            }

            if record.is_secondary() {
                if !cli.use_secondary {
                    skipped_records += 1;
                    continue;
                }
            }

            if record.mapq() < cli.mapq {
                skipped_records += 1;
                continue;
            }

            n_record += 1;

            // dont use record.tid() before it is filled by bam.read()
            // otherwise it will be random number and cuase the error with random memory access.
            let mut chrom =
                String::from_utf8(header.tid2name(record.tid() as u32).to_vec()).unwrap();
            if chrom.len() >= 3 && chrom.starts_with("chr") {
                chrom = chrom[3..].to_owned();
            }
            chrom_set.insert(chrom.clone());

            let strand = record.strand().to_owned(); // MUST move out of match, Because of mutable borrow by strand().
            let mut left = record.pos() as u32;
            let mut right = record.pos() as u32;
            let mut isoform = SingleRead::new(
                chrom.clone(),
                record.mapq(),
                1,
                record.seq().len() as u64,
                left as u64,
            );
            // process the main CIGAR

            record.cigar().iter().for_each(|cigar| {
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
                        right += *n;
                    }
                    Cigar::RefSkip(n) => {
                        isoform.add_segment(
                            chrom.clone(),
                            left as u64,
                            right as u64,
                            &strand,
                            false,
                        );
                        left = right + n;
                        right = left;
                    }
                    _ => {}
                }
            });
            // handle the last segment in main CIGAR
            // record.tname()
            isoform.add_segment(chrom.clone(), left as u64, right as u64, &strand, false);
            // update ref span
            isoform.update_right(right as u64);
            // process the supplementary alignment
            match record.aux("SA".as_bytes()) {
                Err(_) => {}
                Ok(sa) => {
                    if let Aux::String(sa_string) = sa {
                        sa_string
                            .split(";")
                            .collect::<Vec<&str>>()
                            .into_iter()
                            .filter(|x| x.len() > 0)
                            .for_each(|s| {
                                let sa_vec: Vec<&str> = s.split(',').collect();
                                let chrom2 = sa_vec[0].to_string();
                                let mut left = sa_vec[1].parse::<u64>().unwrap();
                                let mut right = left;
                                let strand = sa_vec[2].to_string();
                                let mut tmp_num_string = String::new();
                                // process cigar string
                                sa_vec[3].chars().into_iter().for_each(|c| {
                                    if c.is_ascii_digit() {
                                        tmp_num_string.push(c);
                                    } else if c == 'M' || c == 'D' || c == '=' || c == 'X' {
                                        right += tmp_num_string.parse::<u64>().unwrap();
                                        tmp_num_string.clear();
                                    } else if c == 'N' {
                                        isoform.add_supp_segment(
                                            chrom2.clone(),
                                            left,
                                            right,
                                            &strand,
                                            true,
                                        );
                                        left = right + tmp_num_string.parse::<u64>().unwrap();
                                        right = left;
                                        tmp_num_string.clear();
                                    } else {
                                        tmp_num_string.clear();
                                    }
                                    // handle the last segment
                                });
                                isoform.add_supp_segment(
                                    chrom2.clone(),
                                    left,
                                    right,
                                    &strand,
                                    true,
                                );
                            });
                    } else {
                        error!("SA tag is not a string: {:?}", sa);
                    }
                }
            }

            // add read name to single read
            if cli.rname {
                isoform.add_info(
                    "RN".to_string(),
                    String::from_utf8(record.qname().to_vec()).expect("can not parse read name"),
                );
            }

            isoform.process();

            match agg_isoform_map.get_mut(&isoform.signature) {
                Some(agg_isoform) => {
                    agg_isoform.add_read(isoform);
                }
                None => {
                    let agg_isoform = AggrRead::new(isoform);
                    agg_isoform_map.insert(agg_isoform.signature, agg_isoform);
                }
            }
        }
    } else if cli.gtf.is_some() && cli.bam.is_none() {
        info!("Using GTF file for isoform profiling.");

        let gtf_path = cli.gtf.as_ref().unwrap();

        info!("Loading GTF file...");
        // let gtfreader: noodles_gtf::Reader<BufReader<std::fs::File>> = noodles_gtf::io::Reader::new(
        //     BufReader::new(std::fs::File::open(gtf_path).expect("can not read gtf")),
        // );

        let gtfreader = open_gtf_reader(gtf_path.to_str().unwrap())?;

        let mut gtf = TranscriptChunker::new(gtfreader);
        let gtf_vec = gtf.get_all_transcripts_vec();
        info!("Processing transcripts from GTF...");

        for rec in gtf_vec.iter() {
            chrom_set.insert(rec.chrom.clone());

            total_count += 1;
            if total_count % batch_size == 0 {
                info!(
                    "Processing {} records",
                    total_count.to_formatted_string(&Locale::en)
                );
            }

            let mut isoform = rec.to_single_read();

            if cli.tid {
                isoform.add_info("TID".to_string(), rec.trans_id.clone());
            }
            if cli.gid {
                isoform.add_info("GID".to_string(), rec.gene_id.clone());
            }

            match agg_isoform_map.get_mut(&isoform.signature) {
                Some(agg_isoform) => {
                    agg_isoform.add_read(isoform);
                }
                None => {
                    let agg_isoform = AggrRead::new(isoform);
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
        .expect("can not write headers..");

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
            .expect("can not write record...");
    }
    info!(
        "Total records processed: {}",
        total_count.to_formatted_string(&Locale::en)
    );

    info!(
        "Total valid records: {}",
        n_record.to_formatted_string(&Locale::en)
    );
    // info!("Total records : {}", n_record);
    info!(
        "Total records skipped: {}",
        skipped_records.to_formatted_string(&Locale::en)
    );
    info!("The skipped records are due to unmapped, low mapping quality(--mapq), or secondary alignments(--use-secondary).");
    info!("Finished");
    Ok(())
}
