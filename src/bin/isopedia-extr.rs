use clap::Parser;
use core::panic;
use indexmap::IndexSet;
use log::{error, info};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar},
    Read, Record,
};
use rustc_hash::FxHashMap;
use std::env;
use std::path::PathBuf;

use isopedia::{
    reads::{AggrRead, SingleRead},
    writer::MyGzWriter,
};
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
#[command(name = "isopedia-extr")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com,Xinchang.Zheng@bcm.edu>
", long_about = None)]
struct Cli {
    /// Input file in BAM/CRAM format
    #[arg(short, long)]
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

    /// Include secondary reads in the analysis
    #[arg(long = "use-secondary")]
    pub use_secondary: bool,

    /// Debug mode
    #[arg(long)]
    pub debug: bool,
}

impl Cli {
    fn validate(&self) {
        let mut is_ok = true;
        if !self.bam.exists() {
            error!("--bam: input file {} does not exist", self.bam.display());
            is_ok = false;
        }

        // let output_dir = self
        //     .output
        //     .parent()
        //     .unwrap_or_else(|| std::path::Path::new("."));

        // if !output_dir.is_dir() {
        //     error!(
        //         "--output: parent dir {} does not exist",
        //         output_dir.display()
        //     );
        //     is_ok = false;
        // }

        // check output is not a directory
        if self.output.is_dir() {
            error!(
                "--output: {} is a directory, not a file",
                self.output.display()
            );
            is_ok = false;
        }

        // optional: prevent silent overwrite unless forced
        // if self.output.exists() {
        //     error!(
        //         "--output: file {} already exists. Consider removing it or adding --force",
        //         self.output.display()
        //     );
        //     is_ok = false;
        // }

        if is_ok != true {
            panic!("Invalid arguments, please check the error messages above.");
        }
    }
}

fn greetings(args: &Cli) {
    println!("\nIsopedia: [Extract raw isoform singals from BAM/CRAM]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => println!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let bam_path = &cli.bam;

    let mut bam_reader = bam::Reader::from_path(bam_path).expect("Failed to open BAM file");

    if &cli.bam.ends_with(".cram") & &cli.reference.is_some() {
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
    let mut chrom_set: IndexSet<String> = IndexSet::new();
    let mut agg_isoform_map: FxHashMap<u64, AggrRead> = FxHashMap::default();

    let mut mywriter = MyGzWriter::new(&cli.output).expect(&format!(
        "Can not create output file {} .",
        &cli.output.display()
    ));

    let mut batches = 0;
    let mut curr_count = 0;
    let batch_size = 1000000; // 1 million records per batch
    let mut skipped_records = 0;

    while let Some(result) = bam_reader.read(&mut record) {
        
        curr_count += 1;

        match result {
            Err(_) => {
                eprintln!("Can not read record");
                break;
            } //exit if the last one was processed.
            Ok(()) => {}
        }

        if curr_count % batch_size == 0 {
            info!("Processing {} records", batches * batch_size);
            batches += 1;
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

        // dont use record.tid() before it is filled by bam.read()
        // otherwise it will be random number and cuase the error with random memory access.
        let mut chrom = String::from_utf8(header.tid2name(record.tid() as u32).to_vec()).unwrap();
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
                    isoform.add_segment(chrom.clone(), left as u64, right as u64, &strand, false);
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
        isoform.update_ref_span(right as u64);
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
                            isoform.add_supp_segment(chrom2.clone(), left, right, &strand, true);
                        });
                } else {
                    error!("SA tag is not a string: {:?}", sa);
                }
            }
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

    // write the chrom_set
    info!("Process chromsomes");
    let mut chrom_str = Vec::new();
    chrom_set.into_iter().for_each(|x| {
        chrom_str.extend_from_slice(x.as_bytes());
        chrom_str.push(b'\t');
        // writer.write_all(x.as_bytes()).unwrap();
        // writer.write_all(b"\t").unwrap();
    });
    chrom_str.push(b'\n');
    // writer.write_all(b"\n").unwrap();

    info!("Write to output file");
    // writer
    //     .write(b"signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n")
    //     .unwrap();

    chrom_str.extend_from_slice(b"signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n");
    mywriter
        .write_all_bytes(&chrom_str)
        .expect("can not write headers..");

    if cli.debug {
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
        batches * batch_size + curr_count
    );
    info!("Total records skipped: {}", skipped_records);
    info!("Finished");
}
