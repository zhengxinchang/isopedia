use anyhow::Result;
use clap::{arg, Parser};
use isopedia::bptree::BPForest;
use isopedia::breakpoints::{self, BreakPointsPair};
use isopedia::dataset_info::DatasetInfo;
use isopedia::utils::{get_total_memory_bytes, warmup};
use isopedia::writer::MyGzWriter;
use isopedia::{constants::*, utils};
use log::{error, info};
use serde::Serialize;
use std::path::PathBuf;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-ann-splice")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>, <xinchang.zheng@bcm.edu>
", long_about = None)]
#[clap(after_long_help = "
Use 
")]
struct Cli {
    /// Path to the index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// Splice junction in 'chr1:pos1,chr2:pos2' format
    #[arg(short = 'i', long = "splice")]
    pub splice: Option<String>,

    /// Path to splice junction bed file
    #[arg(short = 'I', long = "splice-bed")]
    pub splice_bed: Option<PathBuf>,

    /// Flanking size (in bases) before and after the position
    #[arg(short, long, default_value_t = 10)]
    pub flank: u64,

    /// Minimum number of reads required to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// Output file for search results
    #[arg(short, long)]
    pub output: PathBuf,

    /// Memory size to use for warming up (in gigabytes).  
    /// Example: 1GB. Increasing this will significantly improve performance;  
    /// set it as large as your system allows.
    #[arg(short, long, default_value_t = 4)]
    pub warmup_mem: usize,

    /// Maximum number of cached nodes per tree
    #[arg(short = 'c', long = "cached_nodes", default_value_t = 100_000)]
    pub lru_size: usize,
}

impl Cli {
    fn validate(&self) {
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
                "--idxdir: Aggr file {} does not exist in {}, please run `stix-isoform aggr` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join("bptree_0.idx").exists() {
            error!(
                "--idxdir: index files {} does not exist in {}, please run `stix-isoform idx` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if let Some(ref bed_path) = self.splice_bed {
            if !bed_path.exists() {
                error!(
                    "--splice-bed: splice bed file {} does not exist",
                    bed_path.display()
                );
                is_ok = false;
            }
        }

        let max_mem_bytes = match get_total_memory_bytes() {
            Some(bytes) => bytes,
            None => {
                error!("Failed to get total memory bytes");
                std::process::exit(1);
            }
        };

        if let Some(splice_str) = &self.splice {
            if utils::parse_breakpoint_str(&splice_str).is_err() {
                error!("Failed to parse splice junction: {}", splice_str);
                std::process::exit(1);
            }
        }

        if self.warmup_mem == 0 {
            error!("--warmup-mem: must be greater than 0");
            is_ok = false;
        } else {
            let warmup_bytes = (self.warmup_mem as u64) * 1024 * 1024 * 1024;
            if warmup_bytes > max_mem_bytes {
                error!(
                    "--warmup-mem: {} GB larger than system memory, please set it to less than {} GB",
                    self.warmup_mem,
                    max_mem_bytes / (1024 * 1024 * 1024)
                );
                is_ok = false;
            }
        }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &Cli) {
    eprintln!("\nIsopedia: [Annotate provided gtf file]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    info!("Warmup index file");
    let max_gb = cli.warmup_mem * 1024 * 1024 * 1024;
    warmup(&cli.idxdir.clone().join(MERGED_FILE_NAME), max_gb)?;

    info!("loading indexes");
    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    let mut hit_count = 0u32;
    let mut miss_count = 0u32;

    // init the output writer
    let mut mywriter = MyGzWriter::new(&cli.output)?;
    let mut header_str =
        String::from("chr1\tpos1\tchr2\tpos2\tid\trest_info\tpositive/sample_size");
    let sample_name = dataset_info.get_sample_names();
    for name in sample_name {
        header_str.push_str(&format!("\t{}", name));
    }
    header_str.push('\n');
    mywriter.write_all_bytes(header_str.as_bytes())?;

    let mut queries: Vec<BreakPointsPair> = if cli.splice.is_some() {
        info!("parse breakpoints pair from command line...");
        let splice_str = cli.splice.clone().unwrap();
        let mut bp =
            BreakPointsPair::parse_string(&splice_str).expect("Can not parse splice junction...");
        bp.sort();
        vec![bp]
    } else if cli.splice_bed.is_some() {
        let splice_bed_path = cli.splice_bed.unwrap();
        info!(
            "parse breakpoints pairs from file {}",
            &splice_bed_path.display()
        );

        let bpv = breakpoints::bed2breakpointsvec(&splice_bed_path).expect(&format!(
            "Can not parse splice junctions from {}",
            &splice_bed_path.display()
        ));
        bpv
    } else {
        error!("Please provide either --splice or --splice-bed");
        std::process::exit(1);
    };

    // validate the breakpoints from same chromsome
    for breakpoint_pair in &queries {
        if !breakpoint_pair.is_same_chr() {
            error!(
                "splice junction should be in same chromosome: {}:{} vs {}:{}",
                &breakpoint_pair.left_chr,
                &breakpoint_pair.left_pos,
                &breakpoint_pair.right_chr,
                &breakpoint_pair.right_pos
            );
            std::process::exit(1);
        }
    }

    for query in &queries {}

    Ok(())
}
