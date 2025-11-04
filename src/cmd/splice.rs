use crate::bptree::BPForest;
use crate::breakpoints::{self, BreakPointPair};
use crate::dataset_info::DatasetInfo;
use crate::io::{self, DBInfos, Line};
use crate::isoform::MergedIsoform;
use crate::results::TableOutput;
use crate::{constants::*, meta, utils};
use anyhow::Result;
use clap::{arg, Parser};
use log::{error, info, warn};
use memmap2::Mmap;
use serde::Serialize;
use std::fs::File;
use std::path::PathBuf;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia splice")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
[Query/visualization] Search all isoforms overlapping with the provided splice junction(s) and visualize the results.
", long_about = None)]
#[clap(after_long_help = "

Note that if you are using the coordinates from GFF/GTF, please convert the 1-based to 0-based coordinates.

")]
pub struct AnnSpliceCli {
    /// Path to the index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// Splice junction in 'chr1:pos1,chr2:pos2' format
    #[arg(short = 's', long = "splice")]
    pub splice: Option<String>,

    /// Path to splice junction bed file
    #[arg(short = 'S', long = "splice-bed")]
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

    /// Maximum number of cached nodes per tree
    #[arg(short = 'c', long = "cached_nodes", default_value_t = 100_000)]
    pub lru_size: usize,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_number: usize,

    /// Cached isoform chunk size in Mb
    #[arg(long, default_value_t = 128)]
    pub cached_chunk_size_mb: u64,

    /// Debug mode
    #[arg(short, long, default_value_t = false, help = "Enable debug mode")]
    pub debug: bool,
}

impl AnnSpliceCli {
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

        if let Some(splice_str) = &self.splice {
            if utils::parse_splice_junction_str(&splice_str).is_err() {
                error!("Failed to parse splice junction: {}", splice_str);
                std::process::exit(1);
            }
        }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &AnnSpliceCli) {
    eprintln!("\nIsopedia: [Annotate provided gtf file]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }

    eprintln!("Note that if you are using the coordinates from GFF/GTF, please convert the 1-based to 0-based coordinates.");
}

pub fn run_anno_splice(cli: &AnnSpliceCli) -> Result<()> {
    greetings(&cli);
    cli.validate();

    info!("loading indexes");
    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    info!("loading metadata");
    let meta = meta::Meta::parse(cli.idxdir.join(META_FILE_NAME), None)?;

    // let archive_file_handle = File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
    //     .expect("Can not open aggregated records file...");

    // let archive_mmap = unsafe { Mmap::map(&archive_file_handle).expect("Failed to map the file") };

    // archive_mmap
    //     .advise(memmap2::Advice::Sequential)
    //     .expect("Failed to set mmap advice");

    let mut archive_cache = crate::isoformarchive::ArchiveCache::new(
        cli.idxdir.clone().join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb * 1024 * 1024, // chunk size
        cli.cached_chunk_number,                // number of chunks
    );

    const SPLICE_FORMAT: &str = "COUNT:CPM:START,END,STRAND";

    let mut out_header = io::Header::new();
    out_header.add_column("id")?;
    out_header.add_column("chr1")?;
    out_header.add_column("pos1")?;
    out_header.add_column("chr2")?;
    out_header.add_column("pos2")?;
    out_header.add_column("total_evidence")?;
    out_header.add_column("cpm")?;
    out_header.add_column("matched_sj_idx")?;
    out_header.add_column("dist_to_matched_sj")?;
    out_header.add_column("n_exons")?;
    out_header.add_column("start_pos_left")?;
    out_header.add_column("start_pos_right")?;
    out_header.add_column("end_pos_left")?;
    out_header.add_column("end_pos_right")?;
    out_header.add_column("splice_junctions")?;

    let queries: Vec<BreakPointPair> = if cli.splice.is_some() {
        info!("parse breakpoints pair from command line...");
        let splice_str = cli.splice.clone().unwrap();
        let bp =
            BreakPointPair::parse_string_sj(&splice_str).expect("Can not parse splice junction...");
        vec![bp]
    } else if cli.splice_bed.is_some() {
        let splice_bed_path = cli.splice_bed.clone().unwrap();
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

    info!("Found {} splice junctions for query.", queries.len());

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

    let mut dbinfo = DBInfos::new();

    for (name, evidence_count) in dataset_info.get_sample_evidence_pair_vec() {
        dbinfo.add_sample_evidence(&name, evidence_count);
        out_header.add_sample_name(&name)?;
    }

    let mut tableout = TableOutput::new(
        cli.output.clone(),
        out_header,
        dbinfo,
        meta.clone(),
        SPLICE_FORMAT.to_string(),
    );

    // make output
    let mut out_str = String::with_capacity(1024);
    for query in &queries {
        out_str.clear();
        out_str.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}",
            query.id, query.left_chr, query.left_pos, query.right_chr, query.right_pos,
        ));

        let (isoforms_ptr, _) =
            forest.search2_all_match(&query.to_pos_vec(), cli.flank, cli.lru_size);

        if isoforms_ptr.is_empty() {
            warn!("No isoforms found for query: {}", query.id);
        } else {
            for offset in &isoforms_ptr {
                let mut out_line = Line::new();
                out_line.add_field(&query.id);
                out_line.add_field(&query.left_chr);
                out_line.add_field(&query.left_pos.to_string());
                out_line.add_field(&query.right_chr);
                out_line.add_field(&query.right_pos.to_string());

                out_line.update_format_str(SPLICE_FORMAT);

                let record: MergedIsoform = archive_cache.read_bytes(&offset);
                // read_record_from_mmap(&archive_mmap, &offset, &mut archive_buf);
                if cli.debug {
                    dbg!(&record);
                }
                match record.get_splice_report(query, cli.flank, &dataset_info) {
                    Some((middle_part, sample_vec)) => {
                        for item in &middle_part {
                            out_line.add_field(item);
                        }

                        for samplec in &sample_vec {
                            out_line.add_sample(samplec.clone());
                        }
                        tableout.add_line(&out_line)?;
                    }
                    None => {}
                }
            }
        }
    }

    tableout.finish()?;

    info!("Finished!");

    Ok(())
}
