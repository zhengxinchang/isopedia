use anyhow::Result;
use clap::{Parser, Subcommand};
use isopedia::chromosome::ChromMapping;
use isopedia::constants::*;
use isopedia::isoformarchive::read_record_from_mmap;
use isopedia::reads::AggrRead;
use isopedia::tmpidx::Tmpindex;
use isopedia::writer::MyGzWriter;
use log::{error, info, warn};
use memmap2::Mmap;
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

pub fn inspect_intrim_file(idx: &PathBuf, output: &PathBuf) {
    let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);
    for chrom_id in chromamp.get_chrom_idxs() {
        let chrom_name = chromamp.get_chrom_name(chrom_id);
        let blocks = interim.get_blocks(chrom_id);
        if blocks.len() == 0 {
            continue;
        }
        for block in blocks {
            for record in block {
                writeln!(
                    writer,
                    "{}:{}-{:?}",
                    chrom_name, &record.pos, record.record_ptr_vec
                )
                .unwrap();
            }
        }
    }
}

pub fn inspect_archive(idx: &PathBuf, output: &PathBuf) {
    let mut processed_signautres = std::collections::HashSet::new();

    let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
                                                           // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);

    // let mut aggr_reader = std::io::BufReader::new(
    //     std::fs::File::open(idx.join(MERGED_FILE_NAME))
    //         .expect("Can not open aggregated records file...exit"),
    // );

    let archive_file_handler = File::open(idx.join(MERGED_FILE_NAME))
        .expect("Can not open aggregated records file...");
    let archive_mmap = unsafe { Mmap::map(&archive_file_handler).expect("Failed to map the file") };

    archive_mmap.advise(memmap2::Advice::Sequential).expect("Failed to set mmap advice");

    for chrom_id in chromamp.get_chrom_idxs() {
        // let chrom_name = chromamp.get_chrom_name(chrom_id);
        let blocks = interim.get_blocks(chrom_id);
        if blocks.len() == 0 {
            continue;
        }
        for block in blocks {
            for record_grp in block {
                for record_ptr in record_grp.record_ptr_vec {
                    let rec =
                        read_record_from_mmap(&archive_mmap, &record_ptr, &mut archive_buf);
                    if processed_signautres.contains(&rec.signature) {
                        continue;
                    }
                    processed_signautres.insert(rec.signature);
                    writeln!(writer, "{}", rec.get_string()).unwrap();
                }
            }
        }
    }
}

fn inspect_meta(idx: &PathBuf) {
    let dataset_info =
        isopedia::dataset_info::DatasetInfo::load_from_file(&idx.join(DATASET_INFO_FILE_NAME));
    dbg!(&dataset_info);
}

fn merge_replicates(files: &Vec<PathBuf>, output: &PathBuf) -> Result<()> {
    let mut chrom_set: HashSet<String> = HashSet::new();
    let mut total_count = 0;

    if files.len() == 0 || files.len() == 1 {
        error!("At least two files are required to merge replicates");
        //exit
        std::process::exit(1);
    }

    for file in files {
        if !file.exists() {
            error!("File {:?} does not exist", file);
            std::process::exit(1);
        }
    }

    let base_file = files[0].clone();

    //open gz
    let mut buf = String::new();
    let mut gzrdr = BufReader::new(flate2::read::GzDecoder::new(std::fs::File::open(
        base_file,
    )?));
    // merge chroms
    gzrdr.read_line(&mut buf)?;
    let chroms = buf.trim().split('\t').collect::<Vec<&str>>();

    for chrom in chroms {
        chrom_set.insert(chrom.to_string());
    }

    // skip the second line(header line)
    gzrdr.read_line(&mut buf)?;

    // start to build the map of agg-read

    let mut merged_map = std::collections::HashMap::new();

    buf.clear();
    while let Ok(n) = gzrdr.read_line(&mut buf) {
        if n == 0 {
            break;
        }

        total_count += 1;

        let agg_read = AggrRead::from_string(&buf);
        merged_map.insert(agg_read.signature.clone(), agg_read);
        buf.clear();
    }

    // iteraterlly add agg_read from the rest of files

    for file in &files[1..] {
        let mut buf = String::new();
        let mut gzrdr = BufReader::new(flate2::read::GzDecoder::new(std::fs::File::open(file)?));

        // read first line for choromsome
        gzrdr.read_line(&mut buf)?;
        let chroms = buf.trim().split('\t').collect::<Vec<&str>>();
        for chrom in chroms {
            if !chrom_set.contains(chrom) {
                warn!("Chromosome {} is not present in all input files, this might happen in non-main contigs, please be careful if you see this message for any main chromosomes.", chrom);
            }
            chrom_set.insert(chrom.to_string());
        }

        // skip the header
        gzrdr.read_line(&mut buf)?;
        buf.clear();
        // start to build the map of agg-read
        while let Ok(n) = gzrdr.read_line(&mut buf) {
            if n == 0 {
                break;
            }

            total_count += 1;

            let agg_read = AggrRead::from_string(&buf);

            merged_map
                .entry(agg_read.signature.clone())
                .and_modify(|e| e.merge(agg_read.clone()))
                .or_insert(agg_read);

            buf.clear();
        }
    }

    info!("{} Chromsomes detected and merged.", chrom_set.len());

    // write the chrom_set
    info!("Process chromsomes");
    let mut chrom_str = Vec::new();
    chrom_set.into_iter().for_each(|x| {
        chrom_str.extend_from_slice(x.as_bytes());
        chrom_str.push(b'\t');
    });
    chrom_str.push(b'\n');

    info!("Write to output file");

    let mut mywriter = MyGzWriter::new(output).expect(&format!(
        "Can not create output file {} .",
        output.display()
    ));

    info!(
        "Merged {} records into {} agg-reads from {} files",
        total_count,
        &merged_map.len(),
        files.len()
    );

    chrom_str.extend_from_slice(b"signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n");
    mywriter
        .write_all_bytes(&chrom_str)
        .expect("can not write headers..");

    let mut sorted_agg_isoform_map = merged_map.into_iter().collect::<Vec<(u64, AggrRead)>>();
    sorted_agg_isoform_map.sort_by_key(|x| x.0);
    for (_, agg_isoform) in sorted_agg_isoform_map.iter() {
        // writer.write(agg_isoform.to_record().as_bytes()).unwrap();
        mywriter
            .write_all_bytes(agg_isoform.to_record().as_bytes())
            .expect("can not write record...");
    }

    info!("Finished");

    Ok(())
}

#[derive(Parser, Clone, Debug, Serialize, Deserialize)]
#[command(name = "isopedia-tool")]
#[command(about = "Auxiliary tools for Isopedia")]
#[command(author = "Xinchang Zheng", version)]
pub struct ToolCli {
    #[command(subcommand)]
    pub command: ToolsCommands,
}

#[derive(Subcommand, Clone, Debug, Serialize, Deserialize)]
pub enum ToolsCommands {
    #[allow(non_camel_case_types)]
    #[command(about = "Inspect index files")]
    inspect(InspectArgs),
    #[allow(non_camel_case_types)]
    #[command(about = "Merge replicates file from isopedia-extr")]
    merge(MergeArgs),
}

pub trait Validate {
    fn validate(&self) -> bool;
}

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "

Convert the binary index files into plain text format.

Examples: 

isopedia-tool inspect --idx <IDX> -t [tmpidx|archive|dbinfo] --output <OUTPUT>

For the -f paramter:
tmpidx: Inspect the temporary index files
archive: Inspect the archive files
dbinfo: Inspect the database information file

")]
pub struct InspectArgs {
    /// index directory
    #[arg(short, long)]
    pub idx: PathBuf,

    /// type of file to be inspect can be either 'tmpidx' or 'archive' or "dbinfo"
    #[arg(short = 't', long, name = "type")]
    pub type_f: String,

    /// Output file name
    #[arg(short, long)]
    pub output: PathBuf,
}

impl Validate for InspectArgs {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.idx.exists() {
            error!("--idx: index directory does not exist");
            is_ok = false;
        }

        if self.type_f != "tmpidx" && self.type_f != "archive" && self.type_f != "dbinfo" {
            error!("--type: type must be either 'tmpidx' or 'archive' or 'dbinfo'");
            is_ok = false;
        }

        is_ok
    }
}

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "

Merge replicate files from isopedia-extr into a single file.

This is useful when you have multiple replicate files that you want to combine into a single file.

Example:

isopedia-tool merge -i <INPUT_FILE 1> <INPUT_FILE 2> -o <OUTPUT>


At least two input files must be specified.
")]
pub struct MergeArgs {
    /// input files that generated from isopedia-extr
    #[arg(short, long, num_args = 2..)]
    pub input_files: Vec<PathBuf>,

    /// output merged file in gz format
    #[arg(short, long)]
    pub output: PathBuf,
}

impl MergeArgs {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if self.input_files.len() < 2 {
            error!("--input: at least two input files must be specified");
            is_ok = false;
        }

        is_ok
    }
}

fn greetings(args: &ToolCli) {
    println!("\nIsopedia: [Extract raw isoform singals from BAM/CRAM]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => println!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() {
    env_logger::Builder::from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    )
    .format(|buf, refcord| {
        writeln!(
            buf,
            "[{}] [{}]: {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            refcord.level(),
            refcord.args()
        )
    })
    .init();

    let cli = ToolCli::parse();

    greetings(&cli);

    match cli.command {
        ToolsCommands::inspect(inspec_args) => {
            if !inspec_args.validate() {
                std::process::exit(1);
            }

            if inspec_args.type_f == "tmpidx" {
                inspect_intrim_file(&inspec_args.idx, &inspec_args.output);
            } else if inspec_args.type_f == "archive" {
                inspect_archive(&inspec_args.idx, &inspec_args.output);
            } else if inspec_args.type_f == "dbinfo" {
                inspect_meta(&inspec_args.idx);
            }
        }
        ToolsCommands::merge(merge_args) => {
            if !merge_args.validate() {
                std::process::exit(1);
            }

            merge_replicates(&merge_args.input_files, &merge_args.output)
                .expect("Failed to merge replicates");
        }
    }
}
