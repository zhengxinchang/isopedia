use crate::chromosome::ChromMapping;
use crate::constants::*;
#[allow(dead_code, unused)]
use crate::isoformarchive::read_record_from_mmap;
use crate::logger::init_logger;
use crate::reads::AggrRead;
#[allow(dead_code, unused)]
use crate::tmpidx::Tmpindex;
use crate::tools::ToolCmdValidate;
use anyhow::Result;
use clap::{Parser, Subcommand};
// use isopedia::tools::ToolCmdValidate;
use crate::writer::MyGzWriter;
use log::{error, info, warn};
#[allow(dead_code, unused)]
use memmap2::Mmap;
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

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

impl ToolCmdValidate for InspectArgs {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.idx.exists() {
            error!("--idx: index directory does not exist");
            is_ok = false;
        }

        if self.type_f != "tmpidx"
            && self.type_f != "archive"
            && self.type_f != "dbinfo"
            && self.type_f != "chroms"
        {
            error!("--type: type must be either 'tmpidx' or 'archive' or 'dbinfo'");
            is_ok = false;
        }

        is_ok
    }
}

#[allow(dead_code, unused)]
pub fn inspect_intrim_file(idx: &PathBuf, output: &PathBuf) {
    // let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    // let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    // let chromamp = ChromMapping::decode(&chrom_bytes);
    // // let blocks = idx.get_blocks(chrom_id);
    // let writer = std::fs::File::create(output).unwrap();
    // let mut writer = std::io::BufWriter::new(writer);
    // for chrom_id in chromamp.get_chrom_idxs() {
    //     let chrom_name = chromamp.get_chrom_name(chrom_id);
    //     let blocks = interim.get_blocks(chrom_id);
    //     if blocks.len() == 0 {
    //         continue;
    //     }
    //     for block in blocks {
    //         for record in block {
    //             writeln!(
    //                 writer,
    //                 "{}:{}-{:?}",
    //                 chrom_name, &record.pos, record.record_ptr_vec
    //             )
    //             .unwrap();
    //         }
    //     }
    // }
}

#[allow(dead_code, unused)]
pub fn inspect_archive(idx: &PathBuf, output: &PathBuf) {
    // let mut processed_signautres = std::collections::HashSet::new();

    // let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    // let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    // let chromamp = ChromMapping::decode(&chrom_bytes);
    // let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
    //                                                        // let blocks = idx.get_blocks(chrom_id);
    // let writer = std::fs::File::create(output).unwrap();
    // let mut writer = std::io::BufWriter::new(writer);

    // let archive_file_handler =
    //     File::open(idx.join(MERGED_FILE_NAME)).expect("Can not open aggregated records file...");
    // let archive_mmap = unsafe { Mmap::map(&archive_file_handler).expect("Failed to map the file") };

    // archive_mmap
    //     .advise(memmap2::Advice::Sequential)
    //     .expect("Failed to set mmap advice");

    // for chrom_id in chromamp.get_chrom_idxs() {
    //     // let chrom_name = chromamp.get_chrom_name(chrom_id);
    //     let blocks = interim.get_blocks(chrom_id);
    //     if blocks.len() == 0 {
    //         continue;
    //     }
    //     for block in blocks {
    //         for record_grp in block {
    //             for record_ptr in record_grp.record_ptr_vec {
    //                 let rec = read_record_from_mmap(&archive_mmap, &record_ptr, &mut archive_buf);
    //                 if processed_signautres.contains(&rec.signature) {
    //                     continue;
    //                 }
    //                 processed_signautres.insert(rec.signature);
    //                 writeln!(writer, "{}", rec.get_string()).unwrap();
    //             }
    //         }
    //     }
    // }
}

pub fn inspect_meta(idx: &PathBuf) {
    let dataset_info =
        crate::dataset_info::DatasetInfo::load_from_file(&idx.join(DATASET_INFO_FILE_NAME));
    dbg!(&dataset_info);
}

pub fn inspect_chroms(idx: &PathBuf) {
    let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    dbg!(&chromamp);
}
