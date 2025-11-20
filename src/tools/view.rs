use crate::chromosome::ChromMapping;
use crate::constants::*;
use crate::isoformarchive::ArchiveCache;
use crate::tmpidx::Tmpindex;
use crate::tools::ToolCmdValidate;
use clap::Parser;
// use isopedia::tools::ToolCmdValidate;
use crate::io::MyGzWriter;
use log::error;
// #[allow(dead_code, unused)]
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "

Convert the binary index files into plain text format.

Examples: 

isopedia-tool inspect --idx <IDX> -t [tmpidx|archive|dbinfo|chroms] --output <OUTPUT>

For the -t parameter:
tmpidx: View the temporary index files
archive: View the archive files
dbinfo: View the database information file
chroms: View the chromosome mapping file

")]
pub struct ViewArgs {
    /// index directory
    #[arg(short, long)]
    pub idx: PathBuf,

    /// type of file to be inspect can be either 'tmpidx' or 'archive' or "dbinfo"
    #[arg(short = 't', long, name = "type")]
    pub type_f: String,

    /// Output file name
    #[arg(short, long)]
    pub output: PathBuf,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_number: usize,

    /// Cached isoform chunk size in Mb
    #[arg(long, default_value_t = 128)]
    pub cached_chunk_size_mb: u64,
}

impl ToolCmdValidate for ViewArgs {
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
pub fn view_tmpidx(cli: &ViewArgs) {
    let mut tmpidx = Tmpindex::load(&cli.idx.join(TMPIDX_FILE_NAME));
    let chrom_bytes = std::fs::read(&cli.idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(&cli.output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);
    for chrom_id in chromamp.get_chrom_idxs() {
        let chrom_name = chromamp.get_chrom_name(chrom_id);
        let blocks = tmpidx
            .groups(chrom_id)
            .expect("Failed to get blocks from tmpidx");
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

#[allow(dead_code, unused)]
pub fn view_archive(cli: &ViewArgs) {
    let mut processed_signautres = std::collections::HashSet::new();

    // let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    let chrom_bytes = std::fs::read(&cli.idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    // let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
    //                                                        // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(&cli.output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);

    // let archive_file_handler = File::open(cli.idx.join(MERGED_FILE_NAME))
    //     .expect("Can not open aggregated records file...");
    // let archive_mmap = unsafe { Mmap::map(&archive_file_handler).expect("Failed to map the file") };

    // archive_mmap
    //     .advise(memmap2::Advice::Sequential)
    //     .expect("Failed to set mmap advice");

    let mut archive = ArchiveCache::new(
        &cli.idx.join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb,
        cli.cached_chunk_number,
    );

    let mut mygzwriter = MyGzWriter::new(&cli.output).expect("Failed to create MyGzWriter");

    for chrom_id in chromamp.get_chrom_idxs() {
        let cache_name = &cli.idx.join(format!("bptree_{}.idx", chrom_id));
        let mut cache = crate::bptree::Cache::from_disk(cache_name, 10)
            .expect("Failed to load cache from disk");

        let leafs = cache.get_leaf_notes();

        for leaf in leafs {
            for ptr in &leaf.data.merge_isoform_offset_vec {
                // let rec = read_record_from_mmap(&archive_mmap, &ptr, &mut archive_buf);

                let rec = archive.read_bytes(ptr);
                if processed_signautres.contains(&rec.signature) {
                    continue;
                }
                processed_signautres.insert(rec.signature);
                let mut s = rec.get_string();
                s.push('\n');
                mygzwriter
                    .write_all_bytes(&s.into_bytes())
                    .expect("Failed to write record");
            }
        }
    }
}

pub fn view_dbinfo(view_args: &ViewArgs) {
    let dataset_info = crate::dataset_info::DatasetInfo::load_from_file(
        &view_args.idx.join(DATASET_INFO_FILE_NAME),
    )
    .expect("Failed to load dataset info file");
    dataset_info
        .save_to_file(&view_args.output)
        .expect("Failed to save dataset info file");
}

pub fn view_chroms(view_args: &ViewArgs) {
    let chrom_bytes = std::fs::read(&view_args.idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    dbg!(&chromamp);
}
