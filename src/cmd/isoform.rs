use std::{hash::Hash, path::PathBuf, vec};

use crate::{
    // assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    global_stats::GlobalStats,
    grouped_tx::{ChromGroupedTxManager, TxAbundance, MSJC},
    gtf::{open_gtf_reader, TranscriptChunker},
    io::{DBInfos, Header, Line, SampleChip},
    isoform::MergedIsoform,
    isoformarchive::ArchiveCache,
    meta::Meta,
    results::TableOutput,
    tmpidx::MergedIsoformOffsetPtr,
    utils::log_mem_stats,
};
use ahash::HashSet;
use anyhow::Result;
use clap::{command, Parser};
use log::{error, info};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
// use nix::sys::mman::{posix_madvise, PosixMadvise};
// Removed incorrect import of mmap::{PosixMadvise,posix_madvise}

use once_cell::sync::OnceCell;
use rayon::prelude::*; // 如果需要

pub static GLOBAL_POOL: OnceCell<rayon::ThreadPool> = OnceCell::new();

pub fn init_thread_pool(n_threads: usize) {
    GLOBAL_POOL
        .set(
            rayon::ThreadPoolBuilder::new()
                .num_threads(n_threads)
                .build()
                .unwrap(),
        )
        .ok();
}

#[derive(Parser, Debug, Serialize, Clone)]
#[command(name = "isopedia isoform")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[Query] Annotate provided gtf file(transcripts/isoforms) with the index.
", long_about = None)]
#[clap(after_long_help = "

# Isopedia isoform needs the input gtf file to be sorted. use the following command to sort the gtf file:
gffread -T -o- input.gtf  | sort -k1,1 -k4,4n | gffread - -o sorted.gtf

")]
pub struct AnnIsoCli {
    /// Path to the index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// Path to the GTF file
    #[arg(short, long)]
    pub gtf: PathBuf,

    /// Flanking size (in bases) before and after the position
    #[arg(short, long, default_value_t = 10)]
    pub flank: u64,

    /// Minimum number of reads required to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// Output file for search results
    #[arg(short, long)]
    pub output: PathBuf,

    /// Whether to include additional information in the output
    #[arg(long, default_value_t = false)]
    pub info: bool,

    /// Include in-complete splice junction matches
    #[arg(long, default_value_t = false)]
    pub asm: bool,

    /// Max EM iterations
    #[arg(long, default_value_t = 50)]
    pub em_iter: usize,

    // EM convergence threshold
    #[arg(long, default_value_t = 0.01)]
    pub em_converge: f32,

    /// Maximum number of cached tree nodes in memory
    #[arg(short = 'c', long = "cached-nodes", default_value_t = 10)]
    pub cached_nodes: usize,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_num: usize,

    /// Cached isoform chunk size in Mb
    #[arg(long, default_value_t = 128)]
    pub cached_chunk_size_mb: u64,

    /// Verbose mode
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
}

impl AnnIsoCli {
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
                "--idxdir: merged isoform data {} does not exist in {}, please run `isopedia merge` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join("bptree_0.idx").exists() {
            error!(
                "--idxdir: tree files does not exist in {}, please run `isopedia index` first",
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.gtf.exists() {
            error!("--gtf: gtf file {} does not exist", self.gtf.display());
            is_ok = false;
        }

        if is_ok != true {
            std::process::exit(1);
        }
    }
}

fn greetings(args: &AnnIsoCli) {
    // eprintln!("\nCommand: isoform = \n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

pub fn run_anno_isoform(cli: &AnnIsoCli) -> Result<()> {
    // env::set_var("RUST_LOG", "info");
    // env_logger::init();

    greetings(&cli);
    cli.validate();

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    info!("Loading GTF file...");

    init_thread_pool(4);

    let gtfreader = open_gtf_reader(cli.gtf.to_str().unwrap())?;

    let mut gtf = TranscriptChunker::new(gtfreader);
    // let gtf_vec = gtf.get_all_transcripts_vec();

    let gtf_by_chrom = gtf.get_all_transcripts_by_chrom();

    info!(
        "Loaded {} transcripts from gtf file",
        gtf.trans_count.to_formatted_string(&Locale::en)
    );

    info!("Loading index file");

    let mut archive_cache = ArchiveCache::new(
        cli.idxdir.clone().join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb * 1024 * 1024, // 512MB chunk size
        cli.cached_chunk_num,                   // max 4 chunks in cache ~2GB
    );

    let mut global_stats = GlobalStats::new(dataset_info.get_size());

    let meta = Meta::parse(&cli.idxdir.join(META_FILE_NAME), None)?;
    let mut out_header = Header::new();
    out_header.add_column("chrom")?;
    out_header.add_column("start")?;
    out_header.add_column("end")?;
    out_header.add_column("length")?;
    out_header.add_column("exon_count")?;
    out_header.add_column("trans_id")?;
    out_header.add_column("gene_id")?;
    out_header.add_column("confidence")?;
    out_header.add_column("detected")?;
    out_header.add_column("min_read")?;
    out_header.add_column("positive_count/sample_size")?;
    out_header.add_column("attributes")?;
    let mut db_infos = DBInfos::new();
    for (name, evidence) in dataset_info.get_sample_evidence_pair_vec() {
        db_infos.add_sample_evidence(&name, evidence);
        out_header.add_sample_name(&name)?;
    }

    let mut tableout = TableOutput::new(
        cli.output.clone(),
        out_header.clone(),
        db_infos.clone(),
        meta.clone(),
        "CPM:COUNT:FSM_CPM:FSM_COUNT:EM_CPM:EM_COUNT:INFO".to_string(),
    );

    info!("Initializing transcript groups");

    let mut chrom_grouped_tx_managers: Vec<ChromGroupedTxManager> = gtf_by_chrom
        .into_iter()
        .map(|(chrom, tx_vec)| {
            let mut manager = ChromGroupedTxManager::new(&chrom, dataset_info.get_size());
            manager.add_transcript_by_chrom(&tx_vec);
            manager
        })
        .collect();

    info!("Processing transcripts");
    for chrom_manager in chrom_grouped_tx_managers.iter_mut() {
        chrom_manager.process_tx_groups(
            &mut forest,
            cli,
            &mut archive_cache,
            &dataset_info,
            &mut tableout,
            &mut global_stats,
        );
    }

    /*
    For each chromsome:
    1. obtain all the splice junctions from the transcripts
    2. group all transcripts by their splice junctions into tx_groups

    no need to do this per sj, just find all tx that overlap with each other is fine to group them.

    3. for each tx_group, obtain all the misoform offsets from the index. dedup and dump the offsets into OffsetCache but keep the group_id reference. maintain a mapping from offset --> set[group_id,]
    4. merge tx_groups that share same misoform offsets, update the OffsetCache accordingly.
        dont udpate it, keep it simple, but give the API to improve, global storage of msiformoffset, even by chrom, would take a huge mount of memory
        make an API for further disk-persistent optimization,but make group seperated for now.
    5. for each tx_group, create TxAbundance records, load the misoform from OffsetCache.
        if a misoform matches all the splice junctions of the group,
            load into TxAbundance as FSM
        else
            create MSJC records for those misoforms that partially match the splice junctions, and load into TxAbundance
    6. after all TxAbundance and MSJC are created, run EM algorithm to estimate the abundances
    7. output the TxAbundance results


    multi-threading consideration:

    1. per-chromsome processing, each thread process one chromsome at a time， send a grouptx to  a thread pool.
    2. the thread pool process each group tx independently with EM and send back to main thread for output
    3. the main thread collect results and output

    ========
       splice junction --> inital transcript groups

       init projection TxComponent  按照chrosome 分组，每次处理一个chromsome的内容
       tx --> group_id
       group_id --> tx_ids // used for update gruop id later

       init projection OffsetComponent
       offset --> set[group_id,]

       我既然获得了所有的sj，我其实可以顺序查询，而不是每次都是一个transcript查询了

       一个group的tx 一次查询，返回所有的misoformoffset


       Vec（
           （position，vec（misoform——offset））
       )

       每个transcript，查看这些position，然后找到对应的FSM对应的FSM transcript

       剩下的就是用于创建MSJC的misofromoffset, 这些可以dump到一个临时文件OffsetCache中。
       Txabundance应该包含对MSJC的misoformoffset的index引用，这样做EM的时候可以直接load进来？

       OffsetCache {
           group_id --> offset，length记录
       }


       如果一个offset 出现在多个gorup中，那么需要合并gruop，
       实际的合并发生在对TxGruop的操作,也发生在对OffsetCache的操作

       到这里，就可以获得所有的 MSJC的offset了， 这里可能会存在重复的offset，就是不同的group包含了同一个offset，这种情况会
       多产生几次isoformarchive的IO不影响结果

       此时所有的Txabundance 已经创建好，并且有了


    */
    Ok(())
}
