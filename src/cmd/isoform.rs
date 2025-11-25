use std::{hash::Hash, path::PathBuf, vec};

use crate::{
    assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    em::{TxAbundance, MSJC},
    gtf::{open_gtf_reader, TranscriptChunker},
    io::{DBInfos, Header, Line, SampleChip},
    isoform::MergedIsoform,
    isoformarchive::ArchiveCache,
    meta::Meta,
    results::TableOutput,
    runtime::Runtime,
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

#[derive(Parser, Debug, Serialize)]
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

    /// Maximum number of cached tree nodes in memory
    #[arg(short = 'c', long = "cached-nodes", default_value_t = 1000)]
    pub lru_size: usize,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_num: usize,

    /// Cached isoform chunk size in Mb
    #[arg(long, default_value_t = 128)]
    pub cached_chunk_size_mb: u64,

    /// Debug mode
    #[arg(long, default_value_t = false)]
    pub debug: bool,
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

    let gtfreader = open_gtf_reader(cli.gtf.to_str().unwrap())?;

    let mut gtf = TranscriptChunker::new(gtfreader);
    let gtf_vec = gtf.get_all_transcripts_vec();
    info!(
        "Loaded {} transcripts from gtf file",
        gtf_vec.len().to_formatted_string(&Locale::en)
    );

    info!("Loading sample meta...");

    let meta = Meta::parse(&cli.idxdir.join(META_FILE_NAME), None)?;

    const ISOFORM_FORMAT: &str = "CPM:COUNT:INFO";

    let mut out_header = Header::new();

    let mut all_returned_results_count = 0;

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
        ISOFORM_FORMAT.to_string(),
    );

    info!("Loading index file");

    let mut archive_cache = ArchiveCache::new(
        cli.idxdir.clone().join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb * 1024 * 1024, // 512MB chunk size
        cli.cached_chunk_num,                   // max 4 chunks in cache ~2GB
    );

    info!("Processing transcripts");
    let mut iter_count = 0;
    let mut batch = 0;

    // let mut assembler = Assembler::init(dataset_info.get_size());

    // let mut runtime = Runtime::init(dataset_info.get_size());

    let mut TxAbundance_vec: Vec<TxAbundance> = Vec::new();
    let mut curr_tx_id: usize = 0;
    let mut MSJC_vec: Vec<MSJC> = Vec::new();
    let mut curr_msjc_id: usize = 0;

    /*
    For each chromsome:
    1. obtain all the splice junctions from the transcripts
    2. group all transcripts by their splice junctions into tx_groups
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

    key data structures:

    for iterm 1 and 2:
    TxGroupManager {
        sj -> set[tx_id,]
        tx_id -> group_id
        group_id -> set[tx_id,]
    }

    // for item 3 , 4, and 5:
    // OffsetCache {
    //     offset -> MisoformOffsetRecord {
    //         length: u32,
    //         group_id_set: set[group_id,],
    //     }
    // }



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
    let mut tmp_offset_set = HashSet::default();
    for trans in gtf_vec {
        iter_count += 1;

        // runtime.add_one_total();
        // runtime.reset();

        if iter_count == 10_000 {
            batch += 1;
            info!(
                "Processed {} transcripts",
                (batch * 10_000).to_formatted_string(&Locale::en)
            );
            iter_count = 0;

            if cli.debug {
                let isoform_out_mem = &tableout.get_mem_size();
                let tree_mem = forest.get_mem_size();
                info!(
                    "Querying chromosome: {}, Current output in-memory size: {} MB, tree in-memory size: {} MB",
                    trans.chrom,
                    isoform_out_mem / (1024 * 1024),
                    tree_mem / (1024 * 1024),
                );

                log_mem_stats();
            }
        }

        let mut tx_abundance = TxAbundance::new(curr_tx_id, dataset_info.get_size(), &trans);
        curr_tx_id += 1;

        let mut queries: Vec<(String, u64)> = trans.get_quieries();
        queries.sort_by_key(|x| x.1);

        let (exact_matches, position_based_matches) =
            forest.search2_all_match(&queries, cli.flank, cli.lru_size);

        // load fsm in to tx_abundance
        for ext_match in &exact_matches {
            if ext_match.n_splice_sites as usize != queries.len() {
                continue;
            }

            let record: MergedIsoform = archive_cache.read_bytes(ext_match);

            tx_abundance.add_fsm_misoform(&record);
        }

        let exact_matches_set = exact_matches
            .iter()
            .collect::<std::collections::HashSet<_>>();

        let mut candidates_offsets_set = HashSet::default();

        // only keep the misoform that no more than splice junctions as the query
        for offsets in position_based_matches.iter() {
            for offset in offsets {
                if exact_matches_set.contains(offset) {
                    continue;
                };

                if offset.n_splice_sites >= 2 && offset.n_splice_sites < queries.len() as u32 {
                    candidates_offsets_set.insert(offset);
                }
            }
        }

        let sj_pairs = trans.get_splice_junction_pairs();

        for cand_offset in candidates_offsets_set.iter() {
            // load the isoform from archive
            // let record = archive_cache.read_bytes(cand_offset);

            // let (is_all_matched, first_match, matched_count) =
            //     record.match_splice_junctions(&sj_pairs, cli.flank);

            // if is_all_matched {
            //     let mut msjc = MSJC::new(curr_msjc_id, dataset_info.get_size(), &record);
            //     curr_msjc_id += 1;

            //     tx_abundance.add_msjc(&mut msjc);
            //     // MSJC_vec.push(msjc);
            // }
        }

        TxAbundance_vec.push(tx_abundance);
    }

    // info out how many tx_abundance and msjc created
    info!(
        "Created {} TxAbundance and {} MSJC records",
        curr_tx_id.to_formatted_string(&Locale::en),
        curr_msjc_id.to_formatted_string(&Locale::en)
    );

    info!("Start EM optimization");

    for msjc in MSJC_vec.iter_mut() {
        msjc.prepare_em();
    }

    info!("Finish preparing MSJC records");

    let em_iter = 50; // temporary set the end point...

    for iter in 0..em_iter {
        info!("EM iteration {}", iter + 1);

        for msjc in MSJC_vec.iter_mut() {
            msjc.e_step(&TxAbundance_vec);
        }

        for tx_abundance in TxAbundance_vec.iter_mut() {
            tx_abundance.m_step(&MSJC_vec);
        }

        // print the frist 5 tx_abundance for debug
    }

    Ok(())
}
