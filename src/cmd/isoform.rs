use std::path::PathBuf;

use crate::{
    // assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    global_stats::GlobalStats,
    grouped_tx::{ChromGroupedTxManager, TmpOutputManager},
    gtf::{open_gtf_reader, TranscriptChunker},
    isoformarchive::ArchiveCache,
    meta::Meta,
    myio::{DBInfos, Header},
    results::TableOutput,
    utils::greetings2,
};
use anyhow::Result;
use clap::Parser;
use log::{error, info};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;

use rayon::ThreadPoolBuilder;
use sysinfo::{Pid, ProcessRefreshKind, System};

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

    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    pub num_threads: usize,

    /// Max EM iterations
    #[arg(long, default_value_t = 50)]
    pub em_max_iter: usize,

    // EM convergence threshold
    #[arg(long, default_value_t = 0.01)]
    pub em_conv_min_diff: f32,

    // EM chunk size, reduce it if you have low memory
    #[arg(long, default_value_t = 4)]
    pub em_chunk_size: usize,

    // EM effective length coefficient
    #[arg(long, default_value_t = 1)]
    pub em_effective_len_coef: usize,

    // EM damping factor
    #[arg(long, default_value_t = 0.3)]
    pub em_damping_factor: f32,

    // Minimum EM abundance to report
    #[arg(long, default_value_t = 0.01)]
    pub min_em_abundance: f32,

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

    /// Output temporary shard record size
    #[arg(long, default_value_t = 10_000)]
    pub output_tmp_shard_counts: usize,
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

pub fn run_anno_isoform(cli: &AnnIsoCli) -> Result<()> {
    greetings2(&cli);
    cli.validate();

    ThreadPoolBuilder::new()
        .num_threads(cli.num_threads)
        .build_global()
        .expect("Can not allocate thread pool");

    let mut forest = BPForest::init(&cli.idxdir);
    let index_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    info!("Loading GTF file...");

    let gtfreader = open_gtf_reader(cli.gtf.to_str().unwrap())?;

    let mut gtf = TranscriptChunker::new(gtfreader);

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

    let mut global_stats = GlobalStats::new(index_info.get_size());

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
    for (name, evidence) in index_info.get_sample_evidence_pair_vec() {
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
            let mut manager = ChromGroupedTxManager::new(&chrom, index_info.get_size());
            manager.add_transcript_by_chrom(&tx_vec, cli);
            manager
        })
        .collect();

    let mut tmp_tx_manger = TmpOutputManager::new(&cli.output.with_extension(&"tmp"), cli);

    let mut sys = System::new_all();

    let pid = Pid::from_u32(std::process::id());

    info!("Processing transcripts");
    for chrom_manager in chrom_grouped_tx_managers.iter_mut() {
        sys.refresh_processes_specifics(
            sysinfo::ProcessesToUpdate::Some(&[pid]),
            true,
            ProcessRefreshKind::everything(),
        );
        let process = sys.process(pid).unwrap();
        let mem_before = process.memory() / 1024 / 1024; // MB

        chrom_manager.process_tx_groups(
            &mut forest,
            cli,
            &mut archive_cache,
            &index_info,
            &mut tmp_tx_manger,
            &mut global_stats,
        );

        forest.clear_all_caches();

        archive_cache.clear_cache();

        chrom_manager.clear();

        sys.refresh_processes_specifics(
            sysinfo::ProcessesToUpdate::Some(&[pid]),
            true,
            ProcessRefreshKind::everything(),
        );
        let process = sys.process(pid).unwrap();
        let mem_after = process.memory() / 1024 / 1024; // MB

        info!(
            "Chromosome {}: memory {}MB -> {}MB (delta: {:+}MB)",
            chrom_manager.chrom,
            mem_before,
            mem_after,
            mem_after as i64 - mem_before as i64
        );
    }

    info!("Finalizing temporary output");

    {
        tmp_tx_manger.finish();
        info!("Sorting final output table");
        while let Some(tx_abd_view) = tmp_tx_manger.next() {
            // info!("writeing transcript {}", tx_abd.orig_tx_id);
            // let mut line = tx_abd.to_output_line(&global_stats, &dataset_info, &cli);
            tx_abd_view.write_line_directly(&mut global_stats, &index_info, cli, &mut tableout)?;
            // tableout.add_line(&mut line)?;
            // tx_abd.write_line_directly( &global_stats, &dataset_info, &cli,&mut tableout)?;
            // info!("writen transcript {} done", tx_abd.orig_tx_id);
        }
    }
    info!("Writing final output table");

    info!("Cleaning up temporary files");
    tmp_tx_manger.clean_up()?;

    tableout.finish()?;

    // output global stats each sample is a row, columns are  sample name, fsm_total(pct), em_total(pct), fsm_em_total(pct)
    info!("> Stats summary:");
    info!("> Sample FSM(CPT) EM(CPT) FSM+EM(CPT)");
    for sample_name in index_info.get_sample_names() {
        let sample_idx = index_info
            .get_sample_idx_by_name(&sample_name)
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let fsm_count = global_stats.sample_posi_tx_count_fsm[sample_idx];
        let em_count = global_stats.sample_posi_tx_count_em[sample_idx];
        let fsm_em_count = global_stats.sample_posi_tx_count_fsm_em[sample_idx];
        let total_tx = gtf.trans_count as f32;
        let fsm_pct = (fsm_count as f32 / total_tx) * 100.0;
        let em_pct = (em_count as f32 / total_tx) * 100.0;
        let fsm_em_pct = (fsm_em_count as f32 / total_tx) * 100.0;

        info!(
            "> {}\t{}({:.2}%)\t{}({:.2}%)\t{}({:.2}%)",
            sample_name, fsm_count, fsm_pct, em_count, em_pct, fsm_em_count, fsm_em_pct
        );
    }

    info!("Finished!");
    Ok(())
}
