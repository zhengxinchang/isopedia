use std::path::PathBuf;
use std::sync::mpsc::sync_channel;

use crate::{
    // assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    global_stats::GlobalStats,
    grouped_tx::{ChromGroupedTxManager, TmpOutputManager},
    gtf::{open_gtf_reader, TranscriptChunker},
    meta::Meta,
    myio::{DBInfos, Header},
    ptir_archive::PTIRArchiveCache,
    results::TableOutput,
    utils::{greetings2, log_wall_time, start_wall_timer},
};
use anyhow::Result;
use clap::Parser;
use log::{error, info, warn};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;

use rayon::{prelude::*, ThreadPoolBuilder};
use sysinfo::{Pid, ProcessRefreshKind, System};

pub(crate) const ISOFORM_FORMAT: &str = "RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_FSM_JC_EXACT:RC_FSM_JC_WOBBLE_ONLY:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES";

#[derive(Parser, Debug, Serialize, Clone)]
#[command(name = "isopedia isoform")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[Query] Annotate provided gtf file(transcripts/isoforms) with the index.
", long_about = None)]
#[clap(after_long_help = "

# Isopedia isoform needs the input gtf file to be sorted. use the following command to sort the gtf file:
gffread -T -o- input.gtf  | sort -k1,1 -k4,4n | gffread - -o sorted.gtf

Each sample field is:
RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_FSM_JC_EXACT:RC_FSM_JC_WOBBLE_ONLY:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES
All CPM fields use the sample's total indexed evidence as the denominator.
Direct detection uses RC_FSM_JC only; EM-assigned reads contribute to abundance estimates, not detection.
For multi-exon transcripts, --flank controls candidate search and full-junction-chain matching; exact means every junction coordinate equals the annotation. For mono-exon transcripts, --mono-exon-wobble controls both candidate search and two-end FSM matching; exact/wobble-only describe equality/tolerance of both ends, not a junction chain. --tss-wob and --tes-wob independently score terminal support in either direction.

")]
pub struct AnnIsoCli {
    /// Path to the index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// Path to the GTF file
    #[arg(short, long)]
    pub gtf: PathBuf,

    /// Maximum deviation (bp) for multi-exon junction candidate search and full-chain matching; 0 requires exact coordinates
    #[arg(short, long, default_value_t = 10)]
    pub flank: u64,

    /// Expand mono-exon searches by this many bp on each side; FSM requires both read ends within this distance of the annotated ends
    #[arg(short = 'F', long, default_value_t = 50)]
    pub mono_exon_wobble: u64,

    /// Minimum number of reads required to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// Output file for search results
    #[arg(short, long)]
    pub output: PathBuf,

    /// Whether to include additional information in the output
    #[arg(long, default_value_t = false)]
    pub info: bool,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    pub num_threads: usize,

    /// Max EM iterations
    #[arg(long, default_value_t = 100)]
    pub em_max_iter: usize,

    /// EM convergence threshold
    #[arg(long, default_value_t = 0.01)]
    pub em_conv_min_diff: f32,

    /// EM chunk size, reduce it if you have low memory
    #[arg(long, default_value_t = 4)]
    pub em_chunk_size: usize,

    /// EM effective length coefficient, avoid divide by zero when transcript is very short.
    #[arg(long, default_value_t = 2)]
    pub em_effective_len_coef: usize,

    /// EM damping factor
    #[arg(long, default_value_t = 0.3)]
    pub em_damping_factor: f32,

    /// Minimum EM abundance to report
    #[arg(long, default_value_t = 0.0001)]
    pub min_em_abundance: f32,

    /// Deprecated compatibility flag; TSS/TES support is always counted
    #[arg(long, default_value_t = false)]
    pub no_check_tss_tes: bool,

    /// Maximum absolute deviation between read and annotated TSS positions
    #[arg(
        long,
        default_value_t = 50,
        help = "Maximum allowed absolute deviation (bp) between read and annotated TSS positions. The tolerance applies in both directions."
    )]
    pub tss_wob: u64,

    /// Maximum absolute deviation between read and annotated TES positions
    #[arg(
        long,
        default_value_t = 50,
        help = "Maximum allowed absolute deviation (bp) between read and annotated TES positions. The tolerance applies in both directions."
    )]
    pub tes_wob: u64,

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

        if self.min_read == 0 {
            error!("--min-read must be at least 1");
            is_ok = false;
        }

        if is_ok != true {
            std::process::exit(1);
        }
    }
}

pub fn run_isoform_annotation(cli: &AnnIsoCli) -> Result<()> {
    let started = start_wall_timer();
    greetings2(&cli);
    cli.validate();
    if cli.no_check_tss_tes {
        warn!("--no-check-tss-tes is obsolete: TSS/TES support is always counted for the isoform output");
    }

    ThreadPoolBuilder::new()
        .num_threads(cli.num_threads)
        .build_global()
        .expect("Can not allocate thread pool");

    let forest = BPForest::init(&cli.idxdir);
    let index_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    info!("Loading GTF file...");

    let gtfreader = open_gtf_reader(cli.gtf.to_str().unwrap())?;

    let mut gtf = TranscriptChunker::new(gtfreader);

    let gtf_by_chrom = gtf.get_all_transcripts_by_chrom();

    let missing_chroms: Vec<&String> = gtf_by_chrom
        .iter()
        .map(|(chrom, _)| chrom)
        .filter(|chrom| forest.chrom_mapping.get_chrom_idx(chrom).is_none())
        .collect();

    if !missing_chroms.is_empty() {
        warn!(
            "The following query chromosomes are not present in the index and will return no matches: {}",
            missing_chroms
                .iter()
                .map(|chrom| chrom.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );

        if missing_chroms.len() == gtf_by_chrom.len() {
            warn!(
                "None of the query chromosomes are present in the index; all annotation rates will be 0."
            );
        }
    }

    info!(
        "Loaded {} transcripts from gtf file",
        gtf.trans_count.to_formatted_string(&Locale::en)
    );

    info!("Loading index file");

    drop(forest);

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
    out_header.add_column("ranking_score")?;
    out_header.add_column("detected")?;
    out_header.add_column("min_read")?;
    out_header.add_column("n_pos_samples/sample_size")?;
    out_header.add_column("attributes")?;
    let mut db_infos = DBInfos::new();
    for (name, evidence) in index_info.get_sample_evidence_pair_vec() {
        db_infos.add_sample_evidence(&name, evidence);
        out_header.add_sample_name(&name)?;
    }

    info!("Initializing transcript groups on demand");
    info!("Processing transcripts");
    let tmp_path = cli.output.with_extension("tmp");
    let mut tmp_tx_manger = std::thread::scope(|scope| {
        let (sender, receiver) = sync_channel(2);
        let writer = scope.spawn(|| {
            let mut manager = TmpOutputManager::new(&tmp_path, cli);
            for chunk in receiver {
                for txs in chunk {
                    manager.dump_txs(txs);
                }
            }
            manager
        });

        gtf_by_chrom.into_par_iter().for_each_init(
            || {
                (
                    BPForest::init(&cli.idxdir),
                    PTIRArchiveCache::new(
                        cli.idxdir.join(MERGED_FILE_NAME),
                        cli.cached_chunk_size_mb * 1024 * 1024,
                        cli.cached_chunk_num,
                    ),
                    System::new(),
                )
            },
            |(forest, archive_cache, sys), (chrom, tx_vec)| {
                let pid = Pid::from_u32(std::process::id());
                sys.refresh_processes_specifics(
                    sysinfo::ProcessesToUpdate::Some(&[pid]),
                    true,
                    ProcessRefreshKind::everything(),
                );
                let mem_before = sys.process(pid).unwrap().memory() / 1024 / 1024;

                let mut chrom_manager = ChromGroupedTxManager::new(&chrom, index_info.get_size());
                chrom_manager.process_transcripts(
                    tx_vec,
                    forest,
                    cli,
                    archive_cache,
                    &index_info,
                    |chunk| {
                        sender.send(chunk).expect("temporary output writer stopped");
                    },
                );
                forest.clear_all_caches();
                archive_cache.clear_cache();
                chrom_manager.clear();

                sys.refresh_processes_specifics(
                    sysinfo::ProcessesToUpdate::Some(&[pid]),
                    true,
                    ProcessRefreshKind::everything(),
                );
                let mem_after = sys.process(pid).unwrap().memory() / 1024 / 1024;
                info!(
                    "Chromosome {}: memory {}MB -> {}MB (delta: {:+}MB)",
                    chrom,
                    mem_before,
                    mem_after,
                    mem_after as i64 - mem_before as i64
                );
            },
        );
        drop(sender);
        writer.join().expect("temporary output writer failed")
    });

    let mut tableout = TableOutput::new(
        cli.output.clone(),
        out_header,
        db_infos,
        meta,
        ISOFORM_FORMAT.to_string(),
    );

    info!("Finalizing temporary output");

    {
        tmp_tx_manger.finish();
        info!("Sorting final output table");
        let mut line = Vec::new();
        while let Some(tx_abd_view) = tmp_tx_manger.next() {
            // info!("writeing transcript {}", tx_abd.orig_tx_id);
            // let mut line = tx_abd.to_output_line(&global_stats, &dataset_info, &cli);
            tx_abd_view.write_line_directly(
                &mut global_stats,
                &index_info,
                cli,
                &mut line,
                &mut tableout,
            )?;
            // tableout.add_line(&mut line)?;
            // tx_abd.write_line_directly( &global_stats, &dataset_info, &cli,&mut tableout)?;
            // info!("writen transcript {} done", tx_abd.orig_tx_id);
        }
    }
    info!("Writing final output table");

    info!("Cleaning up temporary files");
    tmp_tx_manger.clean_up()?;

    tableout.finish()?;

    info!("> Stats summary:");
    info!("> Sample\tDirect-supported transcripts (pct)");
    for sample_name in index_info.get_sample_names() {
        let sample_idx = index_info
            .get_sample_idx_by_name(&sample_name)
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let fsm_count = global_stats.get_fsm_tx_by_sample_idx(sample_idx);
        let total_tx = gtf.trans_count as f32;
        let fsm_pct = if total_tx > 0.0 {
            (fsm_count as f32 / total_tx) * 100.0
        } else {
            0.0
        };

        info!("> {}\t{}({:.2}%)", sample_name, fsm_count, fsm_pct);
    }
    info!("Save output to file {:?}", tableout.get_out_path().unwrap());

    info!("Finished!");
    log_wall_time("isoform", started);
    Ok(())
}
