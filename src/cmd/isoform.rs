use std::{path::PathBuf, vec};

use crate::{
    assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
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

    let mut assembler = Assembler::init(dataset_info.get_size());

    let mut runtime = Runtime::init(dataset_info.get_size());

    for trans in gtf_vec {
        iter_count += 1;

        runtime.add_one_total();
        runtime.reset();

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

        let mut queries: Vec<(String, u64)> = trans.get_quieries();
        queries.sort_by_key(|x| x.1);

        let (hits, all_res) = forest.search2_all_match(&queries, cli.flank, cli.lru_size);

        // enable the assembler to assemble fragmented reads into isoforms
        let ism_hit = if cli.asm {
            assembler.reset();
            let is_hit =
                assembler.assemble(&queries, &all_res, &mut archive_cache, cli, &mut runtime);

            if cli.debug {
                // assembler.print_matrix();
            }
            is_hit
        } else {
            false
        };

        // make sure the returned isoform has exactly the same number of splice sites as the query
        let target: Vec<MergedIsoformOffsetPtr> = hits
            .into_iter()
            .filter(|x| x.n_splice_sites == queries.len() as u32)
            .collect();

        if target.len() > 0 {
            runtime.add_one_fsm_hit();
        } else {
            if ism_hit {
                runtime.add_one_ism_hit();
            }
        }

        let mut out_line = Line::new();
        out_line.update_format_str(ISOFORM_FORMAT);

        if target.len() > 0 {
            if cli.debug {
                info!(
                    "Transcript {}:{}-{} has {} matched isoforms in the index.",
                    trans.chrom,
                    trans.start,
                    trans.end,
                    target.len()
                );
            }

            for offset in &target {
                let record: MergedIsoform = archive_cache.read_bytes(offset);

                runtime.update_fsm_record(&record)?;

                if cli.info {
                    runtime.fsm_trigger_add_read_info(&record)?;
                }
            }

            out_line.add_field(&trans.chrom);
            out_line.add_field(&trans.start.to_string());
            out_line.add_field(&trans.end.to_string());
            out_line.add_field(&trans.get_transcript_length().to_string());
            out_line.add_field(&trans.get_exon_count().to_string());
            out_line.add_field(&trans.trans_id);
            out_line.add_field(&trans.gene_id);
            // out_line.add_field(&trans.get_attributes());

            let confidence = match cli.asm {
                true => runtime.get_fsm_ism_confidence(&dataset_info),
                false => runtime.get_fsm_confidence(&dataset_info),
            };
            out_line.add_field(&confidence.to_string());

            out_line.add_field("yes");
            out_line.add_field(&cli.min_read.to_string());

            let positive_count = match cli.asm {
                true => runtime.get_fsm_ism_positive_count_by_min_read(cli),
                false => runtime.get_fsm_positive_count_by_min_read(&cli),
            };
            out_line.add_field(&format!("{}/{}", positive_count, dataset_info.get_size()));

            out_line.add_field(&trans.get_attributes());

            for (idx, data) in runtime
                .get_sample_chip_data(&dataset_info, cli)
                .iter()
                .enumerate()
            {
                let samplechip = SampleChip::new(
                    Some(dataset_info.get_sample_names()[idx].clone()),
                    vec![format!("{}:{}:{}", data.0, data.1, data.2)], // cpm:count:info
                );
                out_line.add_sample(samplechip);
            }
            tableout.add_line(&out_line)?;
        } else {
            out_line.add_field(&trans.chrom);
            out_line.add_field(&trans.start.to_string());
            out_line.add_field(&trans.end.to_string());
            out_line.add_field(&trans.get_transcript_length().to_string());
            out_line.add_field(&trans.get_exon_count().to_string());
            out_line.add_field(&trans.trans_id);
            out_line.add_field(&trans.gene_id);
            // out_line.add_field(&trans.get_attributes());
            out_line.add_field("0");
            out_line.add_field("no");
            out_line.add_field(&cli.min_read.to_string());
            out_line.add_field(&format!("0/{}", dataset_info.get_size()));
            out_line.add_field(&trans.get_attributes());
            // for _ in 0..dataset_info.get_size() {

            let sample_count = dataset_info.get_sample_names().len();
            for i in 0..sample_count {
                let samplechip = SampleChip::new(
                    Some(dataset_info.get_sample_names()[i].clone()),
                    vec!["0:0:NULL".to_string()],
                );
                out_line.add_sample(samplechip);
            }
            // isoform_out.add_line(&out_line)?;
            tableout.add_line(&out_line)?;
        }

        runtime.log_stats();
    }

    // log_mem_stats();

    let total = runtime.total;
    let missed = total - runtime.fsm_hit - runtime.ism_hit;
    info!(
        "Processed {} transcripts",
        (10_000 * batch + iter_count).to_formatted_string(&Locale::en)
    );
    info!("Sample-wide stats: ");
    info!("> Sample\t\ttotal\thit(fsm)\thit(ism)\tmiss\tpct");
    for i in 0..dataset_info.get_size() {
        info!(
            "> {:}\t{}\t{}\t{}\t{}\t{:.2}%",
            dataset_info.get_sample_names()[i],
            total,
            runtime.sample_wide_positive_transcript_fsm[i],
            runtime.sample_wide_positive_transcript_ism[i],
            total - runtime.sample_wide_positive_transcript_both[i],
            (runtime.sample_wide_positive_transcript_both[i]) as f64 / total as f64 * 100f64
        );
    }
    info!(
        "Index-wide stats: total: {}, hit(fsm): {}, hit(ism): {}, miss: {},  pct: {:.2}%",
        (runtime.total).to_formatted_string(&Locale::en),
        runtime.fsm_hit.to_formatted_string(&Locale::en),
        runtime.ism_hit.to_formatted_string(&Locale::en),
        missed.to_formatted_string(&Locale::en),
        (runtime.fsm_hit + runtime.ism_hit) as f64 / (runtime.total) as f64 * 100f64
    );

    tableout.finish()?;
    info!("Finished");
    Ok(())
}
