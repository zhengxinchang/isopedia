use std::{fs::File, path::PathBuf, vec};

use crate::{
    assemble::Assembler,
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    gtf::{open_gtf_reader, TranscriptChunker},
    io::{DBInfos, Header, Line, SampleChip},
    isoform::{self, MergedIsoform},
    isoformarchive::ArchiveCache,
    meta::Meta,
    // output_traits::GeneralTableOutputTrait,
    results::TableOutput,
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, log_mem_stats},
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
#[command(version = "0.1.0")]
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
    pub use_incomplete: bool,

    /// Maximum number of cached tree nodes in memory
    #[arg(long = "cached_nodes", default_value_t = 1000)]
    pub lru_size: usize,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_number: usize,

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
    eprintln!("\nIsopedia: [Annotate provided gtf file]\n");
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
    info!("Loaded {} transcripts from gtf file", gtf_vec.len());

    let mut fsm_hit_count = 0u32;
    let mut ism_hit_count = 0u32;
    let mut miss_count = 0u32;
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

    // let mut isoform_out = IsoformTableOut::new(
    //     cli.output.clone().with_extension("old.gz"),
    //     out_header,
    //     db_infos,
    //     meta.clone(),
    //     ISOFORM_FORMAT.to_string(),
    // );

    // info!("Warmup index file");
    // let max_gb = cli.warmup_mem * 1024 * 1024 * 1024;
    // warmup(&cli.idxdir.clone().join(MERGED_FILE_NAME), max_gb)?;

    info!("Loading index file");

    let mut archive_cache = ArchiveCache::new(
        cli.idxdir.clone().join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb * 1024 * 1024, // 512MB chunk size
        cli.cached_chunk_number,                // max 4 chunks in cache ~2GB
    );

    info!("Processing transcripts");
    let mut iter_count = 0;
    let mut batch = 0;
    let mut acc_pos_count = vec![0u32; dataset_info.get_size()];
    let mut acc_sample_evidence_arr = vec![0u32; dataset_info.get_size()];
    let mut acc_sample_read_info_arr = vec!["".to_string(); dataset_info.get_size()];

    let mut total_acc_evidence_flag_vec = vec![0u32; dataset_info.get_size()];
    let mut assembler = Assembler::init(dataset_info.get_size());
    let mut assembled = (false, Vec::with_capacity(dataset_info.get_size()));

    for trans in gtf_vec {
        iter_count += 1;

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
                let acc_sample_read_info_arr_mem: usize =
                    acc_sample_read_info_arr.iter().map(|s| s.len()).sum();
                info!(
                    "Current output in-memory size: {} MB, tree in-memory size: {} MB, read info array size: {} MB",
                    isoform_out_mem / (1024 * 1024),
                    tree_mem / (1024 * 1024),
                    acc_sample_read_info_arr_mem / (1024 * 1024)
                );

                log_mem_stats();
            }
        }

        let mut queries: Vec<(String, u64)> = trans.get_quieries();
        queries.sort_by_key(|x| x.1);

        let (hits, all_res) = forest.search2_all_match(&queries, cli.flank, cli.lru_size);

        // make sure the returned isoform has exactly the same number of splice sites as the query
        let target: Vec<MergedIsoformOffsetPtr> = hits
            .into_iter()
            .filter(|x| x.n_splice_sites == queries.len() as u32)
            .collect();

        if target.len() == 0 {
            miss_count += 1;
        } else {
            fsm_hit_count += 1;
        }

        // enable the assembler to assemble fragmented reads into isoforms
        if cli.use_incomplete {
            // assembler.reset();
            // assembled = assembler.assemble(&queries, &all_res, &archive_mmap, cli.flank);

            // if assembled.0 {
            //     ism_hit_count += 1;
            // }
        }

        let mut out_line = Line::new();
        out_line.update_format_str(ISOFORM_FORMAT);

        if target.len() > 0 {
            // have hits
            acc_pos_count.fill(0);
            acc_sample_evidence_arr.fill(0);
            acc_sample_read_info_arr.fill("NULL".to_string());

            for offset in &target {
                // let record: MergedIsoform =
                //     read_record_from_mmap(&archive_mmap, offset, &mut archive_buf);

                let record: MergedIsoform = archive_cache.read_bytes(offset);

                acc_pos_count
                    .iter_mut()
                    .zip(record.get_positive_array(&cli.min_read))
                    .for_each(|(a, b)| *a += b);

                let single_sample_evidence_arr = record.get_sample_evidence_arr();

                acc_sample_evidence_arr
                    .iter_mut()
                    .zip(single_sample_evidence_arr.iter())
                    .for_each(|(a, b)| *a += b);

                if cli.info {
                    for (i, ofs) in record.get_sample_offset_arr().iter().enumerate() {
                        let ofs = *ofs as usize;
                        let length = single_sample_evidence_arr[i] as usize;
                        // get readdiffslim
                        if length > 0 {
                            let read_info_str = record.isoform_reads_slim_vec[ofs..ofs + length]
                                .iter()
                                .map(|delta| delta.to_string_no_offsets())
                                .collect::<Vec<String>>()
                                .join(",");
                            if acc_sample_read_info_arr[i] == "NULL" {
                                acc_sample_read_info_arr[i] = read_info_str;
                            } else {
                                acc_sample_read_info_arr[i] =
                                    format!("{},{}", acc_sample_read_info_arr[i], read_info_str);
                            }
                        }
                    }
                }
            }

            acc_sample_evidence_arr
                .iter()
                .enumerate()
                .for_each(|(i, x)| {
                    if *x > 0 {
                        total_acc_evidence_flag_vec[i] += 1;
                    }
                });

            let confidence = isoform::MergedIsoform::get_confidence_value(
                &acc_sample_evidence_arr,
                dataset_info.get_size(),
                &dataset_info.sample_total_evidence_vec,
            );

            out_line.add_field(&trans.chrom);
            out_line.add_field(&trans.start.to_string());
            out_line.add_field(&trans.end.to_string());
            out_line.add_field(&trans.get_transcript_length().to_string());
            out_line.add_field(&trans.get_exon_count().to_string());
            out_line.add_field(&trans.trans_id);
            out_line.add_field(&trans.gene_id);
            // out_line.add_field(&trans.get_attributes());
            out_line.add_field(&confidence.to_string());
            out_line.add_field("yes");
            out_line.add_field(&cli.min_read.to_string());
            out_line.add_field(&format!(
                "{}/{}",
                acc_pos_count.iter().filter(|&&x| x > 0).count(),
                dataset_info.get_size()
            ));
            out_line.add_field(&trans.get_attributes());

            for (i, acc_readc) in acc_sample_evidence_arr.iter().enumerate() {
                let cpm = utils::calc_cpm(acc_readc, &dataset_info.sample_total_evidence_vec[i]);

                let samplechip = if cli.info {
                    // only include in-compelte results when FSM mode is off
                    if cli.use_incomplete && *acc_readc == 0 {
                        let incomp_cpm = utils::calc_cpm(
                            &assembled.1[i],
                            &dataset_info.sample_total_evidence_vec[i],
                        );

                        SampleChip::new(
                            Some(dataset_info.get_sample_names()[i].clone()),
                            vec![format!(
                                "{}:{}:LEVEL=INCOMPLETE;{}",
                                incomp_cpm, assembled.1[i], ""
                            )],
                        )
                    } else {
                        SampleChip::new(
                            Some(dataset_info.get_sample_names()[i].clone()),
                            vec![format!(
                                "{}:{}:LEVEL=COMPLETE;{}",
                                cpm, acc_readc, acc_sample_read_info_arr[i]
                            )],
                        )
                    }
                } else {
                    if cli.use_incomplete && *acc_readc == 0 {
                        let incomp_cpm = utils::calc_cpm(
                            &assembled.1[i],
                            &dataset_info.sample_total_evidence_vec[i],
                        );

                        SampleChip::new(
                            Some(dataset_info.get_sample_names()[i].clone()),
                            vec![format!(
                                "{}:{}:LEVEL=INCOMPLETE",
                                incomp_cpm, assembled.1[i]
                            )],
                        )
                    } else {
                        SampleChip::new(
                            Some(dataset_info.get_sample_names()[i].clone()),
                            vec![format!("{}:{}:LEVEL=COMPLETE", cpm, acc_readc)],
                        )
                    }
                };

                out_line.add_sample(samplechip);
            }
            // isoform_out.add_line(&out_line)?;
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
    }

    // log_mem_stats();

    let total = fsm_hit_count + miss_count;
    info!(
        "Processed {} transcripts",
        (10_000 * batch + iter_count).to_formatted_string(&Locale::en)
    );
    info!("Sample-wide stats: ");
    info!("> Sample\thit\tmiss\tpct");
    for i in 0..dataset_info.get_size() {
        info!(
            "> {:}\t{}\t{}\t {:.2}%",
            dataset_info.get_sample_names()[i],
            total_acc_evidence_flag_vec[i],
            total - total_acc_evidence_flag_vec[i],
            total_acc_evidence_flag_vec[i] as f64 / total as f64 * 100f64
        );
    }
    info!(
        "Index-wide stats: hit: {}, miss: {}, total: {}, pct: {:.2}%",
        fsm_hit_count.to_formatted_string(&Locale::en),
        miss_count.to_formatted_string(&Locale::en),
        (fsm_hit_count + miss_count).to_formatted_string(&Locale::en),
        fsm_hit_count as f64 / (fsm_hit_count + miss_count) as f64 * 100f64
    );
    // info!("Writing output to {}", cli.output.display());
    // isoform_out.finish()?;
    tableout.finish()?;
    info!("Finished");
    Ok(())
}
