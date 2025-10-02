use std::{env, fs::File, io::BufReader, path::PathBuf, vec};

use crate::{
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    gtf::TranscriptChunker,
    io::{DBInfos, GeneralOutputIO, Header, Line, SampleChip},
    isoform::{self, MergedIsoform},
    isoformarchive::read_record_from_mmap,
    meta::Meta,
    output::{GeneralTableOutput, IsoformTableOut},
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, get_total_memory_bytes, warmup},
    writer::MyGzWriter,
};
use anyhow::Result;
use clap::{command, Parser};
use log::{error, info};
use memmap2::Mmap;
use nix::libc::{self, posix_madvise, POSIX_MADV_DONTNEED};
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
// use nix::sys::mman::{posix_madvise, PosixMadvise};
// Removed incorrect import of mmap::{PosixMadvise,posix_madvise}

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-ann-isoform")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
[Annotation] Annotate provided gtf file(transcripts/isoforms) with the index.
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

    /// Memory size to use for warming up (in gigabytes).  
    /// Example: 1GB. Increasing this will significantly improve performance;  
    /// set it as large as your system allows.
    #[arg(short, long, default_value_t = 4)]
    pub warmup_mem: usize,

    /// Maximum number of cached nodes per tree
    #[arg(short = 'c', long = "cached_nodes", default_value_t = 100_000)]
    pub lru_size: usize,
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
                "--idxdir: Aggr file {} does not exist in {}, please run `isopedia-aggr` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join("bptree_0.idx").exists() {
            error!(
                "--idxdir: bptree_0.idx file does not exist in {}, please run `isopedia-idx` first",
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.gtf.exists() {
            error!("--gtf: gtf file {} does not exist", self.gtf.display());
            is_ok = false;
        }

        let max_mem_bytes = match get_total_memory_bytes() {
            Some(bytes) => bytes,
            None => {
                error!("Failed to get total memory bytes");
                std::process::exit(1);
            }
        };

        if self.warmup_mem == 0 {
            error!("--warmup-mem: must be greater than 0");
            is_ok = false;
        } else {
            let warmup_bytes = (self.warmup_mem as u64) * 1024 * 1024 * 1024;
            if warmup_bytes > max_mem_bytes {
                error!(
                    "--warmup-mem: {} GB larger than system memory, please set it to less than {} GB",
                    self.warmup_mem,
                    max_mem_bytes / (1024 * 1024 * 1024)
                );
                is_ok = false;
            }
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
    env::set_var("RUST_LOG", "info");
    // env_logger::init();

    cli.validate();
    greetings(&cli);

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;
    let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer

    info!("Loading GTF file...");
    let gtfreader: noodles_gtf::Reader<BufReader<std::fs::File>> = noodles_gtf::io::Reader::new(
        BufReader::new(std::fs::File::open(cli.gtf.clone()).expect("can not read gtf")),
    );

    let gtf = TranscriptChunker::new(gtfreader);

    let mut hit_count = 0u32;
    let mut miss_count = 0u32;
    info!("Loading sample meta...");

    let meta = Meta::parse(&cli.idxdir.join(META_FILE_NAME), None)?;

    const ISOFORM_FORMAT: &str = "CPM:COUNT";

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

    let mut dbInfos = DBInfos::new();

    for (name, evidence) in dataset_info.get_sample_evidence_pair_vec() {
        dbInfos.add_sample_evidence(&name, evidence);
        out_header.add_sample_name(&name)?;
    }

    let mut isoform_out = IsoformTableOut::new(
        out_header,
        dbInfos,
        meta.clone(),
        ISOFORM_FORMAT.to_string(),
    );

    info!("Warmup index file");
    let max_gb = cli.warmup_mem * 1024 * 1024 * 1024;
    warmup(&cli.idxdir.clone().join(MERGED_FILE_NAME), max_gb)?;

    info!("Loading index file");
    let archive_file_handle = File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
        .expect("Can not open aggregated records file...");

    let archive_mmap = unsafe { Mmap::map(&archive_file_handle).expect("Failed to map the file") };

    archive_mmap
        .advise(memmap2::Advice::Sequential)
        .expect("Failed to set mmap advice");

    info!("Processing transcripts");
    let mut iter_count = 0;
    let mut batch = 0;
    let mut acc_pos_count = vec![0u32; dataset_info.get_size()];
    let mut acc_sample_evidence_arr = vec![0u32; dataset_info.get_size()];

    let mut total_acc_evidence_flag_vec = vec![0u32; dataset_info.get_size()];

    let mut released_bytes: usize = 0; // 已释放的区间长度

    const RELEASE_STEP: usize = 256 * 1024 * 1024; // 每次释放 256MB

    for trans in gtf {
        iter_count += 1;

        if iter_count == 10_000 {
            batch += 1;
            info!(
                "Processed {} transcripts",
                (batch * 10_000).to_formatted_string(&Locale::en)
            );
            iter_count = 0;

            // let release_up_to = released_bytes + RELEASE_STEP;
            // let base_ptr = archive_mmap.as_ptr() as *mut libc::c_void;

            // let ret = unsafe {
            //     posix_madvise(
            //         base_ptr.add(released_bytes),
            //         RELEASE_STEP,
            //         POSIX_MADV_DONTNEED,
            //     )
            // };

            // if ret == 0 {
            //     released_bytes = release_up_to;
            //     info!(
            //         "Released first {} MB from page cache",
            //         released_bytes / 1024 / 1024
            //     );
            // } else {
            //     info!(
            //         "posix_madvise failed at {} MB",
            //         released_bytes / 1024 / 1024
            //     );
            // }
        }

        let mut queries: Vec<(String, u64)> = trans.get_quieries();
        queries.sort_by_key(|x| x.1);
        let res = forest.search_all_match(&queries, cli.flank, cli.lru_size);

        // make sure the returned isoform has exactly the same number of splice sites as the query
        let target: Vec<MergedIsoformOffsetPtr> = res
            .into_iter()
            .filter(|x| x.n_splice_sites == queries.len() as u32)
            .collect();

        if target.len() == 0 {
            miss_count += 1;
        } else {
            hit_count += 1;
        }

        let mut out_line = Line::new();
        out_line.update_format_str(ISOFORM_FORMAT);

        if target.len() > 0 {
            // have hits
            acc_pos_count.fill(0);
            acc_sample_evidence_arr.fill(0);

            for offset in &target {
                let record: MergedIsoform =
                    read_record_from_mmap(&archive_mmap, offset, &mut archive_buf);

                acc_pos_count
                    .iter_mut()
                    .zip(record.get_positive_array(&cli.min_read))
                    .for_each(|(a, b)| *a += b);

                acc_sample_evidence_arr
                    .iter_mut()
                    .zip(record.get_sample_evidence_arr())
                    .for_each(|(a, b)| *a += b);
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
                acc_sample_evidence_arr.clone(),
                dataset_info.get_size(),
                &dataset_info.sample_total_evidence_vec,
            );

            let mut outline = format!(
                "{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}\t{}\tyes\t{}\t{}/{}\t{}\t{}\t",
                trans.chrom,
                trans.start,
                trans.end,
                trans.get_transcript_length(),
                trans.get_exon_count(),
                trans.trans_id,
                trans.gene_id,
                confidence,
                &cli.min_read,
                acc_pos_count.iter().filter(|&&x| x > 0).count(),
                dataset_info.get_size(),
                trans.get_attributes(),
                ISOFORM_FORMAT
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

            for (i, val) in acc_sample_evidence_arr.iter().enumerate() {
                if i > 0 {
                    // write!(writer, "\t")?;
                    outline.push('\t');
                }

                let cpm = utils::calc_cpm(val, &dataset_info.sample_total_evidence_vec[i]);

                let samplechip = SampleChip::new(
                    Some(dataset_info.get_sample_names()[i].clone()),
                    vec![format!("{}:{}", cpm, val)],
                );
                out_line.add_sample(samplechip);
                // write!(writer, "{}:{}", cpm, val)?;
                outline.push_str(&format!("{}:{}", cpm, val));
            }
            // write!(writer, "\n")?;
            isoform_out.add_line(&out_line)?;
            outline.push('\n');
            // mywriter.write_all_bytes(outline.as_bytes())?;
        } else {
            let mut output = format!(
                "{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}\t0\tno\t{}\t0/{}\t{}\t{}\t",
                trans.chrom,
                trans.start,
                trans.end,
                trans.get_transcript_length(),
                trans.get_exon_count(),
                trans.trans_id,
                trans.gene_id,
                &cli.min_read,
                dataset_info.get_size(),
                trans.get_attributes(),
                ISOFORM_FORMAT
            );

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
                if i > 0 {
                    output.push('\t');
                }

                let samplechip = SampleChip::new(
                    Some(dataset_info.get_sample_names()[i].clone()),
                    vec!["0:0".to_string()],
                );
                out_line.add_sample(samplechip);

                output.push_str("0:0");
            }
            isoform_out.add_line(&out_line)?;

            output.push('\n');
            // mywriter.write_all_bytes(output.as_bytes())?;
        }
    }

    isoform_out.save_to_file(&cli.output.with_extension("isoform.gz"))?;

    let total = hit_count + miss_count;
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
        hit_count.to_formatted_string(&Locale::en),
        miss_count.to_formatted_string(&Locale::en),
        (hit_count + miss_count).to_formatted_string(&Locale::en),
        hit_count as f64 / (hit_count + miss_count) as f64 * 100f64
    );

    info!("Finished");
    Ok(())
}
