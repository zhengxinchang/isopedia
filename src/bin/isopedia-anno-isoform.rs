use std::{
    env,
    io::{BufReader, BufWriter},
    path::PathBuf,
    vec,
    fs::File,
};

use anyhow::Result;
use clap::{command, Parser};
use isopedia::{
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    gtf::TranscriptChunker,
    isoform::{self, MergedIsoform},
    isoformarchive::read_record_from_mmap,
    meta::Meta,
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, get_total_memory_bytes, warmup},
};
use log::{error, info};
use memmap2::Mmap;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-ann-isoform")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
#[clap(after_long_help = "
")]
struct Cli {
    /// index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// gtf file
    #[arg(short, long)]
    pub gtf: PathBuf,

    /// flank size for search, before and after the position
    #[arg(short, long, default_value_t = 10)]
    pub flank: u64,

    /// minimal reads to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// output file for search results
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// memory size for warming up, in Gigabytes, example: 1GB
    #[arg(short, long, default_value_t = 4)]
    pub warmup_mem: usize,

    /// number of cached nodes for each tree in maximal
    #[arg(short, long, default_value_t = 1024)]
    pub max_cached_nodes: usize,

}

impl Cli {
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
                "--idxdir: Aggr file {} does not exist in {}, please run `stix-isoform aggr` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join("bptree_0.idx").exists() {
            error!(
                "--idxdir: index files {} does not exist in {}, please run `stix-isoform idx` first",
                MERGED_FILE_NAME,
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
            if warmup_bytes > max_mem_bytes  {
                error!(
                    "--warmup-mem: {} GB larger than system memory, please set it to less than {} GB",
                    self.warmup_mem,
                    max_mem_bytes / (1024 * 1024 * 1024)
                );
                is_ok = false;
            }
        }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &Cli) {
    println!("\nIsopedia: [Annotate provided gtf file]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => println!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() -> Result<()> {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;
    let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer

    info!("Search by gtf/gff file");
    let gtfreader: noodles_gtf::Reader<BufReader<std::fs::File>> = noodles_gtf::io::Reader::new(
        BufReader::new(std::fs::File::open(cli.gtf).expect("can not read gtf")),
    );

    let gtf = TranscriptChunker::new(gtfreader);

    let mut hit_count = 0u32;
    let mut miss_count = 0u32;

    //output file,if cli.output is none then write to stdout otherwise write to file

    let mut writer = match cli.output {
        Some(ref output) => {
            let f = std::fs::File::create(output).expect("can not create output file");
            Box::new(BufWriter::new(f)) as Box<dyn std::io::Write>
        }
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn std::io::Write>,
    };

    let meta = Meta::parse(&cli.idxdir.join(META_FILE_NAME))?;

    info!("Writing sample meta...");
    writer.write(meta.get_meta_table(Some("##")).as_bytes())?;

    writer.write(
        "#chrom\tstart\tend\tlength\texon_count\ttrans_id\tgene_id\tconfidence\tdetected\tmin_read\tpositive_count/sample_size\tattributes\tFORMAT".as_bytes()
    )?;

    dataset_info.get_sample_names().iter().for_each(|x| {
        writer
            .write(format!("\t{}", x).as_bytes())
            .expect("Failed to write sample header");
    });

    writer.write("\n".as_bytes())?;


    info!("Warmup index file");
    let max_gb = cli.warmup_mem * 1024 * 1024 * 1024;
    warmup(&cli.idxdir.clone().join(MERGED_FILE_NAME), max_gb)?;


    info!("Loading index file");
    let archive_file_handle = File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
        .expect("Can not open aggregated records file...");



    let archive_mmap = unsafe { Mmap::map(&archive_file_handle).expect("Failed to map the file") };

    archive_mmap.advise(memmap2::Advice::Sequential).expect("Failed to set mmap advice");

    info!("Start to process transcripts");
    let mut iter_count = 0;
    let mut batch = 0;
    let mut acc_pos_count = vec![0u32; dataset_info.get_size()];
    let mut acc_sample_evidence_arr = vec![0u32; dataset_info.get_size()];

    let mut total_acc_evidence_flag_vec = vec![0u32; dataset_info.get_size()];

    const FORMAT: &str = "CPM:COUNT";
    for trans in gtf {
        iter_count += 1;

        if iter_count == 10_000 {
            batch += 1;
            info!("Processed {} transcripts", (batch * 10_000).to_formatted_string(&Locale::en));
            iter_count = 0;
        }

        let mut queries: Vec<(String, u64)> = trans.get_quieries();
        queries.sort_by_key(|x| x.1);
        let res = forest.search_all_match(&queries, cli.flank);

        if res.is_none() {
            // error!("No results found for queries: {:?}", queries);
            continue;
        }
        let res = res.unwrap();

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
                &dataset_info,
            );

            write!(
                writer,
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
                FORMAT
            )?;

            for (i, val) in acc_sample_evidence_arr.iter().enumerate() {
                if i > 0 {
                    write!(writer, "\t")?;
                }

                let cpm = utils::calc_cpm(val, &dataset_info.sample_total_evidence_vec[i]);

                write!(writer, "{}:{}", cpm, val)?;
            }
            write!(writer, "\n")?;
        } else {
            write!(
                writer,
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
                FORMAT
            )?;
            let sample_count = dataset_info.get_sample_names().len();
            for i in 0..sample_count {
                if i > 0 {
                    write!(writer, "\t")?;
                }
                write!(writer, "0:0")?;
            }

            write!(writer, "\n")?;
        }
    }
    let total = hit_count + miss_count;
    info!("Processed {} transcripts", (100000 * batch + iter_count).to_formatted_string(&Locale::en));
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
