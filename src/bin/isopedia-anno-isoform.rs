use std::{
    env,
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
    vec,
};

use anyhow::Result;
use bincode::de;
use clap::{command, Parser};
use isopedia::{
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    gtf::TranscriptChunker,
    isoform::{self, MergedIsoform},
    isoformarchive::read_record_from_mmap,
    meta::Meta,
    output::OutputWriter,
    tmpidx::MergedIsoformOffsetPtr,
    utils::{self, get_total_memory_bytes, warmup},
    writer::MyGzWriter,
};
use log::{error, info};
use memmap2::Mmap;
use num_format::{Locale, ToFormattedString};
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-ann-isoform")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>, <xinchang.zheng@bcm.edu>
", long_about = None)]
#[clap(after_long_help = "
")]
struct Cli {
    /// Path to the index directories
    #[arg(short, long)]
    pub idxdirs: Option<Vec<PathBuf>>,

    /// Path to the manifest of a list of index directories, each line is a valid path to an index directory
    #[arg(short = 'I', long = "idx-manifest")]
    pub manifest: Option<PathBuf>,

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
    pub output: Option<PathBuf>,

    /// Memory size to use for warming up (in gigabytes).  
    /// Example: 1GB. Increasing this will significantly improve performance;  
    /// set it as large as your system allows.
    #[arg(short, long, default_value_t = 4)]
    pub warmup_mem: usize,

    /// Maximum number of cached nodes per tree
    #[arg(short = 'c', long = "cached_nodes", default_value_t = 100_000)]
    pub lru_size: usize,
}

impl Cli {
    fn validate(&self) {
        let mut is_ok = true;

        if self.idxdirs.is_none() && self.manifest.is_none() {
            error!("Either --idxdirs or --idx-manifest must be provided");
            is_ok = false;
        }

        if self.idxdirs.is_some() {
            let idxdirs = self.idxdirs.as_ref().unwrap();
            for idxdir in idxdirs {
                if !utils::check_index_dir(idxdir) {
                    error!(
                        "--idxdir: index directory {} is not valid",
                        idxdir.display()
                    );
                    is_ok = false;
                    break;
                }
            }
        }

        // check if the manifest file exists
        if self.manifest.is_some() {
            let manifest = self.manifest.as_ref().unwrap();
            if !manifest.exists() {
                is_ok = false;
            }

            // read the manifest file and check each index directory
            let content = std::fs::read_to_string(manifest).expect("Failed to read manifest file");
            let lines: Vec<&str> = content.lines().collect();
            if lines.len() == 0 {
                error!(
                    "--idx-manifest: manifest file {} is empty",
                    manifest.display()
                );
                is_ok = false;
            } else {
                for line in lines {
                    let idxdir = PathBuf::from(line);
                    if !utils::check_index_dir(&idxdir) {
                        is_ok = false;
                    }
                }
            }
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
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &Cli) {
    eprintln!("\nIsopedia: [Annotate provided gtf file]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct IsoformOutputRecord {
    idx: u32,
    chrom: String,
    start: u64,
    end: u64,
    length: u64,
    exon_count: u32,
    trans_id: String,
    gene_id: String,
    confidence: f64,
    detected: String,
    min_read: u32,
    positive_count: u32,
    sample_size: u32,
    attributes: String,
    format: String,
    sample_vec: Vec<String>,
    evidence_vec: Vec<u32>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Output {
    header: Vec<String>,
    sample_names: Vec<String>,
    sample_size: usize,
    records: Vec<IsoformOutputRecord>,
    record_size: usize,
    total_signals_vec: Vec<u32>,
}

impl Output {
    pub fn new(
        header: Vec<String>,
        sample_names: Vec<String>,
        total_signals_vec: Vec<u32>,
    ) -> Self {
        let sample_size = sample_names.len();
        Output {
            header,
            sample_names,
            sample_size,
            records: Vec::new(),
            record_size: 0,
            total_signals_vec,
        }
    }
}

impl OutputWriter<Output, IsoformOutputRecord> for Output {
    fn add_header(&mut self, header: &str) -> Result<()> {
        self.header.push(header.to_string());
        Ok(())
    }

    fn add_one(&mut self, item: IsoformOutputRecord) -> Result<()> {
        self.records.push(item);
        self.record_size += 1;
        Ok(())
    }

    fn merge(&mut self, other: &Self) -> Result<()> {
        if self.header != other.header {
            return Err(anyhow::anyhow!("Headers do not match, cannot merge"));
        }

        if self.record_size != other.record_size {
            return Err(anyhow::anyhow!("Sizes do not match, cannot merge"));
        }

        // extend sample names
        self.sample_names.extend_from_slice(&other.sample_names);

        for i in 0..self.record_size {
            if self.records[i].idx != other.records[i].idx {
                return Err(anyhow::anyhow!(
                    "Record idx at position {} do not match, cannot merge",
                    i
                ));
            }

            self.records[i]
                .sample_vec
                .extend_from_slice(&other.records[i].sample_vec);

            // update positive count

            self.records[i].positive_count =
                self.records[i].positive_count + other.records[i].positive_count;

            // update sample size
            self.records[i].sample_size =
                self.records[i].sample_size + other.records[i].sample_size;

            // update evidence vec
            self.records[i]
                .evidence_vec
                .extend_from_slice(&other.records[i].evidence_vec);
        }

        // update total evidence vec
        self.total_signals_vec
            .extend_from_slice(&other.total_signals_vec);

        Ok(())
    }

    fn finalize(&mut self) -> Result<()> {
        // calcuate confidence for each record
        for record in &mut self.records {
            record.confidence = isoform::MergedIsoform::get_confidence_value(
                record.evidence_vec.clone(),
                self.sample_size,
                &self.total_signals_vec,
            );
        }

        Ok(())
    }

    fn dump(&mut self, output: &PathBuf, meta_string: &str) -> Result<()> {
        self.finalize()?;

        // init writer
        let mut mywriter = MyGzWriter::new(output)?;

        // write meta
        mywriter.write_all_bytes(meta_string.as_bytes())?;

        // write header
        let mut header_line = self.header.join("\t");
        header_line.push_str("\t");
        header_line.push_str(&self.sample_names.join("\t"));
        header_line.push_str("\n");
        mywriter.write_all_bytes(header_line.as_bytes())?;

        // write records
        for record in &self.records {
            let mut line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}/{}\t{}\t{}\t",
                record.chrom,
                record.start,
                record.end,
                record.length,
                record.exon_count,
                record.trans_id,
                record.gene_id,
                record.confidence,
                record.detected,
                record.min_read,
                record.positive_count,
                record.sample_size,
                record.attributes,
                record.format
            );

            // add sample strings
            line.push_str(&record.sample_vec.join("\t"));
            line.push_str("\n");
            mywriter.write_all_bytes(line.as_bytes())?;
        }

        Ok(())
    }
}

fn main() -> Result<()> {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);


    // the index shold be put in a vec and then be merged.
    // build the vec of index directories
    let idx_dirs = if cli.manifest.is_some() {

        let manifest = cli.manifest.as_ref().unwrap();
        let content =
            std::fs::read_to_string(manifest).expect("Failed to read manifest file");
        let lines: Vec<&str> = content.lines().collect();
        let mut idx_dirs: Vec<PathBuf> = Vec::new();
        for line in lines {
            let idxdir = PathBuf::from(line);
            idx_dirs.push(idxdir);
        }
        idx_dirs

    }else {
        cli.idxdirs.as_ref().unwrap().clone()

    };

    // open gtf file

    info!("Search by gtf/gff file");
    let gtfreader: noodles_gtf::Reader<BufReader<std::fs::File>> = noodles_gtf::io::Reader::new(
        BufReader::new(std::fs::File::open(cli.gtf).expect("can not read gtf")),
    );

    // create the base output object

    // load indexes one by one

    for (i, idx_dir) in idx_dirs.iter().enumerate() {
        info!("Annotating index shard {}/{}", i+1, idx_dirs.len());

        let forest = BPForest::init(idx_dir);
        let dataset_info = DatasetInfo::load_from_file(&idx_dir.join(DATASET_INFO_FILE_NAME))?;
        let mut archive_buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
        let meta = Meta::parse(idx_dir.join(META_FILE_NAME))?;

        let gtf: TranscriptChunker<BufReader<File>> = TranscriptChunker::new(gtfreader);

        info!("Warmup index file");
        let max_gb = cli.warmup_mem * 1024 * 1024 * 1024;
        warmup(&idx_dir.clone().join(MERGED_FILE_NAME), max_gb)?;

        
    }


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

    archive_mmap
        .advise(memmap2::Advice::Sequential)
        .expect("Failed to set mmap advice");

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
            info!(
                "Processed {} transcripts",
                (batch * 10_000).to_formatted_string(&Locale::en)
            );
            iter_count = 0;
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
