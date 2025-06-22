use std::{
    env,
    io::{BufReader, BufWriter},
    path::PathBuf,
};

use clap::{command, Parser};
use isopedia::{
    bptree::BPForest, constants::*, gtf::TranscriptChunker, isoform::MergedIsoform,
    isoformarchive::read_record_from_archive, meta::Meta, tmpidx::MergedIsoformOffset,
};
use log::{error, info};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-anno")]
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

    /// positions to be search(-p chr:pos * N)
    #[arg(short, long)]
    pub pos: Vec<String>,

    /// gtf file
    #[arg(short, long)]
    pub gtf: Option<PathBuf>,

    /// flank size for search, before and after the position
    #[arg(short, long, default_value_t = 2)]
    pub flank: u64,

    /// minimal reads to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// output file for search results
    #[arg(short, long)]
    pub output: Option<PathBuf>,
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

        if let Some(ref gtf) = self.gtf {
            if !gtf.exists() {
                error!("--gtf: gtf file {} does not exist", gtf.display());
                is_ok = false;
            }
        }

        if let Some(ref output) = self.output {
            let output_dir = output.parent().unwrap();

            if !output_dir.exists() {
                error!(
                    "--output: parent dir {} does not exist",
                    output_dir.display()
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

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let mut forest = BPForest::init(&cli.idxdir);
    let meta = Meta::load(&cli.idxdir.join(META_FILE_NAME));

    if cli.pos.len() > 0 {
        info!("Search by a list of positions");
    } else if !cli.gtf.is_none() {
        info!("Search by gtf file");
        let gtfreader = noodles_gtf::io::Reader::new(BufReader::new(
            std::fs::File::open(cli.gtf.unwrap()).expect("can not read gtf"),
        ));

        let gtf = TranscriptChunker::new(gtfreader);

        // let mut aggr_file = cli.idxdir.clone();
        // aggr_file.push(MERGED_FILE_NAME);

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

        writer
            .write(
                "chrom\tstart\tend\ttrans_id\tgene_id\thit\tmin_read\tpositive_count/sample_size\tattributes"
                    .as_bytes(),
            )
            .unwrap();
        meta.get_sample_names().iter().for_each(|x| {
            writer.write(format!("\t{}", x).as_bytes()).unwrap();
        });
        writer.write("\n".as_bytes()).unwrap();

        let mut isofrom_archive = std::io::BufReader::new(
            std::fs::File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
                .expect("Can not open aggregated records file...exit"),
        );

        let mut iter_count = 0;
        let mut batch = 0;

        let mut total_acc_evidence_flag_vec = vec![0u32; meta.get_size()];
        for trans in gtf {
            iter_count += 1;

            if iter_count == 100000 {
                batch += 1;
                info!("Processed {} transcripts", iter_count * batch );
                iter_count = 0;
            }

            let mut queries: Vec<(String, u64)> = trans.get_quieries();
            queries.sort_by_key(|x| x.1);
            let res = forest.search_multi(&queries, cli.flank);

            let target: Vec<MergedIsoformOffset> = res
                .into_iter()
                .filter(|x| x.n_splice_sites == queries.len() as u32)
                .collect();

            if target.len() == 0 {
                miss_count += 1;
            } else {
                hit_count += 1;
            }

            if target.len() > 0 {
                // dbg!(target.len());

                let mut acc_pos_count = 0;
                let mut acc_sample_evidence_arr = vec![0u32; meta.get_size()];

                for i in 0..target.len() {
                    let record: MergedIsoform =
                        read_record_from_archive(&mut isofrom_archive, &target[i]);
                    acc_pos_count += record.get_positive_count(&cli.min_read);
                    acc_sample_evidence_arr = acc_sample_evidence_arr
                        .iter()
                        .zip(record.get_sample_evidence_arr().iter())
                        .map(|(a, b)| a + b)
                        .collect();
                }

                acc_sample_evidence_arr
                    .iter()
                    .enumerate()
                    .for_each(|(i, x)| {
                        if *x > 0 {
                            total_acc_evidence_flag_vec[i] += 1;
                        }
                    });

                writer
                    .write(
                        format!(
                            "{}\t{:?}\t{:?}\t{}\t{}\tyes\t{}\t{}\t{}\t{}\n",
                            trans.chrom,
                            trans.start,
                            trans.end,
                            trans.trans_id,
                            trans.gene_id,
                            &cli.min_read,
                            format!(
                                "{}/{}",
                                acc_pos_count,
                                // record.sample_size
                                meta.get_size()
                            ),
                            trans.get_attributes(),
                            acc_sample_evidence_arr
                                .into_iter()
                                .map(|x| { x.to_string() })
                                .collect::<Vec<String>>()
                                .join("\t")
                        )
                        .as_bytes(),
                    )
                    .unwrap();
            } else {
                writer
                    .write(
                        format!(
                            "{}\t{:?}\t{:?}\t{}\t{}\tno\t{}\tNA\t{}\t{}\n",
                            trans.chrom,
                            trans.start,
                            trans.end,
                            trans.trans_id,
                            trans.gene_id,
                            &cli.min_read,
                            trans.get_attributes(),
                            meta.get_sample_names()
                                .iter()
                                .map(|_| "0")
                                .collect::<Vec<&str>>()
                                .join("\t")
                        )
                        .as_bytes(),
                    )
                    .unwrap();
            }
        }
        let total = hit_count + miss_count;
        info!("Processed {} transcripts", 100000 * batch + iter_count);
        info!("Sample-wide stats: ");
        info!("> Sample\thit\tmiss\tpct");
        for i in 0..meta.get_size() {
            info!(
                "> {:}\t{}\t{}\t {:.2}%",
                meta.get_sample_names()[i],
                total_acc_evidence_flag_vec[i],
                total - total_acc_evidence_flag_vec[i],
                total_acc_evidence_flag_vec[i] as f64 / total as f64 * 100f64
            );
        }
        info!(
            "Index-wide stats: hit: {}, miss: {}, total: {}, pct: {:.2}%",
            hit_count,
            miss_count,
            hit_count + miss_count,
            hit_count as f64 / (hit_count + miss_count) as f64 * 100f64
        );
    }

    info!("Finished");
}
