use std::{env, path::PathBuf};

use crate::{
    bptree::BPTree, chromosome::ChromMapping, constants::*, dataset_info::DatasetInfo, meta::Meta,
    tmpidx::Tmpindex,
};
use anyhow::Result;
use clap::Parser;
use log::{error, info, warn};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-index")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
[build index, step3] Build the tree index on the aggregated index folder.
", long_about = None)]
#[clap(after_long_help = "

Example: isopedia index -i /path/to/index/dir [--meta /path/to/meta.tsv]

The meta file is optional.

The format of the meta file is tab-separated with the first line as header.

Example of meta file:
(required) (required)  (optional)  (optional)
---------  ----------  ----------  ---------
name       path           meta1         meta2
sample1    /a/b/c         15             A
sample2    /a/b/d         37             B

")]
pub struct IndexCli {
    #[arg(
        short,
        long,
        help = "Index the directory that was generated in the aggr step."
    )]
    pub idxdir: PathBuf,

    #[arg(short, long, help = "Metadata for the samples in the index.")]
    pub meta: Option<PathBuf>,

    #[arg(
        short,
        long,
        help = "Number of threads to use. Default: 1",
        default_value_t = 4
    )]
    pub threads: usize,
}

impl IndexCli {
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
                "--idxdir: Aggr file {} does not exist in {}",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if let Some(meta_path) = &self.meta {
            if !meta_path.exists() {
                error!("--meta: meta file {} does not exist", meta_path.display());
                is_ok = false;
            }
        }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &IndexCli) {
    eprintln!("\nIsopedia: [Indexing merged isoform signals]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

pub fn run_idx(cli: &IndexCli) -> Result<()> {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    cli.validate();
    greetings(&cli);

    let chrom_bytes =
        std::fs::read(cli.idxdir.join(CHROM_FILE_NAME)).expect("Cannot read chrom map file");

    let chrom_map = ChromMapping::decode(&chrom_bytes);

    info!("Loading...");
    let mut tmpidx = Tmpindex::load(&cli.idxdir.join(TMPIDX_FILE_NAME));
    // info!("Indexing...");
    let nchrs = chrom_map.get_size();

    // info!("")

    dbg!(&tmpidx.meta);
    // change to multi-threaded building of B+ tree

    if cli.threads == 1 {
        info!("Using single thread for indexing");
        let mut cur_chr = 0;
        for chrom_id in chrom_map.get_chrom_idxs() {
            cur_chr += 1;
            info!("Indexing {}/{}th chromosome", cur_chr, nchrs);

            let tmpidx_chunker = tmpidx
                .groups(chrom_id)
                .expect("Failed to get tmpidx chunker");
            BPTree::build_tree(tmpidx_chunker, &cli.idxdir, chrom_id)?;
        }
    } else {
        info!("Using {} threads for indexing", cli.threads);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build()
            .unwrap();

        pool.install(|| {
            chrom_map
                .get_chrom_idxs()
                .par_iter()
                .enumerate()
                .for_each(|(i, &chrom_id)| {
                    info!("Indexing chromosome {}/{} ", i + 1, nchrs);

                    match tmpidx.groups(chrom_id) {
                        Some(tmpidx_chunker) => {
                            if let Err(e) =
                                BPTree::build_tree(tmpidx_chunker, &cli.idxdir, chrom_id)
                            {
                                error!("Chrom {} failed: {:?}", chrom_id, e);
                                std::process::exit(1);
                            }
                        }
                        None => {
                            info!("Skip chromosome {}: {} with zero records after filtering", chrom_id, chrom_map.get_chrom_name(chrom_id));
                            /*
                             * Explanation:
                             * In some cases, after filtering, certain chromosomes may have zero records. but the chromsome name was 
                             * still recorded in the chrom file. 
                             * After merging, the tmpidx will not have offset and count for this chromosome(because there is no record for this chromosome),
                             * thus when we try to build the B+ tree for this chromosome, it will fail. 
                             * Here we just skip the chromosome with zero records, which is a safe operation.
                             */
                            return;
                        }
                    };
                });
        });
    }

    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    if let Some(meta_path) = &cli.meta {
        info!("Integrating meta data from {}", meta_path.display());
        let meta = Meta::parse(meta_path).expect("Failed to parse meta file");
        let meta_samples = meta.get_samples();

        let dataset_samples = dataset_info.get_sample_names();
        if meta_samples.len() != dataset_samples.len() {
            error!(
                "The number of samples in the meta file ({}) does not match the number of samples in the dataset info file ({})",
                meta_samples.len(),
                dataset_samples.len()
            );
            std::process::exit(1);
        }
        for sample in meta_samples {
            if !dataset_samples.contains(&sample) {
                error!(
                    "Sample {} in the meta file does not exist in the dataset info file",
                    sample
                );
                std::process::exit(1);
            }
        }
        meta.save_to_file(&cli.idxdir.join(META_FILE_NAME))
            .expect("Failed to save meta file");
    } else {
        warn!("Meta file is not provided, Write empty meta file, you can redo index to update the metadata later.");
        let empty_meta = Meta::new_empty(dataset_info.get_sample_names());
        empty_meta
            .save_to_file(&cli.idxdir.join(META_FILE_NAME))
            .expect("Failed to save empty meta file");
    }

    info!("Finished!");
    Ok(())
}
