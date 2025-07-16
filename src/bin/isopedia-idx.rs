use std::{env, path::PathBuf};

use anyhow::Result;
use clap::Parser;
use isopedia::{
    bptree::BPTree, chromosome::ChromMapping, constants::*, dataset_info::DatasetInfo, meta::Meta,
    tmpidx::Tmpindex,
};
use log::{error, info, warn};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-aggr")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
#[clap(after_long_help = "

Example: isopedia-idx -i /path/to/index/dir [--meta /path/to/meta.tsv]

The meta file is optional.

The format of the meta file is tab-separated with the first line as header.

Example of meta file:
(required)  (optional)  (optional)
---------  ----------  ---------
name       meta1         meta2
sample1     15             A
sample2     37             B

")]
struct Cli {
    #[arg(
        short,
        long,
        help = "Index the directory that was generated in the aggr step."
    )]
    pub idxdir: PathBuf,

    #[arg(short, long, help = "Metadata for the samples in the index.")]
    pub meta: Option<PathBuf>,
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

fn greetings(args: &Cli) {
    println!("\nIsopedia: [Indexing merged isoform signals]\n");
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

    let chrom_bytes =
        std::fs::read(cli.idxdir.join(CHROM_FILE_NAME)).expect("Cannot read chrom map file");

    let chrom_map = ChromMapping::decode(&chrom_bytes);

    let mut tmpidx = Tmpindex::load(&cli.idxdir.join(TMPIDX_FILE_NAME));

    for chrom_id in chrom_map.get_chrom_idxs() {
        let blocks = tmpidx.get_blocks(chrom_id);
        if blocks.len() == 0 {
            continue;
        }
        BPTree::build_tree(blocks, &cli.idxdir, chrom_id);
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
