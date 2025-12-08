use std::{path::PathBuf, sync::Mutex};

use crate::{
    bptree::BPTree,
    chromosome::{ChromMapping, ChromMappingHelper},
    constants::*,
    dataset_info::DatasetInfo,
    meta::Meta,
    myio::GeneralOutputIO,
    tmpidx::Tmpindex,
};
use anyhow::Result;
use clap::Parser;
use log::{error, info};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia index")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[build index, step3] Build the tree index on the aggregated index folder.
", long_about = None)]
#[clap(after_long_help = "

Example: isopedia index -i /path/to/index/dir --manifest /path/to/meta.tsv

The manifest file is the file that was used in the merge step, columns more than 2 will be treated as metadata columns.

The format of the manifest file is tab-separated with the first line as header.

Example of manifest file:
(required) (required)  (optional)  (optional) ....
---------  ----------  ----------  ---------  ----
name       path           meta1         meta2 ....
sample1    /a/b/c         15             A    ....
sample2    /a/b/d         37             B    ....

")]
pub struct IndexCli {
    #[arg(
        short,
        long,
        help = "Index the directory that was generated in the merge step."
    )]
    pub idxdir: PathBuf,

    #[arg(short, long, help = "Manifest file for the samples in the index.")]
    pub manifest: PathBuf,

    #[arg(short, long, help = "Number of threads to use.", default_value_t = 4)]
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

        // if let Some(meta_path) = &self.meta {
        if !self.manifest.exists() {
            error!(
                "--meta: meta file {} does not exist",
                self.manifest.display()
            );
            is_ok = false;
        } else {
            let meta = Meta::from_file(&self.manifest, None, None)
                .expect("Cannot parse the manifest file");
            meta.validate_sample_names();
            meta.validate_path_column();
        }
        // }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn greetings(args: &IndexCli) {
    // eprintln!("\nIsopedia: [Indexing merged isoform signals]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

pub fn run_idx(cli: &IndexCli) -> Result<()> {
    // env::set_var("RUST_LOG", "info");
    // env_logger::init();

    greetings(&cli);
    cli.validate();

    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    // if let Some(meta_path) = &cli.meta {
    info!("Integrating meta data from {}", cli.manifest.display());
    let mut meta = Meta::parse(&cli.manifest, None).expect("Failed to parse meta file");
    meta.remove_path_column();
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
    // } else {
    //     warn!("Meta file is not provided, Write empty meta file, you can redo index to update the metadata later.");
    //     let empty_meta = Meta::new_empty(dataset_info.get_sample_names());
    //     empty_meta
    //         .save_to_file(&cli.idxdir.join(META_FILE_NAME))
    //         .expect("Failed to save empty meta file");
    // }

    let chrom_bytes =
        std::fs::read(cli.idxdir.join(CHROM_FILE_NAME)).expect("Cannot read chrom map file");

    let mut chrom_map = ChromMapping::decode(&chrom_bytes);

    info!("Loading...");
    let tmpidx = Tmpindex::load(&cli.idxdir.join(TMPIDX_FILE_NAME));
    // info!("Indexing...");
    let nchrs = chrom_map.get_size();

    if cli.threads == 1 {
        info!("Using single thread for indexing");
        let mut chrom_helper = ChromMappingHelper::new();
        let mut cur_chr = 0;
        for chrom_id in chrom_map.get_chrom_idxs() {
            info!("Indexing {}/{}th chromosome", cur_chr, nchrs);

            match tmpidx.groups(chrom_id) {
                Some(tmpidx_chunker) => {
                    if let Err(e) = BPTree::build_tree(tmpidx_chunker, &cli.idxdir, chrom_id) {
                        error!("Chrom {} failed: {:?}", chrom_id, e);
                        std::process::exit(1);
                    }
                    cur_chr += 1;
                    chrom_helper.add_record2chrom(chrom_id);
                }
                None => {
                    info!(
                        "Skip chromosome {}: {} with zero records after filtering",
                        chrom_id,
                        chrom_map.get_chrom_name(chrom_id)
                    );
                    /*
                     * Explanation:
                     * In some cases, after filtering, certain chromosomes may have zero records. but the chromsome name was
                     * still recorded in the chrom file.
                     * After merging, the tmpidx will not have offset and count for this chromosome(because there is no record for this chromosome),
                     * thus when we try to build the B+ tree for this chromosome, it will fail.
                     * Here we just skip the chromosome with zero records, which is a safe operation.
                     */
                    continue;
                }
            }
        }

        chrom_helper.drop_zero(&mut chrom_map);
    } else {
        info!("Using {} threads for indexing", cli.threads);
        let chrom_helper = Mutex::new(ChromMappingHelper::new());

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build()?;

        pool.install(|| {
            chrom_map
                .get_chrom_idxs()
                .par_iter()
                .enumerate()
                .for_each(|(i, &chrom_id)| {
                    info!("Indexing chromosome #{}, total {}", i, nchrs);

                    match tmpidx.groups(chrom_id) {
                        Some(tmpidx_chunker) => {
                            if let Err(e) =
                                BPTree::build_tree(tmpidx_chunker, &cli.idxdir, chrom_id)
                            {
                                error!("Chrom {} failed: {:?}", chrom_id, e);
                                std::process::exit(1);
                            }
                            let mut helper = chrom_helper.lock().unwrap();
                            helper.add_record2chrom(chrom_id);
                        }
                        None => {
                            info!(
                                "Skip chromosome {}: {} with zero records after filtering",
                                chrom_id,
                                chrom_map.get_chrom_name(chrom_id)
                            );
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

        let mut helper = chrom_helper.lock().unwrap();
        helper.drop_zero(&mut chrom_map);
    }

    // rewrite the chrom file to drop zero-record chromosomes

    let chrom_encoded = chrom_map.encode();
    std::fs::write(cli.idxdir.join(CHROM_FILE_NAME), &chrom_encoded)
        .expect("Failed to write updated chrom file");

    info!("Finished!");
    Ok(())
}
