use anyhow::Result;
use clap::Parser;
use log::{error, info, warn};
#[allow(dead_code, unused)]
// use memmap2::Mmap;
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::PathBuf;

use crate::myio::MyGzWriter;
use crate::reads::{profile_header, AggrRead, SingleSampleReader};

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "

Merge replicate files from isopedia-extr into a single file.

This is useful when you have multiple replicate files that you want to combine into a single file.

Example:

isopedia-tool merge -i <INPUT_FILE 1> <INPUT_FILE 2> -o <OUTPUT>


At least two input files must be specified.
")]
pub struct MergeArgs {
    /// input files that generated from isopedia-extr
    #[arg(short, long, num_args = 2..)]
    pub input_files: Vec<PathBuf>,

    /// output merged file in gz format
    #[arg(short, long)]
    pub output: PathBuf,
}

impl MergeArgs {
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if self.input_files.len() < 2 {
            error!("--input: at least two input files must be specified");
            is_ok = false;
        }

        is_ok
    }
}

pub fn merge_replicates(files: &Vec<PathBuf>, output: &PathBuf) -> Result<()> {
    let mut chrom_set: HashSet<String> = HashSet::new();
    let mut total_count = 0;

    if files.len() == 0 || files.len() == 1 {
        error!("At least two files are required to merge replicates");
        //exit
        std::process::exit(1);
    }

    for file in files {
        if !file.exists() {
            error!("File {:?} does not exist", file);
            std::process::exit(1);
        }
    }

    let mut merged_map: std::collections::HashMap<u128, AggrRead> = std::collections::HashMap::new();
    for (i, file) in files.iter().enumerate() {
        let mut reader = SingleSampleReader::new(file)?;
        for chrom in &reader.chroms {
            if i > 0 && !chrom_set.contains(chrom) {
                warn!("Chromosome {} is not present in all input files", chrom);
            }
            chrom_set.insert(chrom.clone());
        }
        while let Some(agg_read) = reader.next_rec()? {
            total_count += 1;
            match merged_map.entry(agg_read.signature) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().merge(agg_read)
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(agg_read);
                }
            }
        }
    }

    info!("{} Chromosomes detected and merged.", chrom_set.len());

    // write the chrom_set
    info!("Processing chromosomes");
    let chrom_str = profile_header(chrom_set.iter().map(String::as_str));
    info!("Write to output file");
    let mut mywriter = MyGzWriter::new(output)?;

    info!(
        "Merged {} records into {} agg-reads from {} files",
        total_count,
        &merged_map.len(),
        files.len()
    );

    mywriter.write_all_bytes(chrom_str.as_bytes())?;

    let mut sorted_agg_isoform_map = merged_map.into_iter().collect::<Vec<(u128, AggrRead)>>();
    sorted_agg_isoform_map.sort_by_key(|x| x.0);
    for (_, agg_isoform) in sorted_agg_isoform_map.iter() {
        // writer.write(agg_isoform.to_record().as_bytes()).unwrap();
        mywriter.write_all_bytes(agg_isoform.to_record().as_bytes())?;
    }

    mywriter.finish()?;
    info!("Finished");

    Ok(())
}
