use anyhow::Result;
use clap::Parser;
use log::{error, info, warn};
#[allow(dead_code, unused)]
// use memmap2::Mmap;
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::myio::MyGzWriter;
use crate::reads::AggrRead;

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

    let base_file = files[0].clone();

    //open gz
    let mut buf = String::new();
    let mut gzrdr = BufReader::new(flate2::read::GzDecoder::new(std::fs::File::open(
        base_file,
    )?));
    // merge chroms
    gzrdr.read_line(&mut buf)?;
    let chroms = buf.trim().split('\t').collect::<Vec<&str>>();

    for chrom in chroms {
        chrom_set.insert(chrom.to_string());
    }

    // skip the second line(header line)
    gzrdr.read_line(&mut buf)?;

    // start to build the map of agg-read

    let mut merged_map = std::collections::HashMap::new();

    buf.clear();
    while let Ok(n) = gzrdr.read_line(&mut buf) {
        if n == 0 {
            break;
        }

        total_count += 1;

        let agg_read = AggrRead::from_string(&buf);
        merged_map.insert(agg_read.signature.clone(), agg_read);
        buf.clear();
    }

    // iteraterlly add agg_read from the rest of files

    for file in &files[1..] {
        let mut buf = String::new();
        let mut gzrdr = BufReader::new(flate2::read::GzDecoder::new(std::fs::File::open(file)?));

        // read first line for choromsome
        gzrdr.read_line(&mut buf)?;
        let chroms = buf.trim().split('\t').collect::<Vec<&str>>();
        for chrom in chroms {
            if !chrom_set.contains(chrom) {
                warn!("Chromosome {} is not present in all input files, this might happen in non-main contigs, please be careful if you see this message for any main chromosomes.", chrom);
            }
            chrom_set.insert(chrom.to_string());
        }

        // skip the header
        gzrdr.read_line(&mut buf)?;
        buf.clear();
        // start to build the map of agg-read
        while let Ok(n) = gzrdr.read_line(&mut buf) {
            if n == 0 {
                break;
            }

            total_count += 1;

            let agg_read = AggrRead::from_string(&buf);

            merged_map
                .entry(agg_read.signature.clone())
                .and_modify(|e| e.merge(agg_read.clone()))
                .or_insert(agg_read);

            buf.clear();
        }
    }

    info!("{} Chromosomes detected and merged.", chrom_set.len());

    // write the chrom_set
    info!("Processing chromosomes");
    let mut chrom_str = Vec::new();
    chrom_set.into_iter().for_each(|x| {
        chrom_str.extend_from_slice(x.as_bytes());
        chrom_str.push(b'\t');
    });
    chrom_str.push(b'\n');

    info!("Write to output file");

    let mut mywriter = MyGzWriter::new(output).expect(&format!(
        "Can not create output file {} .",
        output.display()
    ));

    info!(
        "Merged {} records into {} agg-reads from {} files",
        total_count,
        &merged_map.len(),
        files.len()
    );

    chrom_str.extend_from_slice(b"signature\tevidence\tchrom\tsplice_junctions\tisoform_diffs\n");
    mywriter
        .write_all_bytes(&chrom_str)
        .expect("can not write headers..");

    let mut sorted_agg_isoform_map = merged_map.into_iter().collect::<Vec<(u64, AggrRead)>>();
    sorted_agg_isoform_map.sort_by_key(|x| x.0);
    for (_, agg_isoform) in sorted_agg_isoform_map.iter() {
        // writer.write(agg_isoform.to_record().as_bytes()).unwrap();
        mywriter
            .write_all_bytes(agg_isoform.to_record().as_bytes())
            .expect("can not write record...");
    }

    info!("Finished");

    Ok(())
}
