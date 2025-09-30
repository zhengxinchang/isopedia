use crate::tools::ToolCmdValidate;
use anyhow::Result;
use clap::Parser;
use log::{error, info};
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(after_long_help = "

Split a large manifest file into smaller shards.
Examples:
isopedia-tool split-manifest --input <INPUT> --output_prefix <OUTPUT_PREFIX> --each_size <EACH_SIZE>
")]
pub struct ManifestSplit {
    #[arg(short, long)]
    pub input: PathBuf,
    #[arg(short, long)]
    pub output_prefix: PathBuf,
    #[arg(short = 'c', long)]
    pub each_size: usize,
}

impl ToolCmdValidate for ManifestSplit {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.input.exists() {
            error!("--input: input file does not exist");
            is_ok = false;
        }

        if self.each_size == 0 {
            error!("--each_size: each_size must be greater than 0");
            is_ok = false;
        }

        is_ok
    }
}

/// merge seprated manifest files into one
pub fn split_manifest(file: &PathBuf, outprefix: &PathBuf, each_size: usize) -> Result<()> {
    // open the file
    let mut reader = BufReader::new(File::open(file)?);

    let mut header = String::new();
    reader.read_line(&mut header)?;
    let header = header.trim_end().to_string();

    let mut file_count = 0;
    let mut shards: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if file_count % each_size == 0 {
            if file_count > 0 {
                // write the previous shard

                let mut out_file = outprefix.clone();
                out_file.set_extension(format!("shard_{}.tsv", file_count / each_size));
                let mut writer = BufWriter::new(File::create(&out_file)?);
                writer.write_all(header.as_bytes())?;
                writer.write_all(b"\n")?;
                for shard_line in &shards {
                    writer.write_all(shard_line.as_bytes())?;
                    writer.write_all(b"\n")?;
                }
                info!("Wrote shard file: {:?}", out_file);
                shards.clear();
            }
        }
        shards.push(line);
        file_count += 1;
    }
    // write the last shard
    if shards.len() > 0 {
        let mut out_file = outprefix.clone();
        out_file.set_extension(format!("shard_{}.tsv", file_count / each_size + 1));
        let mut writer = BufWriter::new(File::create(&out_file)?);
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
        for shard_line in &shards {
            writer.write_all(shard_line.as_bytes())?;
            writer.write_all(b"\n")?;
        }
        info!("Wrote shard file: {:?}", out_file);
        shards.clear();
    }

    Ok(())
}

#[allow(dead_code)]
fn merge_manifest(files: &Vec<PathBuf>, out_path: &PathBuf) -> Result<()> {
    todo!("Implement merge_manifest");
    Ok(())
}
