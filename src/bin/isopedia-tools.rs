use anyhow::Result;
use clap::{Parser, Subcommand};
use isopedia::chromosome::ChromMapping;
use isopedia::constants::*;
#[allow(dead_code, unused)]
use isopedia::isoformarchive::read_record_from_mmap;
use isopedia::logger::init_logger;
use isopedia::reads::AggrRead;
#[allow(dead_code, unused)]
use isopedia::tmpidx::Tmpindex;
use isopedia::tools::inspection::InspectArgs;
// use isopedia::tools::ToolCmdValidate;
// use isopedia::writer::MyGzWriter;
// use log::{error, info, warn};
// #[allow(dead_code, unused)]
// use memmap2::Mmap;
// use noodles_fasta::fai::read;
use serde::{Deserialize, Serialize};
// use std::collections::HashSet;
// use std::fs::File;
// use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use isopedia::tools::inspection::*;
use isopedia::tools::merge_replicates::*;
use isopedia::tools::{manifest::*, ToolCmdValidate};

#[derive(Parser, Clone, Debug, Serialize, Deserialize)]
#[command(name = "isopedia-tool")]
#[command(
    about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>, <xinchang.zheng@bcm.edu>
",
    long_about = "This subtool provides several utilities for working with Isopedia."
)]
#[command(author = "Xinchang Zheng", version)]
pub struct ToolCli {
    #[command(subcommand)]
    pub command: ToolsCommands,
}

#[derive(Subcommand, Clone, Debug, Serialize, Deserialize)]
pub enum ToolsCommands {
    #[allow(non_camel_case_types)]
    #[command(about = "Inspect index files")]
    inspect(InspectArgs),

    #[allow(non_camel_case_types)]
    #[command(about = "Merge replicates file from isopedia-extr")]
    merge(MergeArgs),

    #[allow(non_camel_case_types)]
    #[command(about = "Split a large manifest file into smaller shards")]
    split_manifest(ManifestSplit),
}

fn greetings(args: &ToolCli) {
    eprintln!("\nIsopedia: [Extract raw isoform singals from BAM/CRAM]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() {
    init_logger();

    let cli = ToolCli::parse();

    greetings(&cli);

    match cli.command {
        ToolsCommands::inspect(inspec_args) => {
            if !inspec_args.validate() {
                std::process::exit(1);
            }

            if inspec_args.type_f == "tmpidx" {
                inspect_intrim_file(&inspec_args.idx, &inspec_args.output);
            } else if inspec_args.type_f == "archive" {
                inspect_archive(&inspec_args.idx, &inspec_args.output);
            } else if inspec_args.type_f == "dbinfo" {
                inspect_meta(&inspec_args.idx);
            } else if inspec_args.type_f == "chroms" {
                inspect_chroms(&inspec_args.idx);
            }
        }
        ToolsCommands::merge(merge_args) => {
            if !merge_args.validate() {
                std::process::exit(1);
            }

            merge_replicates(&merge_args.input_files, &merge_args.output)
                .expect("Failed to merge replicates");
        }
        ToolsCommands::split_manifest(split_args) => {
            if !split_args.validate() {
                std::process::exit(1);
            }

            split_manifest(
                &split_args.input,
                &split_args.output_prefix,
                split_args.each_size,
            )
            .expect("Failed to split manifest");
        }
    }
}
