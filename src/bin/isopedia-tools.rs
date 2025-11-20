// use anyhow::Result;
use clap::{Parser, Subcommand};
// use isopedia::chromosome::ChromMapping;
// use isopedia::constants::*;
use anyhow::Result;
#[allow(dead_code, unused)]
// use isopedia::isoformarchive::read_record_from_mmap;
use isopedia::logger::init_logger;
use isopedia::tools::output;
use isopedia::tools::profile::*;
use isopedia::tools::view::ViewArgs;
use isopedia::tools::view::*;
use isopedia::tools::{manifest::*, output::*, ToolCmdValidate};
use serde::{Deserialize, Serialize};

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
    #[command(about = "View index files")]
    view(ViewArgs),

    #[allow(non_camel_case_types)]
    #[command(about = "Merge replicates file from isopedia profile")]
    profile(MergeArgs),

    #[allow(non_camel_case_types)]
    #[command(about = "Split a large manifest file into smaller shards")]
    manifest(ManifestArg),

    #[allow(non_camel_case_types)]
    #[command(about = "Split a large manifest file into smaller shards")]
    output(OutputArg),
}

fn greetings(args: &ToolCli) {
    eprintln!("\nIsopedia: [Extract raw isoform singals from BAM/CRAM]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() -> Result<()> {
    init_logger();

    let cli = ToolCli::parse();

    greetings(&cli);

    match cli.command {
        ToolsCommands::view(view_args) => {
            if !view_args.validate() {
                std::process::exit(1);
            }

            if view_args.type_f == "tmpidx" {
                view_tmpidx(&view_args);
            } else if view_args.type_f == "archive" {
                view_archive(&view_args);
            } else if view_args.type_f == "dbinfo" {
                view_dbinfo(&view_args);
            } else if view_args.type_f == "chroms" {
                view_chroms(&view_args);
            }
        }
        ToolsCommands::profile(merge_args) => {
            if !merge_args.validate() {
                std::process::exit(1);
            }

            merge_replicates(&merge_args.input_files, &merge_args.output)?;
        }
        ToolsCommands::manifest(split_args) => {
            if !split_args.validate() {
                std::process::exit(1);
            }

            split_manifest(
                &split_args.input,
                &split_args.output_prefix,
                split_args.each_size,
            )?;
        }
        ToolsCommands::output(output_arg) => {
            if !output_arg.validate() {
                std::process::exit(1);
            }
            // Call the appropriate function based on mode
            output::merge(&output_arg)?;
        }
    }
    Ok(())
}
