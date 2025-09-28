use clap::{Parser, Subcommand};

use isopedia::logger::init_logger;
use isopedia::subcmd::anno_fusion::run_anno_fusion;
use isopedia::subcmd::anno_fusion::AnnFusionCli;
use isopedia::subcmd::anno_isoform::run_anno_isoform;
use isopedia::subcmd::anno_isoform::AnnIsoCli;
use isopedia::subcmd::anno_splice::run_anno_splice;
use isopedia::subcmd::anno_splice::AnnSpliceCli;
use isopedia::subcmd::index::run_idx;
use isopedia::subcmd::index::IndexCli;
use isopedia::subcmd::merge::run_merge;
use isopedia::subcmd::merge::MergeCli;
use isopedia::subcmd::profile::run_extr;
use isopedia::subcmd::profile::ProfileCli;

#[derive(Parser)]
#[command(
    name = "isopedia",
    about = "[Isopedia] Simultaneous exploration of thousands of long-read transcriptomes by read-level indexing\n\nRepository: https://github.com/zhengxinchang/isopedia \nContact: Xinchang Zheng <zhengxc93@gmail.com>",
    author,
    version,
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Isoform(AnnIsoCli),
    Splice(AnnSpliceCli),
    Fusion(AnnFusionCli),
    Profile(ProfileCli),
    Merge(MergeCli),
    Index(IndexCli),
}

fn main() {
    init_logger();

    let cli = Cli::parse();

    match cli.command {
        Commands::Isoform(ref anno_isoform_cli) => {
            if let Err(e) = run_anno_isoform(anno_isoform_cli) {
                eprintln!("Error running annotation isoform: {}", e);
            }
        }
        Commands::Splice(ref anno_splice_cli) => {
            if let Err(e) = run_anno_splice(anno_splice_cli) {
                eprintln!("Error running annotation splice: {}", e);
            }
        }
        Commands::Fusion(ref anno_fusion_cli) => {
            if let Err(e) = run_anno_fusion(anno_fusion_cli) {
                eprintln!("Error running annotation fusion: {}", e);
            }
        }
        Commands::Profile(ref profile_cli) => {
            if let Err(e) = run_extr(profile_cli) {
                eprintln!("Error running extraction: {}", e);
            }
        }
        Commands::Merge(ref merge_cli) => {
            if let Err(e) = run_merge(merge_cli) {
                eprintln!("Error running merge: {}", e);
            }
        }
        Commands::Index(ref idx_cli) => {
            if let Err(e) = run_idx(idx_cli) {
                eprintln!("Error running indexing: {}", e);
            }
        }
    }
}
