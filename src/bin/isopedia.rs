use clap::{Parser, Subcommand};

use isopedia::cmd::fusion::run_anno_fusion;
use isopedia::cmd::fusion::AnnFusionCli;
use isopedia::cmd::index::run_idx;
use isopedia::cmd::index::IndexCli;
use isopedia::cmd::isoform::run_anno_isoform;
use isopedia::cmd::isoform::AnnIsoCli;
use isopedia::cmd::merge::run_merge;
use isopedia::cmd::merge::MergeCli;
use isopedia::cmd::profile::run_profile;
use isopedia::cmd::profile::ProfileCli;
use isopedia::cmd::splice::run_anno_splice;
use isopedia::cmd::splice::AnnSpliceCli;
use isopedia::logger::init_logger;

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
            if let Err(e) = run_profile(profile_cli) {
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
