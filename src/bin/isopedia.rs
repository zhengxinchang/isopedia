use clap::{Parser, Subcommand};

use isopedia::subcmd::aggr::run_aggr;
use isopedia::subcmd::aggr::AggrCli;
use isopedia::subcmd::extr::run_extr;
use isopedia::subcmd::extr::ExtrCli;
use isopedia::subcmd::idx::run_idx;
use isopedia::subcmd::idx::IdxCli;
use isopedia::subcmd::anno_fusion::run_anno_fusion;
use isopedia::subcmd::anno_fusion::AnnFusionCli;
use isopedia::subcmd::anno_isoform::run_anno_isoform;
use isopedia::subcmd::anno_isoform::AnnIsoCli;
use isopedia::subcmd::anno_splice::run_anno_splice;
use isopedia::subcmd::anno_splice::AnnSpliceCli;


#[derive(Parser)]
#[command(
    name = "isopedia",
    about = "Simultaneous exploration of thousands of long-read transcriptomes by read-level indexing",
    author,
    version,
    long_about = None,
    // Show authors in --help:
    help_template = "\n{about}\n\n{usage-heading}\n{usage}\n\n{all-args}\n"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Extr(ExtrCli),
    Aggr(AggrCli),
    Idx(IdxCli),
    Fusion(AnnFusionCli),
    Isoform(AnnIsoCli),
    Splice(AnnSpliceCli),
    // AnnoFusion(AnnoFusionCli),
    // AnnoSplice(AnnoSpliceCli),
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Extr(ref extr_cli) => {
            if let Err(e) = run_extr(extr_cli) {
                eprintln!("Error running extraction: {}", e);
            }
        }
        Commands::Aggr(ref aggr_cli) => {
            if let Err(e) = run_aggr(aggr_cli) {
                eprintln!("Error running aggregation: {}", e);
            }
        }
        Commands::Idx(ref idx_cli) => {
            if let Err(e) = run_idx(idx_cli) {
                eprintln!("Error running indexing: {}", e);
            }
        }
        Commands::Fusion(ref anno_fusion_cli) => {
            if let Err(e) = run_anno_fusion(anno_fusion_cli) {
                eprintln!("Error running annotation fusion: {}", e);
            }
        }
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
    }
}
