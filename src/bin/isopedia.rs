use clap::{Parser, Subcommand};

use isopedia::cmd::download;
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
    about = "Isopedia: Simultaneous exploration of thousands of long-read transcriptomes by read-level indexing\n\nRepository: https://github.com/zhengxinchang/isopedia \nContact: Xinchang Zheng <zhengxc93@gmail.com>",
    author,
    version,
    after_long_help = "
    
    > [Query examples]

    # query isoform
    isopedia isoform -i /path/to/index -g query.gtf -o output_isoform.tsv.gz

    # query fusion breakpoints
    isopedia fusion --idxdir /path/to/index --pos chr1:1000,chr2:2000 -o single.tsv.gz
    isopedia fusion --idxdir /path/to/index --pos-bed /path/to/fusion_breakpoints.bed -o batch.tsv.gz # batch query

    # query gene regions to idenitify potential fusion events
    isopedia fusion --idxdir /path/to/index --gene-gtf /path/to/gene.gtf -o out.tsv.gz

    # query splice junctions
    isopedia splice --idxdir /path/to/index -s chr1:1000,chr2:2000 -o single.tsv.gz
    isopedia splice --idxdir /path/to/index -S /path/to/sj.bed -o batch.tsv.gz # batch query
    isopedia-splice-viz.py -i single.tsv.gz -g gencode.basic.gtf -o splice_viz.html # visualize splice junctions


    > [Build your own index]

    # step1: profile single sample
    isopedia profile -i /path/to/sample.bam -o /path/to/profile.tsv.gz

    # step2: merge multiple profiles
    printf \"name \\t path \\n\" > manifest.tsv
    printf \"sample1 \\t /path/to/sample1_profile.tsv.gz \\n\" >> manifest.tsv
    printf \"sample2 \\t /path/to/sample2_profile.tsv.gz \\n\" >> manifest.tsv
    isopedia merge -i manifest.tsv -o index/ 

    # step3: build index
    isopedia index -i index/ -m manifest.tsv

    > [Full documentation]

    https://github.com/zhengxinchang/isopedia
    "
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

    // output the tool version
    println!("Isopedia version: {}", env!("CARGO_PKG_VERSION"));

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
