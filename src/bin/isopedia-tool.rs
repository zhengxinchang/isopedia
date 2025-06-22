use clap::{Parser, Subcommand};
use isopedia::isoformarchive::read_record_from_archive;
use log::error;
use serde::Serialize;
use std::os::unix::process;
use std::path::PathBuf;
use std::io::Write;
use isopedia::tmpidx::Tmpindex;
use isopedia::chromosome::ChromMapping;
use isopedia::constants::*;

#[derive(Parser, Clone, Debug)]
#[command(name = "isopedia-tool")]
#[command(about = "Auxiliary tools for Isopedia")]
#[command(author = "Xinchang Zheng", version)]
pub struct ToolCli {
    #[command(subcommand)]
    pub command: ToolsCommands,
}

#[derive(Subcommand, Clone, Debug)]
pub enum ToolsCommands {
    #[allow(non_camel_case_types)]
    #[command(about = "Inspect index files")]
    inspect(InspectArgs),
}

pub trait Validate {
    fn validate(&self) -> bool;
}

#[derive(Parser, Debug, Clone)]
pub struct InspectArgs {
    /// index directory
    #[arg(short, long)]
    pub idx: PathBuf,

    /// type of file to be inspect can be either 'tmpidx' or 'archive'
    #[arg(short, long, name = "type")]
    pub type_f: String,

    /// Output file name
    #[arg(short, long)]
    pub output: PathBuf,
}

impl Validate for InspectArgs {
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.idx.exists() {
            error!("--idx: index directory does not exist");
            is_ok = false;
        }

        if self.type_f != "tmpidx" && self.type_f != "archive" {
            error!("--type: type must be either 'tmpidx' or 'archive'");
            is_ok = false;
        }

        if !self.output.parent().unwrap().exists() {
            error!(
                "--output: output dir {} does not exist",
                self.output.parent().unwrap().display()
            );
            is_ok = false;
        }
        is_ok
    }
}



pub fn inspect_intrim_file(idx: &PathBuf, output: &PathBuf) {
    let mut interim = Tmpindex::load(&idx.join(TMPIDX_FILE_NAME));
    let chrom_bytes = std::fs::read(&idx.join(CHROM_FILE_NAME)).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);
    for chrom_id in chromamp.get_chrom_idxs() {
        let chrom_name = chromamp.get_chrom_name(chrom_id);
        let blocks = interim.get_blocks(chrom_id);
        if blocks.len() == 0 {
            continue;
        }
        for block in blocks {
            for record in block {
                writeln!(
                    writer,
                    "{}:{}-{:?}",
                    chrom_name, &record.pos, record.record_ptr_vec
                )
                .unwrap();
            }
        }
    }
}

pub fn inspect_aggr_dat(idx: &PathBuf, output: &PathBuf) {

    let mut processed_signautres = std::collections::HashSet::new();

    let mut interim = Tmpindex::load(&idx.join("interim.idx"));
    let chrom_bytes = std::fs::read(&idx.join("chrom.map")).unwrap();
    let chromamp = ChromMapping::decode(&chrom_bytes);
    // let blocks = idx.get_blocks(chrom_id);
    let writer = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(writer);

    let mut aggr_reader = std::io::BufReader::new(
        std::fs::File::open(idx.join(MERGED_FILE_NAME))
            .expect("Can not open aggregated records file...exit"),
    );

    for chrom_id in chromamp.get_chrom_idxs() {
        // let chrom_name = chromamp.get_chrom_name(chrom_id);
        let blocks = interim.get_blocks(chrom_id);
        if blocks.len() == 0 {
            continue;
        }
        for block in blocks {
            for record_grp in block {
                for record_ptr in record_grp.record_ptr_vec {
                    let rec  = read_record_from_archive(&mut aggr_reader, &record_ptr);
                    if processed_signautres.contains(&rec.signature) {
                        continue;
                    }
                    processed_signautres.insert(rec.signature);
                    writeln!(writer, "{}", rec.get_string()).unwrap();
                }
            }
        }
    }
}


fn main() {
    env_logger::Builder::from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    )
    .format(|buf, refcord| {
        writeln!(
            buf,
            "[{}] [{}]: {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            refcord.level(),
            refcord.args()
        )
    })
    .init();

    let cli = ToolCli::parse();
    match cli.command {
        ToolsCommands::inspect(inspec_args) => {
            if !inspec_args.validate() {
                std::process::exit(1);
            }

            if inspec_args.type_f == "tmpidx" {
                inspect_intrim_file(&inspec_args.idx, &inspec_args.output);
            } else if inspec_args.type_f == "archive" {
                inspect_aggr_dat(&inspec_args.idx, &inspec_args.output);
            }
        }
    }
}
