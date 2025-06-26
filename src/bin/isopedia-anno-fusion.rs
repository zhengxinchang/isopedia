use std::{
    collections::{HashMap, HashSet},
    env,
    fs::File,
    hash::Hash,
    io::{BufReader, BufWriter, Write},
    path::PathBuf,
    vec,
};

use clap::{command, Parser};
use isopedia::{
    bptree::{self, BPForest},
    constants::*,
    error::MyError,
    isoformarchive,
    meta::Meta,
    utils,
};
use itertools::Unique;
use log::{error, info};
use rust_htslib::bcf::record::Buffer;
use serde::{de::value, Serialize};

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-anno-fusion")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
#[clap(after_long_help = "
")]
struct Cli {
    /// index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// two breakpoints for gene fusion to be search(-p chr1:pos1,chr2:pos2)
    #[arg(short, long)]
    pub pos: String,

    /// flank size for search, before and after the position
    #[arg(short, long, default_value_t = 2)]
    pub flank: u64,

    /// minimal reads to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// output file for search results
    #[arg(short, long)]
    pub output: PathBuf,
}

impl Cli {
    fn validate(&self) {
        let mut is_ok = true;
        if !self.idxdir.exists() {
            error!(
                "--idxdir: index directory {} does not exist",
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join(MERGED_FILE_NAME).exists() {
            error!(
                "--idxdir: Aggr file {} does not exist in {}, please run `stix-isoform aggr` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        if !self.idxdir.join("bptree_0.idx").exists() {
            error!(
                "--idxdir: index files {} does not exist in {}, please run `stix-isoform idx` first",
                MERGED_FILE_NAME,
                self.idxdir.display()
            );
            is_ok = false;
        }

        match process_fusion_positions(&self.pos) {
            Ok(_) => {}
            Err(e) => {
                error!("Error parsing --pos: {}", e);
                is_ok = false;
            }
        }

        // if let Some(ref output) = self.output {
        //     let output_dir = output.parent().unwrap();

        //     if !output_dir.exists() {
        //         error!(
        //             "--output: parent dir {} does not exist",
        //             output_dir.display()
        //         );
        //         is_ok = false;
        //     }
        // }

        let output_dir = self.output.parent().unwrap();

        if !output_dir.exists() {
            error!(
                "--output: parent dir {} does not exist",
                output_dir.display()
            );
            is_ok = false;
        }

        if is_ok != true {
            // error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn process_fusion_positions(pos: &str) -> Result<((String, u64), (String, u64)), MyError> {
    let parts: Vec<&str> = pos.split(',').collect();

    if parts.len() != 2 {
        return Err(MyError::InvalidInput(
            "Invalid format for --pos. Expected format: chr1:pos1,chr2:pos2".to_string(),
        ));
    }
    let mut breakpoints = Vec::new();
    for part in parts {
        let subparts: Vec<&str> = part.split(':').collect();
        if subparts.len() != 2 {
            return Err(MyError::InvalidInput(
                "Invalid format for --pos. Expected format: chr1:pos1,chr2:pos2".to_string(),
            ));
        }
        let chr = subparts[0].to_string();
        let chr = utils::trim_chr_prefix_to_upper(&chr);

        let pos = match subparts[1].parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                return Err(MyError::InvalidInput(
                    "Position must be a valid integer".to_string(),
                ))
            }
        };
        breakpoints.push((chr, pos));
    }

    Ok((breakpoints[0].clone(), breakpoints[1].clone()))
}

fn greetings(args: &Cli) {
    println!("\nIsopedia: [Search provided gene fusion]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => println!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let breakpoints = match process_fusion_positions(&cli.pos) {
        Ok(bp) => bp,
        Err(e) => {
            error!("Error parsing --pos: {}", e);
            std::process::exit(1);
        }
    };

    info!(
        "Searching for fusion, breakpoints: {:?} and {:?}",
        breakpoints.0, breakpoints.1
    );

    let mut forest = BPForest::init(&cli.idxdir);
    let meta = Meta::load(&cli.idxdir.join(META_FILE_NAME));

    let mut isofrom_archive = std::io::BufReader::new(
        std::fs::File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
            .expect("Can not open aggregated records file...exit"),
    );

    let mut writer = File::create(cli.output)
        .and_then(|file| {
            let mut writer = BufWriter::new(file);

            Ok(writer)
        })
        .unwrap_or_else(|e| {
            error!("Failed to create output file: {}", e);
            std::process::exit(1);
        });

    writer
        .write_all(b"chr1\tpos1\tchr2\tpos2\tmin_read\tsample_size\tpositive_sample_count")
        .expect("Failed to write header");
    meta.get_sample_names().iter().for_each(|x| {
        writer
            .write_all(format!("\t{}", x).as_bytes())
            .expect("Failed to write sample name");
    });
    writer.write_all(b"\n").expect("Failed to write newline");

    let left_target = forest.search_one_range(&breakpoints.0 .0, breakpoints.0 .1, cli.flank);
    let right_target = forest.search_one_range(&breakpoints.1 .0, breakpoints.1 .1, cli.flank);

    if left_target.is_none() || right_target.is_none() {
        info!("No candidates found for the given breakpoints");
        return;
    }
    info!(
        "Found {} candidates",
        left_target.clone().unwrap().len() + right_target.clone().unwrap().len()
    );

    let mut fusion_evidence_vec = vec![0u32; MAX_SAMPLE_SIZE];

    // process the left targets
    let unique_left = left_target
        .clone()
        .unwrap()
        .into_iter()
        .collect::<HashSet<_>>();
    for target in unique_left {
        let merged_isoform =
            isoformarchive::read_record_from_archive(&mut isofrom_archive, &target);
        let evidence_vec =
            merged_isoform.find_fusion(&breakpoints.1 .0, breakpoints.1 .1, cli.flank);
        // dbg!(&evidence_vec);

        fusion_evidence_vec = fusion_evidence_vec
            .iter()
            .zip(evidence_vec.iter())
            .map(|(a, b)| a + b)
            .collect();
    }

    // process the right targets
    let unique_right = right_target
        .clone()
        .unwrap()
        .into_iter()
        .collect::<HashSet<_>>();
    for target in unique_right {
        let merged_isoform =
            isoformarchive::read_record_from_archive(&mut isofrom_archive, &target);
        let evidence_vec =
            merged_isoform.find_fusion(&breakpoints.0 .0, breakpoints.0 .1, cli.flank);

        fusion_evidence_vec = fusion_evidence_vec
            .iter()
            .zip(evidence_vec.iter())
            .map(|(a, b)| a + b)
            .collect();
    }

    info!(
        "Fusion evidence found in {} samples with minimal read support {}",
        fusion_evidence_vec
            .iter()
            .filter(|&&x| x >= cli.min_read)
            .count(),
        cli.min_read
    );

    // writer.write_all(b"\n").expect("Failed to write newline");

    writer
        .write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                breakpoints.0 .0,
                breakpoints.0 .1,
                breakpoints.1 .0,
                breakpoints.1 .1,
                cli.min_read,
                meta.get_size(),
                fusion_evidence_vec
                    .iter()
                    .filter(|&&x| x >= cli.min_read)
                    .count(),
            )
            .as_bytes(),
        )
        .expect("Failed to write fusion evidence");

    for idx in 0..meta.get_size() {
        if fusion_evidence_vec[idx] >= cli.min_read {
            writer
                .write_all(format!("\t{}", fusion_evidence_vec[idx]).as_bytes())
                .expect("Failed to write sample evidence");
        } else {
            writer.write_all(b"\t0").expect("Failed to write sample evidence");
        }
    }
    writer.write_all(b"\n").expect("Failed to write newline");

}
