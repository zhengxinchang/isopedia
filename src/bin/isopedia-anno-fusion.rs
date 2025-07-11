use std::{collections::HashSet, env, fs::File, io::BufRead, path::PathBuf, vec};

use anyhow::{anyhow, Context, Result};
use clap::{command, Parser};
use isopedia::{
    bptree::BPForest, constants::*, dataset_info::DatasetInfo, isoformarchive, utils,
    writer::MyGzWriter,
};
use log::{error, info};
use serde::Serialize;

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
    pub pos: Option<String>,

    /// two breakpoints for gene fusion to be search(-p chr1:pos1,chr2:pos2)
    #[arg(short = 'P', long)]
    pub pos_bed: Option<PathBuf>,

    /// flank size for search, before and after the position
    #[arg(short, long, default_value_t = 10)]
    pub flank: u64,

    /// minimal reads to define a positive sample
    #[arg(short, long, default_value_t = 1)]
    pub min_read: u32,

    /// output file for search results
    #[arg(short, long)]
    pub output: PathBuf,

    /// debug mode
    #[arg(long, default_value_t = false)]
    pub debug: bool,
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

        if self.pos.is_none() && self.pos_bed.is_none() {
            error!("Please provide either --pos or --pos-bed");
            is_ok = false;
        }

        if let Some(pos) = &self.pos {
            match process_fusion_positions(pos) {
                Ok(_) => (),
                Err(e) => {
                    error!("Error parsing --pos: {}", e);
                    is_ok = false;
                }
            }
        }

        if let Some(pos_bed2) = &self.pos_bed {
            if !pos_bed2.exists() {
                error!("--pos-bed: file {} does not exist", pos_bed2.display());
                is_ok = false;
            }
        }

        let output_dir = self.output.parent().unwrap();

        if !output_dir.exists() {
            error!(
                "--output: parent dir {} does not exist",
                output_dir.display()
            );
            is_ok = false;
        }

        if is_ok != true {
            error!("Please check the input arguments!");
            std::process::exit(1);
        }
    }
}

fn process_fusion_positions(pos: &str) -> Result<((String, u64), (String, u64))> {
    let parts: Vec<&str> = pos.split(',').collect();

    if parts.len() != 2 {
        return Err(anyhow::anyhow!(
            "Invalid format for --pos. Expected format: chr1:pos1,chr2:pos2"
        ));
    }
    let mut breakpoints = Vec::new();
    for part in parts {
        let subparts: Vec<&str> = part.split(':').collect();
        if subparts.len() != 2 {
            return Err(anyhow::anyhow!(
                "Invalid format for --pos. Expected format: chr1:pos1,chr2:pos2"
            ));
        }
        let chr = subparts[0].to_string();
        let chr = utils::trim_chr_prefix_to_upper(&chr);

        let pos = match subparts[1].parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                return Err(anyhow::anyhow!("Position must be a valid integer"));
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

type BreakpointType = ((String, u64), (String, u64), String);

fn anno_single_fusion(
    breakpoints: BreakpointType,
    mywriter: &mut MyGzWriter,
    cli: &Cli,
    forest: &mut BPForest,
    isofrom_archive: &mut std::io::BufReader<File>,
    archive_buf: &mut Vec<u8>,
    meta: &DatasetInfo,
) -> Result<()> {
    info!(
        "Processing breakpoints: {}:{}-{}:{}",
        breakpoints.0 .0, breakpoints.0 .1, breakpoints.1 .0, breakpoints.1 .1
    );
    let left_target = forest.search_one_range(&breakpoints.0 .0, breakpoints.0 .1, cli.flank);
    let right_target = forest.search_one_range(&breakpoints.1 .0, breakpoints.1 .1, cli.flank);

    if left_target.is_none() || right_target.is_none() {
        if cli.debug {
            error!(
                "No candidates found for breakpoints: {}:{}-{}:{}",
                breakpoints.0 .0, breakpoints.0 .1, breakpoints.1 .0, breakpoints.1 .1
            );
        }
        return Ok(());
    }

    let left_target = left_target.unwrap();
    let right_target = right_target.unwrap();
    if cli.debug {
        info!(
            "Found {} candidates",
            left_target.len() + right_target.len()
        );
    }

    let mut fusion_evidence_vec = vec![0u32; MAX_SAMPLE_SIZE];

    // process the left targets
    let unique_left = left_target.into_iter().collect::<HashSet<_>>();
    for target in unique_left {
        let merged_isoform =
            isoformarchive::read_record_from_archive(isofrom_archive, &target, archive_buf);
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
    let unique_right = right_target.into_iter().collect::<HashSet<_>>();
    for target in unique_right {
        let merged_isoform: isopedia::isoform::MergedIsoform =
            isoformarchive::read_record_from_archive(isofrom_archive, &target, archive_buf);
        let evidence_vec =
            merged_isoform.find_fusion(&breakpoints.0 .0, breakpoints.0 .1, cli.flank);

        fusion_evidence_vec = fusion_evidence_vec
            .iter()
            .zip(evidence_vec.iter())
            .map(|(a, b)| a + b)
            .collect();
    }

    let mut record_parts = vec![
        breakpoints.0 .0.to_string(),
        breakpoints.0 .1.to_string(),
        breakpoints.1 .0.to_string(),
        breakpoints.1 .1.to_string(),
        breakpoints.2.to_string(),
        cli.min_read.to_string(),
        meta.get_size().to_string(),
        fusion_evidence_vec
            .iter()
            .filter(|&&x| x >= cli.min_read)
            .count()
            .to_string(),
    ];

    for idx in 0..meta.get_size() {
        if fusion_evidence_vec[idx] >= cli.min_read {
            record_parts.push(fusion_evidence_vec[idx].to_string());
        } else {
            record_parts.push("0".to_string());
        }
    }

    let record_string = record_parts.join("\t") + "\n";
    mywriter
        .write_all_bytes(record_string.as_bytes())
        .context("Failed to write record string")?;
    Ok(())
}

fn parse_bed_line(line: &str) -> Result<BreakpointType> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Err(anyhow!(
            "Invalid bed line: {}. Expected at least 5 fields.",
            line
        ));
    }
    // trim the chr prefix and convert to uppercase
    let chr1 = fields[0].to_string();
    let chr1 = utils::trim_chr_prefix_to_upper(&chr1);
    let pos1 = fields[1]
        .parse::<u64>()
        .map_err(|_| anyhow!("Invalid position"))?;
    let chr2 = fields[2].to_string();
    let chr2 = utils::trim_chr_prefix_to_upper(&chr2);
    let pos2 = fields[3]
        .parse::<u64>()
        .map_err(|_| anyhow!("Invalid position"))?;
    let fusion_id = fields[4].to_string();

    Ok(((chr1, pos1), (chr2, pos2), fusion_id))
}

fn main() -> Result<()> {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load(&cli.idxdir.join(META_FILE_NAME));

    let mut isofrom_archive = std::io::BufReader::new(
        std::fs::File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
            .context("Failed to open merged file")?,
    );

    let mut archive_buf = Vec::<u8>::with_capacity(1024 * 1024); // 1MB buffer

    let mut mywriter = MyGzWriter::new(&cli.output)?;

    let mut header_str =
        String::from("chr1\tpos1\tchr2\tpos2\tid\tmin_read\tsample_size\tpositive_sample_count");
    let sample_name = dataset_info.get_sample_names();
    for name in sample_name {
        header_str.push_str(&format!("\t{}", name));
    }
    header_str.push('\n');
    mywriter.write_all_bytes(header_str.as_bytes())?;

    if !cli.pos.is_none() {
        let pos = &cli.pos.clone().unwrap();
        let breakpoints: BreakpointType = match process_fusion_positions(pos) {
            Ok(bp) => (bp.0, bp.1, "SingleQuery".to_string()),
            Err(e) => {
                error!("Error parsing --pos: {}", e);
                std::process::exit(1);
            }
        };

        anno_single_fusion(
            breakpoints,
            &mut mywriter,
            &cli,
            &mut forest,
            &mut isofrom_archive,
            &mut archive_buf,
            &dataset_info,
        )?;
    } else {
        // open the bed file
        let pos_bed = cli.pos_bed.clone().unwrap();
        let bed_file = File::open(&pos_bed).context("Failed to open the provided bed file")?;
        let reader = std::io::BufReader::new(bed_file);

        for line in reader.lines() {
            let line = line.context("Failed to read line from bed file")?;
            let breakpoints: BreakpointType = parse_bed_line(&line)?;

            anno_single_fusion(
                breakpoints,
                &mut mywriter,
                &cli,
                &mut forest,
                &mut isofrom_archive,
                &mut archive_buf,
                &dataset_info,
            )?;
        }
    }

    mywriter.finish()?;
    info!("Results written to {}", cli.output.display());
    info!("Finished!");
    Ok(())
}
