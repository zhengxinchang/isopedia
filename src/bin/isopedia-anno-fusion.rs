use std::{
    collections::{HashMap, HashSet},
    env,
    fs::File,
    io::BufRead,
    path::PathBuf,
    vec,
};

use anyhow::{anyhow, Context, Result};
use clap::{command, Parser};
use isopedia::{
    bptree::BPForest,
    constants::*,
    dataset_info::DatasetInfo,
    fusion::{FusionAggrReads, FusionCluster},
    gene_index::GeneIntervalTree,
    isoform::MergedIsoform,
    isoformarchive::{self, read_record_from_mmap},
    utils,
    writer::MyGzWriter,
};

use log::{debug, error, info, warn};
use memmap2::{Advice, Mmap};
use noodles_gtf::io::Reader as gtfReader;
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-anno-fusion")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
#[clap(after_long_help = r#"
Examples:

# Annotate a fusion by providing the fusion breakpoints
isopedia-anno-fusion --idxdir /path/to/index --pos chr1:1000,chr2:2000 -o out.txt

# Annotate a list of fusions using a bed file

isopedia-anno-fusion --idxdir /path/to/index --pos-bed /path/to/fusion_breakpoints.bed -o out.txt

Format of fusion_breakpoints.bed:
chr2 \t 50795173 \t chr17 \t 61368325 \t BCAS4:BCAS3,UHR(optional)

# Discover any potential fusions using a GTF file
isopedia-anno-fusion --idxdir /path/to/index --gene-gtf /path/to/gene.gtf -o out.txt

"#)]
struct Cli {
    /// index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// two breakpoints for gene fusion to be search(-p chr1:pos1,chr2:pos2)
    #[arg(short, long)]
    pub pos: Option<String>,

    /// bed file that has the breakpoints for gene fusions. First four columns are chr1, pos1, chr2, pos2, and starts from the fifth column is the fusion id.
    #[arg(short = 'P', long)]
    pub pos_bed: Option<PathBuf>,

    /// bed file that has the start-end positions of the genes, used to find any possible gene fusions within the provided gene regions.
    #[arg(short = 'G', long)]
    pub gene_gtf: Option<PathBuf>,

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

    /// number of cached nodes for each tree in maximal
    #[arg(short='c', long="cached_nodes", default_value_t = 1_000_000)]
    pub lru_size: usize,
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

        if self.pos.is_none() && self.pos_bed.is_none() && self.gene_gtf.is_none() {
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

        if let Some(gene_gtf) = &self.gene_gtf {
            if !gene_gtf.exists() {
                error!("--gene-gtf: file {} does not exist", gene_gtf.display());
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
    archive_mmap: &Mmap,
    archive_buf: &mut Vec<u8>,
    dbinfo: &DatasetInfo,
) -> Result<()> {
    info!(
        "Processing breakpoints: {}:{}-{}:{}",
        breakpoints.0 .0, breakpoints.0 .1, breakpoints.1 .0, breakpoints.1 .1
    );
    let left_target = forest.search_one_range(&breakpoints.0 .0, breakpoints.0 .1, cli.flank, cli.lru_size);
    let right_target = forest.search_one_range(&breakpoints.1 .0, breakpoints.1 .1, cli.flank, cli.lru_size);

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

    let mut fusion_evidence_vec = vec![0u32; dbinfo.get_size()];

    // process the left targets
    let unique_left = left_target.into_iter().collect::<HashSet<_>>();
    for target in unique_left {
        let merged_isoform =
            isoformarchive::read_record_from_mmap(archive_mmap, &target, archive_buf);
        let evidence_vec = merged_isoform.find_fusion_by_breakpoints(
            &breakpoints.1 .0,
            breakpoints.1 .1,
            cli.flank,
            &dbinfo,
        );
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
            isoformarchive::read_record_from_mmap(&archive_mmap, &target, archive_buf);
        let evidence_vec = merged_isoform.find_fusion_by_breakpoints(
            &breakpoints.0 .0,
            breakpoints.0 .1,
            cli.flank,
            &dbinfo,
        );

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
        format!(
            "{}/{}",
            fusion_evidence_vec
                .iter()
                .filter(|&&x| x >= cli.min_read)
                .count()
                .to_string(),
            dbinfo.get_size()
        ),
    ];

    for idx in 0..dbinfo.get_size() {
        // if fusion_evidence_vec[idx] >= cli.min_read {
        //     record_parts.push(fusion_evidence_vec[idx].to_string());
        // } else {
        //     record_parts.push("0".to_string());
        // }
        record_parts.push(fusion_evidence_vec[idx].to_string());
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
    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    if cli.debug {
        env::set_var("RUST_LOG", "debug");
    } else {
        env::set_var("RUST_LOG", "info");
    }
    env_logger::init();

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    // let mut isofrom_archive = std::io::BufReader::new(
    //     std::fs::File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
    //         .context("Failed to open merged file")?,
    // );

    let archive_file_handle = File::open(cli.idxdir.clone().join(MERGED_FILE_NAME))
        .context("Failed to open merged file for mmap")?;

    let archive_mmap = unsafe { Mmap::map(&archive_file_handle).context("Failed to map merged file")? };
    archive_mmap.advise(Advice::Random).context("Failed to set mmap advice")?;

    let mut archive_buf = Vec::<u8>::with_capacity(1024 * 1024); // 1MB buffer

    if cli.pos.is_some() {
        let pos = &cli.pos.clone().unwrap();
        let breakpoints: BreakpointType = match process_fusion_positions(pos) {
            Ok(bp) => (bp.0, bp.1, "SingleQuery".to_string()),
            Err(e) => {
                error!("Error parsing --pos: {}", e);
                std::process::exit(1);
            }
        };

        // let mut mywriter = MyGzWriter::new(&cli.output)?;
        let mut mywriter = MyGzWriter::new(&cli.output)?;

        let mut header_str =
            String::from("chr1\tpos1\tchr2\tpos2\tid\tmin_read\tpositive/sample_size");
        let sample_name = dataset_info.get_sample_names();
        for name in sample_name {
            header_str.push_str(&format!("\t{}", name));
        }
        header_str.push('\n');
        mywriter.write_all_bytes(header_str.as_bytes())?;

        anno_single_fusion(
            breakpoints,
            &mut mywriter,
            &cli,
            &mut forest,
            &archive_mmap,
            &mut archive_buf,
            &dataset_info,
        )?;
        mywriter.finish()?;
    } else if cli.pos_bed.is_some() {
        // open the bed file
        let pos_bed = cli.pos_bed.clone().unwrap();
        let bed_file = File::open(&pos_bed).context("Failed to open the provided bed file")?;
        let reader = std::io::BufReader::new(bed_file);

        // not use the common writer since the header is different in gtf mode(discovery mode)
        let mut mywriter = MyGzWriter::new(&cli.output)?;

        let mut header_str =
            String::from("chr1\tpos1\tchr2\tpos2\tid\tmin_read\tpositive/sample_size");
        let sample_name = dataset_info.get_sample_names();
        for name in sample_name {
            header_str.push_str(&format!("\t{}", name));
        }
        header_str.push('\n');
        mywriter.write_all_bytes(header_str.as_bytes())?;

        for line in reader.lines() {
            let line = line.context("Failed to read line from bed file")?;
            let breakpoints: BreakpointType = parse_bed_line(&line)?;

            anno_single_fusion(
                breakpoints,
                &mut mywriter,
                &cli,
                &mut forest,
                &archive_mmap,
                &mut archive_buf,
                &dataset_info,
            )?;
        }
        mywriter.finish()?;
    } else if cli.gene_gtf.is_some() {
        let gtf_path = cli.gene_gtf.clone().unwrap();
        let file = std::fs::File::open(gtf_path).expect("Failed to open GTF file");
        let reader = std::io::BufReader::new(file);
        let mut gtf_reader = gtfReader::new(reader);

        let mut skipped_genes = 0;

        let gene_indexing =
            GeneIntervalTree::new(&mut gtf_reader).expect("Failed to create GeneIntervalTree");
        info!("loaded {} genes", gene_indexing.count);

        let mut mywriter = MyGzWriter::new(&cli.output)?;
        mywriter.write_all_bytes(FusionCluster::get_table_header(&dataset_info).as_bytes())?;

        let mut fusion_cluster_map = HashMap::new();

        for chromosome in gene_indexing.chroms.clone() {
            info!("Processing chromosome {}", chromosome);

            let gene_list_to_be_queried = gene_indexing.tree.get(&chromosome);
            if gene_list_to_be_queried.is_none() {
                warn!(
                    "No intervals were indexed for chromosome {}, skipping",
                    chromosome
                );
                continue;
            }

            for queried_gene_interval in gene_list_to_be_queried.unwrap() {
                let gene_interval = &queried_gene_interval.val;

                let quried_positions = gene_interval
                    .splice_sites
                    .iter()
                    .map(|x| (chromosome.clone(), *x))
                    .collect::<Vec<(String, u64)>>();

                let target = forest.search_partial_match(&quried_positions, cli.flank, 1, cli.lru_size);

                if target.is_none() {
                    skipped_genes += 1;
                    continue;
                }

                let targets = target.unwrap();

                if targets.is_empty() {
                    skipped_genes += 1;
                    continue;
                }

                for rec_ptr in targets {
                    let isoform: MergedIsoform =
                        read_record_from_mmap(&archive_mmap, &rec_ptr, &mut archive_buf);

                    let candidates = isoform.to_fusion_candidates();

                    if candidates.is_none() {
                        continue;
                    }

                    let candidates: Vec<FusionAggrReads> = candidates.unwrap();

                    for mut candidate in candidates {
                        if candidate.match_gene(&gene_indexing, cli.flank) {
                            // if gene_interval.gene_name.contains("RUNX1") {
                            //     debug!("matched RUNX1 gene: candidates: {:?}", &candidate);
                            // }
                            let record_string = candidate.get_string(&dataset_info);
                            debug!("{}", &record_string);
                            // mywriter
                            //     .write_all_bytes(record_string.as_bytes())
                            //     .context("Failed to write record string")?;

                            fusion_cluster_map
                                .entry(candidate.get_gene_hash().0)
                                .or_insert_with(|| FusionCluster::new(&candidate))
                                .add(&candidate);
                        }
                    }
                }
            }
        }

        for fusion_cluster in fusion_cluster_map.values_mut() {
            if fusion_cluster.evaluate() {
                let record_string = fusion_cluster.get_string(&dataset_info);
                mywriter
                    .write_all_bytes(record_string.as_bytes())
                    .context("Failed to write record string")?;
            }
        }

        info!("Total processed genes: {}", gene_indexing.count);
        info!("Total skipped genes: {}", skipped_genes);

        mywriter.finish()?;
    }

    info!("Total processed samples: {}", dataset_info.get_size());

    info!("Results written to {}", cli.output.display());
    info!("Finished!");
    Ok(())
}
