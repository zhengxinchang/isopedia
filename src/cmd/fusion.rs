use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufRead,
    path::PathBuf,
    vec,
};

use crate::{
    bptree::BPForest,
    breakpoints,
    constants::*,
    dataset_info::DatasetInfo,
    fusion::{FusionAggrReads, FusionCluster},
    gene_index::GeneIntervalTree,
    isoform::MergedIsoform,
    isoformarchive::ArchiveCache,
    meta::Meta,
    myio::*,
    results::TableOutput,
    utils::{self, greetings2},
};
use anyhow::{anyhow, Context, Result};
use clap::Parser;

use log::{error, info, warn};
use noodles_gtf::io::Reader as gtfReader;
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia fusion")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[Query/discovery] Annotate provided fusion breakpoints or discover potential gene fusion events(from provided GTF file) with the index.
", long_about = None)]
#[clap(after_long_help = r#"
Examples:

# Annotate a fusion by providing the fusion breakpoints
isopedia fusion --idxdir /path/to/index --pos chr1:1000,chr2:2000 -o out.txt

# Annotate a list of fusions using a bed file

isopedia fusion --idxdir /path/to/index --pos-bed /path/to/fusion_breakpoints.bed -o out.txt

Format of fusion_breakpoints.bed:
chr2 \t 50795173 \t 50795273 \t chr17 \t 61368325 \t 61368425 \t BCAS4:BCAS3,UHR(optional)

# Discover any potential fusions using a GTF file
isopedia fusion --idxdir /path/to/index --gene-gtf /path/to/gene.gtf -o out.txt

"#)]
pub struct AnnFusionCli {
    /// index directory
    #[arg(short, long)]
    pub idxdir: PathBuf,

    /// two breakpoints for gene fusion to be search(-p chr1:pos1,chr2:pos2)
    #[arg(short, long)]
    pub pos: Option<String>,

    /// bed file that has the breakpoints for gene fusions. First six columns are chr1, start, end, chr2, start, end, and starts from the seventh column is the fusion id.
    #[arg(short = 'P', long)]
    pub pos_bed: Option<PathBuf>,

    /// query two regions for possible gene fusions, return any reads that has breakpoints within the two regions. Format: chr1:start-end,chr2:start-end
    #[arg(short, long)]
    pub region: Option<String>,

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

    /// number of cached nodes for each tree in maximal
    #[arg(short = 'c', long = "cached_nodes", default_value_t = 1_000_000)]
    pub lru_size: usize,

    /// Maximum number of cached isoform chunks in memory
    #[arg(long, default_value_t = 4)]
    pub cached_chunk_number: usize,

    /// Cached isoform chunk size in Mb
    #[arg(long, default_value_t = 128)]
    pub cached_chunk_size_mb: u64,

    /// Verbose mode
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
}

impl AnnFusionCli {
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

        if self.pos.is_none()
            && self.pos_bed.is_none()
            && self.gene_gtf.is_none()
            && self.region.is_none()
        {
            error!("Please provide either --pos or --pos-bed or --region or --gene-gtf");
            is_ok = false;
        }

        if let Some(pos) = &self.pos {
            match FusionBreakPointPair::from_pos_str(pos) {
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

pub struct FusionBreakPointPair {
    pub left_chr: String,
    pub left_start: u64,
    pub left_end: u64,
    pub right_chr: String,
    pub right_start: u64,
    pub right_end: u64,
    pub idstring: String,
}

impl FusionBreakPointPair {
    pub fn from_pos_str(pos: &str) -> Result<Self> {
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

        Ok(FusionBreakPointPair {
            left_chr: breakpoints[0].0.clone(),
            left_start: breakpoints[0].1,
            left_end: breakpoints[0].1,
            right_chr: breakpoints[1].0.clone(),
            right_start: breakpoints[1].1,
            right_end: breakpoints[1].1,
            idstring: String::from("SingleQuery"),
        })
    }

    pub fn from_region_str(region: &str) -> Result<Self> {
        let parts = region.split(',').collect::<Vec<&str>>();
        if parts.len() != 2 {
            panic!("Invalid format for --region. Expected format: chr1:start-end,chr2:start-end");
        }
        let l_parts = parts[0].split(':').collect::<Vec<&str>>();
        let left_chr = l_parts[0].to_string();
        let left_chr = utils::trim_chr_prefix_to_upper(&left_chr);
        let left_start = l_parts[1].split('-').collect::<Vec<&str>>()[0]
            .parse::<u64>()
            .expect("Invalid left region start position");
        let left_end = l_parts[1].split('-').collect::<Vec<&str>>()[1]
            .parse::<u64>()
            .expect("Invalid left region end position");
        let r_parts = parts[1].split(':').collect::<Vec<&str>>();
        let right_chr = r_parts[0].to_string();
        let right_chr = utils::trim_chr_prefix_to_upper(&right_chr);
        let right_start = r_parts[1].split('-').collect::<Vec<&str>>()[0]
            .parse::<u64>()
            .expect("Invalid right region start position");
        let right_end = r_parts[1].split('-').collect::<Vec<&str>>()[1]
            .parse::<u64>()
            .expect("Invalid right region end position");
        Ok(FusionBreakPointPair {
            left_chr,
            left_start,
            left_end,
            right_chr,
            right_start,
            right_end,
            idstring: String::from("RegionQuery"),
        })
    }

    pub fn is_single_bp(&self) -> bool {
        self.left_start == self.left_end && self.right_start == self.right_end
    }

    // cacluate the query position and flank bp by region
    // for example:
    // chr1:100-200 --> chr1,150,50
    // the query position is the center position, and flank is half of the region size
    pub fn get_region_query_positions(&self) -> (String, u64, u64, String, u64, u64) {
        let left_center = (self.left_start + self.left_end) / 2;
        let left_flank = (self.left_end - self.left_start) / 2;

        let right_center = (self.right_start + self.right_end) / 2;
        let right_flank = (self.right_end - self.right_start) / 2;

        (
            self.left_chr.clone(),
            left_center,
            left_flank,
            self.right_chr.clone(),
            right_center,
            right_flank,
        )
    }
}

fn anno_single_fusion(
    breakpoints: FusionBreakPointPair,
    // mywriter: &mut MyGzWriter,
    fusionbrkpt_out: &mut TableOutput,
    cli: &AnnFusionCli,
    forest: &mut BPForest,
    archive_cache: &mut ArchiveCache,
    dbinfo: &DatasetInfo,
) -> Result<()> {
    info!(
        "Processing breakpoints: {}:{}-{}:{}",
        breakpoints.left_chr,
        breakpoints.left_start,
        breakpoints.right_chr,
        breakpoints.right_start
    );
    let left_target = forest.base_single_search_flank(
        &breakpoints.left_chr,
        breakpoints.left_start,
        cli.flank,
        cli.lru_size,
    );
    let right_target = forest.base_single_search_flank(
        &breakpoints.right_chr,
        breakpoints.right_start,
        cli.flank,
        cli.lru_size,
    );

    if left_target.is_empty() || right_target.is_empty() {
        if cli.verbose {
            error!(
                "No candidates found for breakpoints: {}:{}-{}:{}",
                breakpoints.left_chr,
                breakpoints.left_start,
                breakpoints.right_chr,
                breakpoints.right_start
            );
        }
        return Ok(());
    }

    if cli.verbose {
        info!(
            "Found {} candidates",
            left_target.len() + right_target.len()
        );
    }

    let mut fusion_evidence_vec = vec![0u32; dbinfo.get_size()];

    let mut left_isoform_count = 0;
    let mut right_isoform_count = 0;

    // process the left targets
    // check if the right portion is supported
    let unique_left = left_target.into_iter().collect::<HashSet<_>>();
    for target in unique_left {
        let merged_isoform = archive_cache.load_from_disk(&target);
        let evidence_vec = merged_isoform.check_fusion_mate_breakpoint(
            &breakpoints.right_chr,
            breakpoints.right_start,
            cli.flank,
            &dbinfo,
        );

        if evidence_vec.iter().sum::<u32>() > 0 {
            left_isoform_count += 1;
        }

        fusion_evidence_vec = fusion_evidence_vec
            .iter()
            .zip(evidence_vec.iter())
            .map(|(a, b)| a + b)
            .collect();
    }

    // process the right targets
    // check if the left portion is supported
    let unique_right = right_target.into_iter().collect::<HashSet<_>>();
    for target in unique_right {
        let merged_isoform: MergedIsoform = archive_cache.load_from_disk(&target);
        let evidence_vec = merged_isoform.check_fusion_mate_breakpoint(
            &breakpoints.left_chr,
            breakpoints.left_start,
            cli.flank,
            &dbinfo,
        );

        if evidence_vec.iter().sum::<u32>() > 0 {
            right_isoform_count += 1;
        }

        fusion_evidence_vec = fusion_evidence_vec
            .iter()
            .zip(evidence_vec.iter())
            .map(|(a, b)| a + b)
            .collect();
    }

    let mut out_line = Line::with_capacity(dbinfo.get_size());
    out_line.add_field(&breakpoints.left_chr.to_string());
    out_line.add_field(&breakpoints.left_start.to_string());
    out_line.add_field(&breakpoints.right_chr.to_string());
    out_line.add_field(&breakpoints.right_start.to_string());
    out_line.add_field(&breakpoints.idstring);
    out_line.add_field(&cli.min_read.to_string());
    out_line.add_field(&format!(
        "{}/{}",
        fusion_evidence_vec
            .iter()
            .filter(|&&x| x >= cli.min_read)
            .count()
            .to_string(),
        dbinfo.get_size()
    ));

    out_line.add_field(&left_isoform_count.to_string());
    out_line.add_field(&right_isoform_count.to_string());

    for idx in 0..dbinfo.get_size() {
        // record_parts.push(fusion_evidence_vec[idx].to_string());
        let samplechip = SampleChip::new(None, fusion_evidence_vec[idx].to_string());
        out_line.add_sample(samplechip);
    }

    fusionbrkpt_out.add_line(&mut out_line)?;

    Ok(())
}

/// Annotate a single fusion breakpoint but with each raw read deails
/// it generate a different output format as anno_single_fusion
/// breakpoints is not a single bp but two regions for each side,
/// this allows to find all combiantion for two regions, this is useful when query
/// all possible fusions withint the regions(eg. intron) from each gene
pub fn anno_single_fusion_detail(
    breakpoints: FusionBreakPointPair,
    // mywriter: &mut MyGzWriter,
    // fusionbrkpt_out: &mut TableOutput,
    cli: &AnnFusionCli,
    forest: &mut BPForest,
    archive_cache: &mut ArchiveCache,
    meta: &Meta,
    dbinfo: &DatasetInfo,
) -> Result<()> {
    // const FUSION_BRKPT_FORMAT_STR: &str = "COUNT";

    /*

        pub fusion_hash: u64,
        pub left_chrom: String,
        pub left_start: u64,
        pub left_end: u64,
        pub right_chrom: String,
        pub right_start: u64,
        pub right_end: u64,
        pub left_exons_junctions: Vec<u64>, // exon positions on the left side not splice junctions
        pub right_exons_junctions: Vec<u64>, // exon positions on the right side not splice junctions
        pub sample_name: String,
    */

    let mut out_header = Header::new();
    out_header.add_column(&"chr1")?;
    out_header.add_column(&"start1")?;
    out_header.add_column(&"end1")?;
    out_header.add_column(&"chr2")?;
    out_header.add_column(&"start2")?;
    out_header.add_column(&"end2")?;
    out_header.add_column(&"main_exon_count1")?;
    out_header.add_column(&"supp_segment_count2")?;
    out_header.add_column(&"query_part")?;
    out_header.add_column(&"main_isoforms")?;
    out_header.add_column(&"supp_aln_regions")?;
    out_header.add_column(&"sample_name")?;

    // let sample_name = dbinfo.get_sample_names();
    // for name in sample_name {
    //     out_header.add_sample_name(&name)?;
    // }

    let mut dbinfo2 = DBInfos::new();
    for (name, evidence) in dbinfo.get_sample_evidence_pair_vec() {
        dbinfo2.add_sample_evidence(&name, evidence);
    }

    let mut fusionbrkpt_out = TableOutput::new(
        cli.output.clone(),
        out_header,
        dbinfo2,
        meta.clone(),
        "".to_string(),
    );

    // sart to query the breakpoints
    // the query is different to the exact breakpoints

    info!(
        "Searching fusion in regions: {}:{}-{},{}:{}-{}",
        breakpoints.left_chr,
        breakpoints.left_start,
        breakpoints.left_end,
        breakpoints.right_chr,
        breakpoints.right_start,
        breakpoints.right_end
    );

    let (q_l_chr, q_l_pos, q_l_flank, q_r_chr, q_r_pos, q_r_flank) =
        breakpoints.get_region_query_positions();

    if cli.verbose {
        info!(
            "Query left region: {}:{} +/- {}, right region: {}:{} +/- {}",
            q_l_chr, q_l_pos, q_l_flank, q_r_chr, q_r_pos, q_r_flank
        );
    }

    let left_target = forest.base_single_search_flank(&q_l_chr, q_l_pos, q_l_flank, cli.lru_size);
    let right_target = forest.base_single_search_flank(&q_r_chr, q_r_pos, q_r_flank, cli.lru_size);

    if left_target.is_empty() || right_target.is_empty() {
        if cli.verbose {
            error!(
                "No candidates found for regions: {}:{}-{},{}:{}-{}",
                q_l_chr,
                q_l_pos - q_l_flank,
                q_l_pos + q_l_flank,
                q_r_chr,
                q_r_pos - q_r_flank,
                q_r_pos + q_r_flank
            );
        }
        return Ok(());
    }

    // check the mate support for each side
    if cli.verbose {
        info!(
            "Found {} candidates",
            left_target.len() + right_target.len()
        );
    }

    let mut total_fusion_records = Vec::new();

    // process the left targets
    // check if the right portion is supported
    let unique_left = left_target.into_iter().collect::<HashSet<_>>();
    for target in unique_left {
        let merged_isoform = archive_cache.load_from_disk(&target);

        let fusion_read_records = merged_isoform.cast_to_fusion_records(
            &q_r_chr,
            q_r_pos,
            q_r_flank,
            dbinfo,
            &breakpoints,
            "left".to_string(),
        );
        total_fusion_records.extend(fusion_read_records);
    }

    for target_right in right_target {
        let merged_isoform: MergedIsoform = archive_cache.load_from_disk(&target_right);

        let fusion_read_records = merged_isoform.cast_to_fusion_records(
            &q_l_chr,
            q_l_pos,
            q_l_flank,
            dbinfo,
            &breakpoints,
            "right".to_string(),
        );
        total_fusion_records.extend(fusion_read_records);
    }

    info!(
        "Total {} fusion read records found for the provided regions.",
        total_fusion_records.len()
    );

    // generate output for each fusion read record
    for record in total_fusion_records {
        fusionbrkpt_out.write_bytes(&record.get_string().as_bytes())?;
    }
    fusionbrkpt_out.finish()?;

    Ok(())
}

fn parse_bed2(line: &str) -> Result<FusionBreakPointPair> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Err(anyhow!(
            "Invalid bed line: {}. Expected at least 5 fields.",
            line
        ));
    }
    // trim the chr prefix and convert to uppercase
    let left_chr = fields[0].to_string();
    let left_chr = utils::trim_chr_prefix_to_upper(&left_chr);
    let left_start = fields[1]
        .parse::<u64>()
        .map_err(|_| anyhow!("Invalid position"))?;
    let right_chr = fields[2].to_string();
    let right_chr = utils::trim_chr_prefix_to_upper(&right_chr);
    let right_start = fields[3]
        .parse::<u64>()
        .map_err(|_| anyhow!("Invalid position"))?;
    let id = if fields.len() >= 5 {
        fields[4].to_string()
    } else {
        "NA".to_string()
    };

    Ok(FusionBreakPointPair {
        left_chr,
        left_start,
        left_end: left_start,
        right_chr,
        right_start,
        right_end: right_start,
        idstring: id.clone(),
    })
}

pub fn run_anno_fusion(cli: &AnnFusionCli) -> Result<()> {
    greetings2(&cli);
    cli.validate();

    let mut forest = BPForest::init(&cli.idxdir);
    let dataset_info = DatasetInfo::load_from_file(&cli.idxdir.join(DATASET_INFO_FILE_NAME))?;

    let meta = Meta::parse(&cli.idxdir.join(META_FILE_NAME), None)?;

    let mut archive_cache = ArchiveCache::new(
        cli.idxdir.clone().join(MERGED_FILE_NAME),
        cli.cached_chunk_size_mb * 1024 * 1024, // 512MB chunk size
        cli.cached_chunk_number,                // max 4 chunks in cache ~2GB
    );

    if cli.pos.is_some() || cli.pos_bed.is_some() {
        const FUSION_BRKPT_FORMAT_STR: &str = "COUNT";

        let mut out_header = Header::new();
        out_header.add_column(&"chr1")?;
        out_header.add_column(&"pos1")?;
        out_header.add_column(&"chr2")?;
        out_header.add_column(&"pos2")?;
        out_header.add_column(&"id")?;
        out_header.add_column(&"min_read")?;
        out_header.add_column(&"positive/sample_size")?;
        out_header.add_column(&"left_isoforms")?;
        out_header.add_column(&"right_isoforms")?;
        let sample_name = dataset_info.get_sample_names();
        for name in sample_name {
            out_header.add_sample_name(&name)?;
        }

        let mut dbinfo = DBInfos::new();
        for (name, evidence) in dataset_info.get_sample_evidence_pair_vec() {
            dbinfo.add_sample_evidence(&name, evidence);
        }

        let mut fusionbrkpt_out = TableOutput::new(
            cli.output.clone(),
            out_header,
            dbinfo,
            meta.clone(),
            FUSION_BRKPT_FORMAT_STR.to_string(),
        );

        if cli.pos.is_some() {
            let pos = &cli.pos.clone().unwrap();
            // let breakpoints: BreakPointPair =  take_from_cli(pos)?;

            let breakpoints = FusionBreakPointPair::from_pos_str(pos)?;

            anno_single_fusion(
                breakpoints,
                // &mut mywriter,
                &mut fusionbrkpt_out,
                &cli,
                &mut forest,
                &mut archive_cache,
                &dataset_info,
            )?;
            // mywriter.finish()?;
        } else if cli.pos_bed.is_some() {
            // open the bed file
            let pos_bed = cli.pos_bed.clone().unwrap();
            let bed_file = File::open(&pos_bed).context("Failed to open the provided bed file")?;
            let reader = std::io::BufReader::new(bed_file);

            // not use the common writer since the header is different in gtf mode(discovery mode)

            for line in reader.lines() {
                let line = line.context("Failed to read line from bed file")?;
                let breakpoints: FusionBreakPointPair = parse_bed2(&line)?;

                anno_single_fusion(
                    breakpoints,
                    // &mut mywriter,
                    &mut fusionbrkpt_out,
                    &cli,
                    &mut forest,
                    &mut archive_cache,
                    &dataset_info,
                )?;
            }
        }
        info!("Total processed samples: {}", dataset_info.get_size());
        fusionbrkpt_out.finish()?;
    } else if cli.gene_gtf.is_some() {
        const FUGEN_GENE_REGION_FORMAT_STR: &str = "COUNT";

        let mut out_header = Header::new();
        out_header.add_column(&"geneA_name")?;
        out_header.add_column(&"geneB_name")?;
        out_header.add_column(&"total_evidences")?;
        out_header.add_column(&"total_samples")?;
        out_header.add_column(&"is_two_strand")?;
        out_header.add_column(&"AtoB_primary_start")?;
        out_header.add_column(&"AtoB_primary_end")?;
        out_header.add_column(&"AtoB_supp_start")?;
        out_header.add_column(&"AtoB_supp_end")?;
        out_header.add_column(&"BtoA_primary_start")?;
        out_header.add_column(&"BtoA_primary_end")?;
        out_header.add_column(&"BtoA_supp_start")?;
        out_header.add_column(&"BtoA_supp_end")?;

        let mut dbinfo = DBInfos::new();
        for (name, evidence) in dataset_info.get_sample_evidence_pair_vec() {
            dbinfo.add_sample_evidence(&name, evidence);
            out_header.add_sample_name(&name)?;
        }

        let mut fusiondiscovery_out = TableOutput::new(
            cli.output.clone(),
            out_header,
            dbinfo,
            meta.clone(),
            FUGEN_GENE_REGION_FORMAT_STR.to_string(),
        );

        let gtf_path = cli.gene_gtf.clone().unwrap();
        let file = std::fs::File::open(gtf_path).expect("Failed to open GTF file");
        let reader = std::io::BufReader::new(file);
        let mut gtf_reader = gtfReader::new(reader);

        let mut skipped_genes = 0;

        let gene_indexing =
            GeneIntervalTree::new(&mut gtf_reader).expect("Failed to create GeneIntervalTree");
        info!("loaded {} genes", gene_indexing.count);

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

                let targets = forest.fusion_search(&quried_positions, cli.flank, 1, cli.lru_size);

                if targets.is_empty() {
                    skipped_genes += 1;
                    continue;
                }

                for rec_ptr in targets {
                    let isoform: MergedIsoform = archive_cache.load_from_disk(&rec_ptr);

                    let candidates = isoform.to_fusion_candidates();

                    if candidates.is_none() {
                        continue;
                    }

                    let candidates: Vec<FusionAggrReads> = candidates.unwrap();

                    for mut candidate in candidates {
                        if candidate.match_gene(&gene_indexing, cli.flank) {
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
                fusion_cluster.generate_record_line(&dataset_info, &mut fusiondiscovery_out)?;
            }
        }

        info!("Total processed genes: {}", gene_indexing.count);
        info!("Total skipped genes: {}", skipped_genes);
        info!("Total processed samples: {}", dataset_info.get_size());
        fusiondiscovery_out.finish()?;
    } else if cli.region.is_some() {
        let pos = &cli.region.clone().unwrap();
        let breakpoints = FusionBreakPointPair::from_region_str(pos)?;

        anno_single_fusion_detail(
            breakpoints,
            // &mut mywriter,
            // &mut fusionbrkpt_out,
            &cli,
            &mut forest,
            &mut archive_cache,
            &meta,
            &dataset_info,
        )?;
    }

    info!("Finished!");
    Ok(())
}
