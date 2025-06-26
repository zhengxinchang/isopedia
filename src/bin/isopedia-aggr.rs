use std::{cmp::Reverse, collections::BinaryHeap, env, fmt::Display, io::Write, path::PathBuf};

use clap::Parser;
use isopedia::{
    chromosome::ChromMapping,
    constants::*,
    isoform::MergedIsoform,
    isoformarchive::IsoformArchive,
    meta::Meta,
    reads::{AggrRead, SingleSampleReader},
    tmpidx::{MergedIsoformOffset, MergedIsoformOffsetPlusGenomeLoc, Tmpindex},
};
use log::{error, info};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::collections::HashSet;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia-aggr")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "
Contact: Xinchang Zheng <zhengxc93@gmail.com,Xinchang.Zheng@bcm.edu>
", long_about = None)]
#[clap(after_long_help = "Exmaple:

The input file example(tab-separated)
(reqired)            (reqired)  (optional) (optional) ...
-------------------- ---------- ---------- ----------
path                 name       meta1      meta2     <- header
/path/to/sample1.bam sample1    value1     value2
/path/to/sample2.bam sample2    value1     value2

Run aggregation:

stix-iso aggr --input input.tsv --outdir index_dir

")]
struct Cli {
    #[arg(
        short,
        long,
        help = "Input manifest(tabs-separated) file.\
                             \nFrist 2 cols are required: sample path and sample name.\
                             \nThe rest of cols are optional and will be used as sample meta.\
                             \nFirst row is header."
    )]
    pub input: PathBuf,

    /// Output index directory
    #[arg(short, long)]
    pub outdir: PathBuf,
}

impl Cli {
    fn validate(&self) {
        let mut is_ok = true;

        if !self.input.exists() {
            error!(
                "--input: input file {} does not exists.",
                self.input.display()
            );
            is_ok = false;
        } else {
            // check the format of the input file
            let content = std::fs::read_to_string(&self.input).unwrap();
            let header = content.lines().next().unwrap();
            let field_size: usize = header.split('\t').count();
            let mut uniq_path = HashSet::new();
            let mut uniq_name = HashSet::new();
            for (idx, line) in content.lines().enumerate() {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 2 {
                    error!("--input: line {}: at least 2 columns required.", idx + 1);
                    is_ok = false;
                }

                if fields.len() != field_size {
                    error!(
                        "--input: line {}: field size does not match header.",
                        idx + 1
                    );
                    is_ok = false;
                }

                if idx > 0 {
                    if !PathBuf::from(fields[0]).exists() {
                        error!(
                            "--input: line {}: file {} does not exist in the manifest file.",
                            idx + 1,
                            fields[0]
                        );
                        is_ok = false;
                    }

                    if uniq_name.contains(fields[1]) {
                        error!(
                            "--input: line {}: sample name {} is duplicated.",
                            idx + 1,
                            fields[1]
                        );
                        is_ok = false;
                    } else {
                        uniq_name.insert(fields[1]);
                    }

                    if uniq_path.contains(fields[0]) {
                        error!(
                            "--input: line {}: sample path {} is duplicated.",
                            idx + 1,
                            fields[0]
                        );
                        is_ok = false;
                    } else {
                        uniq_path.insert(fields[0]);
                    }
                }
            }
        }

        let output_dir = self.outdir.parent().unwrap();

        if !output_dir.exists() {
            error!(
                "--output: parent dir {} does not exist",
                output_dir.display()
            );
            is_ok = false;
        }

        if !self.outdir.exists() {
            // create the directory
            std::fs::create_dir_all(&self.outdir).expect("Can not create output directory");
        }

        if !self.outdir.is_dir() {
            error!("output path is not a directory")
        }

        if is_ok != true {
            panic!("Invalid arguments, please check the error messages above.");
        }
    }
}

fn greetings(args: &Cli) {
    println!("\nIsopedia: [Aggregate multiple samples]\n");
    match serde_json::to_string_pretty(&args) {
        Ok(json) => println!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
}

#[derive(Debug)]
struct HeapItem<T> {
    rec: T,
    signature: u64,
    file_idx: usize,
}

impl PartialEq for HeapItem<AggrRead> {
    fn eq(&self, other: &Self) -> bool {
        self.signature == other.signature
    }
}

impl Eq for HeapItem<AggrRead> {}

impl PartialOrd for HeapItem<AggrRead> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.signature.cmp(&other.signature))
    }
}

impl Ord for HeapItem<AggrRead> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.signature.cmp(&other.signature)
    }
}

impl Display for HeapItem<AggrRead> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "AggrIsoform{}", self)
    }
}

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let cli = Cli::parse();
    cli.validate();
    greetings(&cli);

    let sample_meta = Meta::parse(&cli.input);

    let file_list = sample_meta.get_path_list();
    if file_list.len() > MAX_SAMPLE_SIZE {
        error!(
            "The number of samples is larger than the maximum number of samples allowed: {}",
            MAX_SAMPLE_SIZE
        );
        std::process::exit(1);
    }

    let mut chroms = ChromMapping::new();
    info!("Reading files...");

    // init the file reader for each sample
    let mut file_readers: Vec<SingleSampleReader> = file_list
        .iter()
        .map(|sample| {
            let mut file_reader = SingleSampleReader::new(sample.to_str().unwrap());
            chroms.udpate_from_string(&file_reader.next_line_str().expect(
                format!("Cannot read chromsomes from file: {}", sample.display()).as_str(),
            ));
            file_reader.skip_lines(1); // skip the header
            file_reader
        })
        .collect();

    info!("Merging...");

    // new the heap
    let mut heap = BinaryHeap::new();
    // init the merge buffer
    let mut merged_map: FxHashMap<u64, MergedIsoform> = FxHashMap::default();
    // init the output file

    // let mut merged_isoform_writer = std::io::BufWriter::new(
    //     std::fs::File::create(&cli.outdir.join(MERGED_FILE_NAME))
    //       .expect("Can not create merged records file...exit")  ,
    // );

    let mut isoform_archive = IsoformArchive::create(&cli.outdir.join(MERGED_FILE_NAME));

    // let mut out_intrim_path = ag_args.outdir.clone();
    // out_intrim_path.push("interim.idx");
    let mut tmpidx = Tmpindex::create(&cli.outdir.join(TMPIDX_FILE_NAME));

    let mut merged_offset = 0;
    // init the heap
    for (idx, reader) in file_readers.iter_mut().enumerate() {
        if let Some(rec) = reader.next_rec() {
            let sig = rec.signature;
            heap.push(Reverse(HeapItem {
                rec: rec,
                signature: sig,
                file_idx: idx,
            }));
        }
    }

    let mut need_deleted = Vec::new();
    // let mut buf = Vec::new();

    // let mut processed = 0;
    let mut curr_batch = 0;
    let batch_size = 1000000;
    let mut batches = 0;
    while let Some(Reverse(HeapItem {
        rec,
        file_idx,
        signature,
    })) = heap.pop()
    {
        curr_batch += 1;

        if let Some(merged_isoform) = merged_map.get_mut(&signature) {
            merged_isoform.add(rec, file_idx as u32);
        } else {
            let merged_isoform = MergedIsoform::init(rec, sample_meta.get_size(), file_idx as u32);
            merged_map.insert(signature, merged_isoform);
        }
        if let Some(rec) = file_readers[file_idx].next_rec() {
            let _sig = rec.signature;
            heap.push(Reverse(HeapItem {
                rec: rec,
                signature: _sig,
                file_idx: file_idx,
            }));
        } else {
            continue;
        }

        // dump the buffer to disk
        if curr_batch >= batch_size {
            for (_, merged_isoform_rec) in merged_map.iter() {
                //ignore the current aggr_record which might not be the final one
                if merged_isoform_rec.signature == signature {
                    continue;
                }
                // calculate the bytes length of the merged_isoform record
                let bytes_len = isoform_archive.dump_to_disk(&merged_isoform_rec);

                need_deleted.push(merged_isoform_rec.signature);
                // buf.clear();

                let sjs = merged_isoform_rec.get_common_splice_sites();
                sjs.iter().for_each(|sj| {
                    let offset_plus_genomeloc = MergedIsoformOffsetPlusGenomeLoc {
                        chrom_id: chroms
                            .get_chrom_idx(merged_isoform_rec.chrom.as_str())
                            .unwrap(),
                        pos: *sj,
                        record_ptr: MergedIsoformOffset {
                            offset: merged_offset,
                            length: bytes_len,
                            n_splice_sites: sjs.len() as u32,
                        },
                    };

                    tmpidx.add_one(offset_plus_genomeloc);
                });

                // add the left and right position of a single read
                let read_ref_spans = merged_isoform_rec.get_read_ref_span_vec();
                read_ref_spans.iter().for_each(|position| {
                    let offset_plus_genomeloc: MergedIsoformOffsetPlusGenomeLoc =
                        MergedIsoformOffsetPlusGenomeLoc {
                            chrom_id: chroms
                                .get_chrom_idx(merged_isoform_rec.chrom.as_str())
                                .unwrap(),
                            pos: *position,
                            record_ptr: MergedIsoformOffset {
                                offset: merged_offset,
                                length: bytes_len,
                                n_splice_sites: 0,
                            },
                        };
                    tmpidx.add_one(offset_plus_genomeloc);
                });

                merged_offset += bytes_len as u64;
            }

            batches += 1;
            info!("Processed: {} records", batch_size * batches);
            curr_batch = 0;
            for sig in need_deleted.iter() {
                merged_map.remove(sig);
            }

            need_deleted.clear();
        }
    }

    for (_, merged_isoform_rec) in &merged_map {
        // let bytes_len = merged_isoform_rec.gz_encode(&mut buf);
        // merged_isoform_writer.write_all(&buf).unwrap();
        let bytes_len = isoform_archive.dump_to_disk(&merged_isoform_rec);
        // buf.clear();
        let sjs = merged_isoform_rec.get_common_splice_sites();
        sjs.iter().for_each(|sj| {
            let interim_record = MergedIsoformOffsetPlusGenomeLoc {
                chrom_id: chroms
                    .get_chrom_idx(merged_isoform_rec.chrom.as_str())
                    .unwrap(),
                pos: *sj,
                record_ptr: MergedIsoformOffset {
                    offset: merged_offset,
                    length: bytes_len,
                    n_splice_sites: sjs.len() as u32,
                },
            };
            tmpidx.add_one(interim_record);
        });


        // add the left and right position of a single read
        let read_ref_spans = merged_isoform_rec.get_read_ref_span_vec();
        read_ref_spans.iter().for_each(|position| {
            let offset_plus_genomeloc: MergedIsoformOffsetPlusGenomeLoc =
                MergedIsoformOffsetPlusGenomeLoc {
                    chrom_id: chroms
                        .get_chrom_idx(merged_isoform_rec.chrom.as_str())
                        .unwrap(),
                    pos: *position,
                    record_ptr: MergedIsoformOffset {
                        offset: merged_offset,
                        length: bytes_len,
                        n_splice_sites: 0,
                    },
                };
            tmpidx.add_one(offset_plus_genomeloc);
        });

        merged_offset += bytes_len as u64;
    }
    info!(
        "Processed: {} records",
        (batch_size * batches + 1) + curr_batch
    );

    info!(
        "Dump merged records to disk: {}",
        &cli.outdir.join(MERGED_FILE_NAME).display()
    );

    info!(
        "Dump interim index to disk: {}",
        &cli.outdir.join(TMPIDX_FILE_NAME).display()
    );

    tmpidx.dump_to_disk();

    // write the chromsome map file
    // let mut out_chrom_map_path = ag_args.outdir.clone();
    // out_chrom_map_path.push("chrom.map");
    let mut out_chrom_map_writer = std::io::BufWriter::new(
        std::fs::File::create(&cli.outdir.join(CHROM_FILE_NAME))
            .expect("Can not create chromsome map file...exit"),
    );
    out_chrom_map_writer.write_all(&chroms.encode()).unwrap();

    // write the sample_meta file
    // let mut out_sample_meta_path = ag_args.outdir.clone();
    // out_sample_meta_path.push("sample.meta");
    let mut out_sample_meta_writer = std::io::BufWriter::new(
        std::fs::File::create(&cli.outdir.join(META_FILE_NAME))
            .expect("Can not create sample meta file...exit"),
    );
    out_sample_meta_writer
        .write_all(sample_meta.get_string().as_bytes())
        .expect("can not write meta data");

    info!("Fnished");
}
