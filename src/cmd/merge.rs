use std::{cmp::Reverse, collections::BinaryHeap, fmt::Display, io::Write, path::PathBuf};

use crate::{
    chromosome::{ChromMapping, ChromMappingHelper},
    constants::*,
    dataset_info::DatasetInfo,
    meta::Meta,
    myio::GeneralOutputIO,
    pnir::PNIR,
    pnir_archive::PNIRArchiveWriter,
    reads::{AggrRead, SingleSampleReader},
    tmpidx::{MergedIsoformOffsetPlusGenomeLoc, PNIROffsetPtr, TmpIndex},
    utils::greetings2,
};
use anyhow::Result;
use clap::Parser;
use log::{error, info};
use num_format::{Locale, ToFormattedString};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

#[derive(Parser, Debug, Serialize)]
#[command(name = "isopedia merge")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(about = "
[build index, step2] Aggregate multiple samples(isoform signal file) into one index(folder).
", long_about = None)]
#[clap(after_long_help = "Example:

The manifest file example(tab-separated)
(required)  (required)                   (optional) (optional) ....
----------  ---------------------------  ---------- ---------  ----
name        path                          meta1      meta2     ....
sample1     /path/to/sample1.isoform.gz    A         15        ....
sample2     /path/to/sample2.isoform.gz    B         37        ....

Run aggregation:

isopedia merge -i manifest.tsv -o index_dir

Note that if the optional metadata columns are present, it will be used in the indexing step directly.

")]
pub struct MergeCli {
    #[arg(
        short,
        long,
        help = "Input manifest(tabs-separated) file.\
                             \nFrist 2 cols are required: sample path and sample name.
                             \nFirst row is header."
    )]
    pub input: PathBuf,

    /// Chunk size for merging
    #[arg(
        short = 'c',
        long,
        default_value_t = 1_000_000,
        help = "Chunk size for merging isoform records.\nLarger chunk size requires more memory but may speed up the merging process\nReduce the chunk size if you encounter out-of-memory error."
    )]
    pub chunk_size: usize,

    /// Output index directory
    #[arg(short, long)]
    pub outdir: PathBuf,
}

impl MergeCli {
    fn validate(&self) {
        let mut is_ok = true;

        if !self.input.exists() {
            error!(
                "--input: input file {} does not exists.",
                self.input.display()
            );
            is_ok = false;
        } else {
            let meta =
                Meta::from_file(&self.input, None, None).expect("Cannot parse the manifest file");
            meta.validate_sample_names();
            meta.validate_path_column();
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
        write!(
            f,
            "AggrIsoform(signature={}, file_idx={})",
            self.signature, self.file_idx
        )
    }
}

pub fn run_merge(cli: &MergeCli) -> Result<()> {
    // env::set_var("RUST_LOG", "info");
    // env_logger::init();

    greetings2(&cli);
    cli.validate();

    let mut dataset_info = DatasetInfo::parse_manifest(&cli.input)?;

    let file_list = dataset_info.get_path_list();

    let mut chroms = ChromMapping::new();
    info!("Reading {} files...", file_list.len());

    // init the file reader for each sample
    let mut file_readers: Vec<SingleSampleReader> = file_list
        .iter()
        .map(|sample| {
            let mut file_reader = SingleSampleReader::new(sample.to_str().unwrap());
            chroms.update_from_string(&file_reader.next_line_str().expect(
                format!("Cannot read chromosomes from file: {}", sample.display()).as_str(),
            ));
            file_reader.skip_lines(1); // skip the header
            file_reader
        })
        .collect();

    let mut chrom_helper = ChromMappingHelper::new();

    info!("Merging...");

    // new the heap
    let mut heap: BinaryHeap<Reverse<HeapItem<AggrRead>>> = BinaryHeap::new();

    // init the merge buffer
    let mut merged_map: FxHashMap<u64, PNIR> = FxHashMap::default();

    let mut tmpidx = TmpIndex::create(&cli.outdir.join(TMPIDX_FILE_NAME), cli.chunk_size);

    let isoform_archive_base = &cli.outdir.join(MERGED_FILE_NAME);

    // init the heap
    for (idx, reader) in file_readers.iter_mut().enumerate() {
        if let Some(rec) = reader.next_rec() {
            let sig = rec.signature;
            // add the evidence of the first record in each sample.
            dataset_info.add_sample_evidence(idx, rec.evidence);
            heap.push(Reverse(HeapItem {
                rec: rec,
                signature: sig,
                file_idx: idx,
            }));
        }
    }

    let mut curr_chunk_size = 0;
    let mut chunks = 0;
    let mut tmp_vec = Vec::new();
    while let Some(Reverse(HeapItem {
        rec,
        file_idx,
        signature,
    })) = heap.pop()
    {
        curr_chunk_size += 1;

        if let Some(merged_isoform) = merged_map.get_mut(&signature) {
            merged_isoform.add(rec, file_idx as u32);
        } else {
            let new_merged_isoform = PNIR::init(
                &rec,
                dataset_info.get_size(),
                file_idx as u32,
                chroms.get_chrom_idx(&rec.chrom.as_str()).unwrap(),
            );
            merged_map.insert(signature, new_merged_isoform);
        }

        if let Some(rec) = file_readers[file_idx].next_rec() {
            // add the evidence of the new record in each sample.
            dataset_info.add_sample_evidence(file_idx, rec.evidence);
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
        if curr_chunk_size >= cli.chunk_size {
            info!("Processing chunk {}...", chunks);

            // First collect records to process and signatures to remove
            let mut signatures_to_remove = Vec::new();

            // First pass: identify records to process and to remove
            for (key, merged_isoform_rec) in &merged_map {
                // Skip the current record which might not be final
                if merged_isoform_rec.signature == signature {
                    continue;
                }

                // Add to processing collection

                chrom_helper.add_record_to_chrom(merged_isoform_rec.chrom_id);

                tmp_vec.push((
                    merged_isoform_rec.chrom_id,
                    merged_isoform_rec.get_start_pos(),
                    merged_isoform_rec.signature,
                ));

                // Mark for removal
                signatures_to_remove.push(*key);
            }

            // Sort the tmp_vec by the first position
            info!(
                "Sorting {} offsets...",
                tmp_vec.len().to_formatted_string(&Locale::en)
            );
            tmp_vec.sort_by(|a, b| {
                if a.0 == b.0 {
                    a.1.cmp(&b.1)
                } else {
                    a.0.cmp(&b.0)
                }
            });
            info!(
                "Dumping {} offsets...",
                tmp_vec.len().to_formatted_string(&Locale::en)
            );

            // init the output file

            let mut isoform_archive_writer = PNIRArchiveWriter::create(
                &isoform_archive_base.with_extension(format!("chunk{}", chunks)),
            );

            // dump the sorted records to disk
            let mut merged_offset = 0;
            for (_, _, signature) in tmp_vec.iter() {
                let merged_isoform_rec = merged_map.get(signature).unwrap();
                let bytes_len = isoform_archive_writer.dump_to_disk(&merged_isoform_rec);

                let sjs = merged_isoform_rec.get_common_splice_sites();
                let mut seen_sj = FxHashSet::default();
                for sj in sjs.iter() {
                    if seen_sj.insert(sj) {
                        let interim_record = MergedIsoformOffsetPlusGenomeLoc {
                            chrom_id: merged_isoform_rec.chrom_id,
                            pos: *sj,
                            record_ptr: PNIROffsetPtr {
                                offset: merged_offset,
                                length: bytes_len,
                                n_splice_sites: sjs.len() as u32,
                            },
                        };
                        tmpidx.add_one(interim_record);
                    }
                }

                let read_ref_spans = merged_isoform_rec.get_read_ref_span_vec();
                let mut seen_span = FxHashSet::default();
                for position in read_ref_spans {
                    if seen_span.insert(position) {
                        let offset_plus_genomeloc = MergedIsoformOffsetPlusGenomeLoc {
                            chrom_id: merged_isoform_rec.chrom_id,
                            pos: position,
                            record_ptr: PNIROffsetPtr {
                                offset: merged_offset,
                                length: bytes_len,
                                n_splice_sites: 0,
                            },
                        };
                        tmpidx.add_one(offset_plus_genomeloc);
                    }
                }

                merged_offset += bytes_len as u64;
            }
            tmp_vec.clear();
            tmp_vec.shrink_to_fit();

            tmpidx.dump_chunk(chunks);

            chunks += 1;
            info!(
                "Processed: {} chunks, {} records",
                chunks.to_formatted_string(&Locale::en),
                (cli.chunk_size * chunks).to_formatted_string(&Locale::en)
            );
            curr_chunk_size = 0;
            for sig in &signatures_to_remove {
                merged_map.remove(sig);
            }

            if chunks % 2 == 0 {
                // info!("Processed {} batches", batches);
                merged_map.shrink_to_fit();
            }
            isoform_archive_writer.close_file()?;
        }
    }

    for (_, merged_isoform_rec) in &merged_map {
        chrom_helper.add_record_to_chrom(merged_isoform_rec.chrom_id);

        tmp_vec.push((
            merged_isoform_rec.chrom_id,
            merged_isoform_rec.get_start_pos(),
            merged_isoform_rec.signature,
        ));
    }

    // dump the sorted records to disk

    tmp_vec.sort_by(|a, b| {
        if a.0 == b.0 {
            a.1.cmp(&b.1)
        } else {
            a.0.cmp(&b.0)
        }
    });

    let mut isoform_archive_writer =
        PNIRArchiveWriter::create(&isoform_archive_base.with_extension(format!("chunk{}", chunks)));

    let mut merged_offset = 0;
    for (_, _, signature) in tmp_vec.iter() {
        let merged_isoform_rec = merged_map.get(signature).unwrap();
        let bytes_len = isoform_archive_writer.dump_to_disk(&merged_isoform_rec);
        let sjs = merged_isoform_rec.get_common_splice_sites();
        sjs.iter().for_each(|sj| {
            let interim_record = MergedIsoformOffsetPlusGenomeLoc {
                chrom_id: merged_isoform_rec.chrom_id,
                pos: *sj,
                record_ptr: PNIROffsetPtr {
                    offset: merged_offset,
                    length: bytes_len,
                    n_splice_sites: sjs.len() as u32,
                },
            };
            tmpidx.add_one(interim_record);
        });

        // add the left and right position of a single read, this ensure the detection of fusion gene.
        let read_ref_spans = merged_isoform_rec.get_read_ref_span_vec();
        read_ref_spans.iter().for_each(|position| {
            let offset_plus_genomeloc: MergedIsoformOffsetPlusGenomeLoc =
                MergedIsoformOffsetPlusGenomeLoc {
                    chrom_id: merged_isoform_rec.chrom_id,
                    pos: *position,
                    record_ptr: PNIROffsetPtr {
                        offset: merged_offset,
                        length: bytes_len,
                        n_splice_sites: 0,
                    },
                };
            tmpidx.add_one(offset_plus_genomeloc);
        });

        merged_offset += bytes_len as u64;
    }

    isoform_archive_writer.close_file()?;

    tmpidx.dump_chunk(chunks);

    info!(
        "Processed: {} chunks, {} records",
        chunks.to_formatted_string(&Locale::en),
        (cli.chunk_size * chunks + 1 + curr_chunk_size).to_formatted_string(&Locale::en)
    );

    info!(
        "Dump merged records to disk: {}",
        &cli.outdir.join(MERGED_FILE_NAME).display()
    );

    info!(
        "Dump tmpidx chunks to disk: {}",
        &cli.outdir.join(TMPIDX_FILE_NAME).display()
    );

    // make mb in two decimal points
    info!(
        "The memory size of tmpidx in MB: {:.2}",
        tmpidx.get_mem_size()
    );

    // info!("Dumping tmpidx to disk...");

    tmpidx.finalize(isoform_archive_base);

    // write the chromsome map file
    let mut out_chrom_map_writer = std::io::BufWriter::new(
        std::fs::File::create(&cli.outdir.join(CHROM_FILE_NAME))
            .expect("Can not create chromsome map file...exit"),
    );

    /*
     * Drop zero records from chromosome mapping
     */
    chrom_helper.drop_zero(&mut chroms);

    out_chrom_map_writer.write_all(&chroms.encode()).unwrap();

    info!(
        "Dump sample info to disk: {}",
        &cli.outdir.join(DATASET_INFO_FILE_NAME).display()
    );

    dataset_info.save_to_file(&cli.outdir.join(DATASET_INFO_FILE_NAME))?;

    info!("Finished");
    Ok(())
}
