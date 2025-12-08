use crate::constants::{BUF_SIZE_4M, BUF_SIZE_64M};
use crate::constants::{MAGIC, ORDER};
use ahash::HashSet;
use indexmap::IndexMap;
use itertools::Itertools;
use log::error;
use log::info;
use memmap2::Mmap;
use num_format::{Locale, ToFormattedString};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
// use std::io::{self, BufRead};
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Seek, Write},
    path::PathBuf,
};
use zerocopy::{FromBytes, IntoBytes};
use zerocopy_derive::{FromBytes, Immutable, IntoBytes};

use lru::LruCache;
use std::hash::Hash;
use std::hash::Hasher;
use std::num::NonZeroUsize;

type ChromIdxStartart = u64;
type ChromIdxLength = u64;
type ChromId = u16;
// use anyhow::Result;

#[derive(Debug, Serialize, Deserialize)]
pub struct TmpindexMeta {
    pub magic: u64,
    pub chrom_offsets: FxHashMap<ChromId, (ChromIdxStartart, ChromIdxLength)>,
    pub data_size: u64,
}

// This is the Temporary Index structure for indexing.
// it use `MergedIsoformOffsetPlusGenomeLoc` to store the genomic location and the record offset of the isforms in the isoform archive file.
#[derive(Debug)]
pub struct Tmpindex {
    pub meta_start: u64,
    pub meta: TmpindexMeta,
    pub offsets: Vec<MergedIsoformOffsetPlusGenomeLoc>,
    pub file_name: PathBuf,
    pub file: File,
    pub chunks: usize,
    pub chunk_size: usize,
}

impl Tmpindex {
    pub fn memory_usage_bytes(&self) -> usize {
        let mut total = 0;
        total += size_of::<Self>();
        total += size_of::<TmpindexMeta>();
        total += self.meta.chrom_offsets.len() * (size_of::<ChromId>() + size_of::<(u64, u64)>());
        total += self.offsets.capacity() * size_of::<MergedIsoformOffsetPlusGenomeLoc>();
        total
    }

    /// Human readable MB
    pub fn get_mem_size(&self) -> f64 {
        self.memory_usage_bytes() as f64 / (1024.0 * 1024.0)
    }

    pub fn create(file_name: &PathBuf, chunk_size: usize) -> Tmpindex {
        let file = fs::File::create(file_name).expect("Can not create file");
        let interim_index = Tmpindex {
            meta_start: 0,
            meta: TmpindexMeta {
                magic: MAGIC,
                chrom_offsets: FxHashMap::default(),
                data_size: 0,
            },
            offsets: Vec::new(),
            file: file,
            file_name: file_name.clone(),
            chunks: 0,
            chunk_size,
        };

        interim_index
    }

    pub fn add_one(&mut self, interim_record: MergedIsoformOffsetPlusGenomeLoc) {
        self.meta.data_size += 1;
        self.offsets.push(interim_record);
    }

    pub fn dump_chunk(&mut self, chunk_idx: usize) {
        if chunk_idx != self.chunks {
            error!(
                "The chunk index of isoform archive does not match the current tmpidx chunk index"
            );
            std::process::exit(1);
        }

        self._sort_records();
        let chunk_file_name_path = self
            .file_name
            .with_extension(format!("chunk{}", self.chunks));

        let mut writer = BufWriter::new(
            fs::File::create(&chunk_file_name_path)
                .expect("Can not create chunk file for interim index file"),
        );
        let mut buffer: [u8; MergedIsoformOffsetPlusGenomeLoc::SIZE];
        for interim_rec in self.offsets.iter() {
            buffer = interim_rec.to_bytes();
            writer.write_all(&buffer).expect("Can not write to file");
        }
        writer.flush().expect("Can not flush the writer");
        self.chunks += 1;
        self.offsets.clear();
        self.offsets.shrink_to_fit();
    }

    fn _sort_records(&mut self) {
        self.offsets
            .sort_by_key(|x: &MergedIsoformOffsetPlusGenomeLoc| (x.chrom_id, x.pos));
    }

    pub fn finalize(&mut self, merged_data_name_base: &PathBuf) {
        info!("{}", &merged_data_name_base.display());

        info!(
            "Merging {} chunk files to final tmpidx file: {}",
            self.chunks,
            self.file_name.display()
        );

        let mut tmpidx_writer = BufWriter::with_capacity(
            64 * 1024 * 1024,
            fs::File::create(&self.file_name).expect("Can not create tmpidx file"),
        );
        tmpidx_writer
            .seek(std::io::SeekFrom::Start(8))
            .expect("Can not move the cursor to start after write interim index file..");

        let mut chrom_count_map = IndexMap::new();

        // k-way merge the chunk files
        let tmp_file_list = (0..self.chunks)
            .map(|i| self.file_name.with_extension(format!("chunk{}", i)))
            .collect::<Vec<_>>();
        let mut readers = tmp_file_list
            .iter()
            .map(|path| {
                let file = fs::File::open(path).expect("Can not open chunk file");
                let reader = BufReader::with_capacity(BUF_SIZE_4M, file);
                reader
            })
            .collect::<Vec<_>>();

        // let mut offset_mapping_vec: Vec<FxHashMap<u64, u64>> = Vec::new();

        // for _ in 0..self.chunks {
        //     offset_mapping_vec.push(FxHashMap::default());
        // }

        let mut offset_managers: Vec<ChunkOffsetManager> = (0..self.chunks)
            .map(|i| {
                ChunkOffsetManager::new(
                    i,
                    self.file_name.clone(),
                    8,  // 256 个桶
                    32, // 每个 chunk 最多 32 个桶在内存
                )
            })
            .collect();

        // mmap the input file

        let merged_data_file_list = (0..self.chunks)
            .map(|i| merged_data_name_base.with_extension(format!("chunk{}", i)))
            .collect::<Vec<_>>();

        // let data_mmaps = merged_data_file_list
        //     .iter()
        //     .map(|path| {
        //         // let file = fs::File::open(path).expect("Can not open chunk file");
        //         // let mmap = unsafe { Mmap::map(&file).expect("mmap failed") };
        //         // mmap.advise(memmap2::Advice::Sequential)
        //         //     .expect("Can not set mmap advice");
        //         // mmap

        //         let file = fs::File::open(path).expect("Can not open chunk file");
        //         let reader = BufReader::with_capacity(BUF_SIZE_1M, file);
        //         reader
        //     })
        //     .collect::<Vec<_>>();

        let data_mmaps = merged_data_file_list
            .iter()
            .map(|path| {
                let file = fs::File::open(path).expect("Can not open chunk file");
                let mmap = unsafe { Mmap::map(&file).expect("mmap failed") };

                mmap.advise(memmap2::Advice::Sequential)
                    .expect("Can not set mmap advice");

                // #[cfg(target_os = "linux")]
                // {
                //     mmap.advise(memmap2::Advice::WillNeed).ok();
                // }

                mmap
            })
            .collect::<Vec<_>>();

        let mut data_writer = BufWriter::with_capacity(
            BUF_SIZE_64M,
            fs::File::create(merged_data_name_base).expect("Can not create new record data file"),
        );

        let mut new_offset: u64 = 0;
        let mut total_idx_n = 0;

        let mut heap: BinaryHeap<Reverse<(MergedIsoformOffsetPlusGenomeLoc, usize)>> =
            BinaryHeap::new();

        // init heap

        for (idx, reader) in readers.iter_mut().enumerate() {
            let mut buffer = [0u8; MergedIsoformOffsetPlusGenomeLoc::SIZE];
            if let Ok(_) = reader.read_exact(&mut buffer) {
                let interim_rec = MergedIsoformOffsetPlusGenomeLoc::from_bytes(&buffer);

                heap.push(Reverse((interim_rec, idx)));
            }
        }

        // let mut buffer = vec![0u8; 10 * 1024 * 1024]; // 10 MB buffer

        // process the heap
        while let Some(Reverse((interim_rec, idx))) = heap.pop() {
            chrom_count_map
                .entry(interim_rec.chrom_id)
                .and_modify(|c| *c += 1)
                .or_insert(1);

            total_idx_n += 1;

            if total_idx_n % (self.chunk_size * 10) as u64 == 0 {
                info!(
                    "write {} offsets...",
                    total_idx_n.to_formatted_string(&Locale::en)
                );
            }

            let mut next_buf = [0u8; MergedIsoformOffsetPlusGenomeLoc::SIZE];
            match readers[idx].read_exact(&mut next_buf) {
                Ok(_) => {
                    let next_rec = MergedIsoformOffsetPlusGenomeLoc::from_bytes(&next_buf);
                    heap.push(Reverse((next_rec, idx)));
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {}
                Err(e) => {
                    error!("read next record failed for chunk {}: {}", idx, e);
                }
            }

            let mut interim_rec = interim_rec.clone();
            let old_offset = interim_rec.record_ptr.offset;

            // if let Some(&mapped_offset) = offset_mapping_vec[idx].get(&old_offset) {
            //     interim_rec.record_ptr.offset = mapped_offset;

            //     tmpidx_writer
            //         .write_all(&interim_rec.to_bytes())
            //         .expect("Can not write to tmpidx file");

            //     continue;
            // }

            if let Some(mapped_offset) = offset_managers[idx].get(old_offset) {
                interim_rec.record_ptr.offset = mapped_offset;
                tmpidx_writer.write_all(&interim_rec.to_bytes()).unwrap();
                continue;
            }

            // first time seeing this record
            // offset_mapping_vec[idx].insert(old_offset, new_offset);
            offset_managers[idx].insert(old_offset, new_offset);

            let length = interim_rec.record_ptr.length as usize;

            // if buffer.len() < length {
            //     buffer.resize(length, 0);
            // }

            // let slice = &data_mmaps[idx][old_offset as usize..old_offset as usize + length];

            let slice = &data_mmaps[idx][old_offset as usize..old_offset as usize + length];

            data_writer
                .write_all(slice)
                .expect("Can not write record to new record data file");

            interim_rec.record_ptr.offset = new_offset;

            // write to tmpidx
            tmpidx_writer
                .write_all(&interim_rec.to_bytes())
                .expect("Can not write to tmpidx file");

            new_offset += interim_rec.record_ptr.length as u64;
        }

        let mut _curr_offset = 0;

        for (chrom_id, count) in chrom_count_map.iter() {
            self.meta
                .chrom_offsets
                .insert(*chrom_id, (_curr_offset, *count));
            _curr_offset += *count;
        }

        info!("Total {} index offsets written", total_idx_n);

        if self.meta.data_size != total_idx_n {
            panic!("The total index records merged does not match the original size, consider the index is corrupted?");
        }

        let meta_bytes = bincode::serialize(&self.meta).expect("Can not serialize meta");
        tmpidx_writer
            .write_all(&meta_bytes)
            .expect("Can not write to file");

        tmpidx_writer
            .seek(std::io::SeekFrom::Start(0))
            .expect("Can not move the cursor to start after write interim index file..");

        let index_bytes = (total_idx_n as u64) * (MergedIsoformOffsetPlusGenomeLoc::SIZE as u64);
        self.meta_start = 8 + index_bytes;
        tmpidx_writer
            .write_all(&self.meta_start.to_le_bytes())
            .expect("Can not update the meta offset");

        tmpidx_writer
            .flush()
            .expect("Can not flush the tmpidx to disk");
        data_writer
            .flush()
            .expect("Can not flush the record data file to disk");

        info!("Clean up temporary files...");

        drop(data_mmaps);
        drop(offset_managers);

        // remove chunk files
        for path in tmp_file_list.iter() {
            fs::remove_file(path).expect("Can not remove chunk file after merge");
        }

        // remove old record data file

        for path in merged_data_file_list.iter() {
            fs::remove_file(path).expect("Can not remove old record data chunk file");
        }
    }

    pub fn load(file_name: &PathBuf) -> Tmpindex {
        let file = fs::File::open(file_name).expect("Can not open file");
        let mut reader = BufReader::with_capacity(BUF_SIZE_64M, file);
        let file2 = fs::File::open(file_name).expect("Can not open file");
        let mut interim_index = Tmpindex {
            meta_start: 0,
            meta: TmpindexMeta {
                magic: 0,
                chrom_offsets: FxHashMap::default(),
                data_size: 0,
            },
            offsets: Vec::new(),
            file: file2,
            file_name: file_name.clone(),
            chunks: 0,
            chunk_size: 0,
        };
        let mut buffer = [0u8; 8];
        reader
            .read_exact(&mut buffer)
            .expect("Can not read meata offset from interim index file...");
        interim_index.meta_start = u64::from_le_bytes(buffer);

        let mut buffer = Vec::new();
        reader
            .seek(std::io::SeekFrom::Start(interim_index.meta_start))
            .expect("Can not move the cursor to start after read header of interim index file...");
        reader
            .read_to_end(&mut buffer)
            .expect("Can not read the rest of the file");
        interim_index.meta = bincode::deserialize(&buffer).expect("Can not deserialize meta");
        interim_index
    }

    pub fn groups(&self, chrom_id: u16) -> Option<TmpIdxChunker> {
        let (chrom_start_idx, chrom_rec_counts) = match self.meta.chrom_offsets.get(&chrom_id) {
            Some((s, l)) => (*s, *l),
            None => return None,
        };

        let start_offset = 8 + chrom_start_idx * (MergedIsoformOffsetPlusGenomeLoc::SIZE as u64);
        let file = fs::File::open(&self.file_name).expect("Can not open file");

        Some(TmpIdxChunker::new(
            BufReader::with_capacity(BUF_SIZE_4M, file),
            ORDER as usize,
            start_offset,
            chrom_rec_counts,
        ))
    }
}

pub struct TmpIdxChunker {
    pub file: BufReader<File>,
    pub order: usize,
    pub chrom_start_pos: u64,
    pub chrom_record_counts: u64,
    pub hold_record: Option<MergedIsoformOffsetPlusGenomeLoc>,
    buffer: [u8; MergedIsoformOffsetPlusGenomeLoc::SIZE],
    records_vec: Vec<MergedIsoformOffsetPlusGenomeLoc>,
    chr_pos_set: HashSet<(u16, u64)>,
    curr_processed: u64,
}

impl TmpIdxChunker {
    pub fn new(
        mut file: BufReader<File>,
        order: usize,
        chrom_start_pos: u64,
        chrom_record_counts: u64,
    ) -> TmpIdxChunker {
        file.seek(std::io::SeekFrom::Start(chrom_start_pos))
            .expect("Can not seek to the chrom start pos");
        TmpIdxChunker {
            file,
            order,
            chrom_start_pos,
            chrom_record_counts,
            hold_record: None,
            buffer: [0u8; MergedIsoformOffsetPlusGenomeLoc::SIZE],
            records_vec: Vec::new(),
            chr_pos_set: HashSet::default(),
            curr_processed: 0,
        }
    }

    pub fn next_chunk(&mut self) -> Option<Vec<MergedIsoformOffsetGroup>> {
        if self.curr_processed == self.chrom_record_counts && self.hold_record.is_none() {
            return None;
        }

        if let Some(held) = self.hold_record.take() {
            self.chr_pos_set.insert((held.chrom_id, held.pos));
            self.records_vec.push(held);
        }

        while self.curr_processed < self.chrom_record_counts {
            if let Err(_) = self.file.read_exact(&mut self.buffer) {
                break;
            }
            // records_read_byte_len += MergedIsoformOffsetPlusGenomeLoc::SIZE as u64;
            self.curr_processed += 1;
            let record = MergedIsoformOffsetPlusGenomeLoc::from_bytes(&self.buffer);
            let key = (record.chrom_id, record.pos);
            if !self.chr_pos_set.contains(&key) && self.chr_pos_set.len() == self.order {
                // hold the last record for next chunk
                self.hold_record = Some(record);
                break;
            }
            self.records_vec.push(record);
            self.chr_pos_set.insert(key);
        }

        let grouped = self
            .records_vec
            .iter()
            .chunk_by(|r| (r.chrom_id, r.pos))
            .into_iter()
            .map(|((chrom, pos), group)| {
                let mut g = MergedIsoformOffsetGroup::new(chrom, pos);
                for r in group {
                    g.record_ptr_vec.push(r.record_ptr.clone());
                }
                g
            })
            .collect::<Vec<_>>();

        self.records_vec.clear();
        self.chr_pos_set.clear();

        self.records_vec.shrink_to_fit();
        self.chr_pos_set.shrink_to_fit();
        Some(grouped)
    }
}

impl Iterator for TmpIdxChunker {
    type Item = Vec<MergedIsoformOffsetGroup>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_chunk()
    }
}

/// this struct is used to write the interim index file
/// it only contaions one record pointer
#[derive(Debug, Clone)]
#[repr(C, align(8))]
pub struct MergedIsoformOffsetPlusGenomeLoc {
    pub pos: u64,
    pub record_ptr: MergedIsoformOffsetPtr,
    pub chrom_id: u16,
}

impl MergedIsoformOffsetPlusGenomeLoc {
    pub const SIZE: usize = size_of::<u16>() + size_of::<u64>() * 3;

    pub fn to_bytes(&self) -> [u8; Self::SIZE] {
        // dbg!(Self::SIZE);
        let mut bytes = [0u8; Self::SIZE];
        let mut pos = 0;

        bytes[pos..pos + 2].copy_from_slice(&self.chrom_id.to_le_bytes());
        pos += 2;
        bytes[pos..pos + 8].copy_from_slice(&self.pos.to_le_bytes());
        pos += 8;
        bytes[pos..pos + 8].copy_from_slice(&self.record_ptr.offset.to_le_bytes());
        pos += 8;
        bytes[pos..pos + 4].copy_from_slice(&self.record_ptr.length.to_le_bytes());
        pos += 4;
        bytes[pos..pos + 4].copy_from_slice(&self.record_ptr.n_splice_sites.to_le_bytes());
        // dbg!(&bytes);
        bytes
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut pos = 0;
        let chrom_id = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap());
        pos += 2;
        let position = u64::from_le_bytes(bytes[pos..pos + 8].try_into().unwrap());
        pos += 8;
        let offset = u64::from_le_bytes(bytes[pos..pos + 8].try_into().unwrap());
        pos += 8;
        let length = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap());
        pos += 4;
        let nsj = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap());
        MergedIsoformOffsetPlusGenomeLoc {
            chrom_id,
            pos: position,
            record_ptr: MergedIsoformOffsetPtr::new(offset, length, nsj),
        }
    }
}

impl PartialEq for MergedIsoformOffsetPlusGenomeLoc {
    fn eq(&self, other: &Self) -> bool {
        self.chrom_id == other.chrom_id && self.pos == other.pos
    }
}

impl Eq for MergedIsoformOffsetPlusGenomeLoc {}

impl PartialOrd for MergedIsoformOffsetPlusGenomeLoc {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.chrom_id == other.chrom_id {
            Some(self.pos.cmp(&other.pos))
        } else {
            Some(self.chrom_id.cmp(&other.chrom_id))
        }
    }
}

impl Ord for MergedIsoformOffsetPlusGenomeLoc {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.chrom_id == other.chrom_id {
            self.pos.cmp(&other.pos)
        } else {
            self.chrom_id.cmp(&other.chrom_id)
        }
    }
}

impl Hash for MergedIsoformOffsetPlusGenomeLoc {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chrom_id.hash(state);

        self.pos.hash(state);

        self.record_ptr.hash(state);
    }
}

/// this struct is used to read the interim index file
/// the record_ptr_vec is a vector of RecordPtr, which
/// is the offset of the record in the original file
#[derive(Debug, Clone)]
pub struct MergedIsoformOffsetGroup {
    pub chrom_id: u16,
    pub pos: u64,
    pub record_ptr_vec: Vec<MergedIsoformOffsetPtr>,
}

impl MergedIsoformOffsetGroup {
    pub fn new(chrom_id: u16, pos: u64) -> MergedIsoformOffsetGroup {
        MergedIsoformOffsetGroup {
            chrom_id,
            pos,
            record_ptr_vec: Vec::new(),
        }
    }

    pub fn add(&mut self, chrom_id: u16, pos: u64, record_ptr_vec: Vec<MergedIsoformOffsetPtr>) {
        assert!(self.chrom_id == chrom_id && self.pos == pos);
        self.record_ptr_vec = record_ptr_vec;
    }
}

#[derive(Debug, Clone, IntoBytes, FromBytes, Immutable, Serialize, Deserialize)]
#[repr(C)]
pub struct MergedIsoformOffsetPtr {
    pub offset: u64,
    pub length: u32,
    pub n_splice_sites: u32, // 0 means the offset is for read terminal positions, otherwise its for splice junction positions
}
impl MergedIsoformOffsetPtr {
    pub fn new(offset: u64, length: u32, n_splice_sites: u32) -> Self {
        Self {
            offset,
            length,
            n_splice_sites,
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.as_bytes().to_vec()
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        if bytes.len() != std::mem::size_of::<MergedIsoformOffsetPtr>() {
            panic!(
                "{}",
                format!(
                    "RecordPtr bytes should be {} bytes long, while the length is {}",
                    std::mem::size_of::<MergedIsoformOffsetPtr>(),
                    bytes.len()
                )
            );
        }
        let record_bytes: &[u8; std::mem::size_of::<MergedIsoformOffsetPtr>()] = bytes
            [0..std::mem::size_of::<MergedIsoformOffsetPtr>()]
            .try_into()
            .unwrap();
        MergedIsoformOffsetPtr::read_from_bytes(record_bytes).unwrap()
    }
}

impl PartialEq for MergedIsoformOffsetPtr {
    fn eq(&self, other: &Self) -> bool {
        self.offset == other.offset
            && self.length == other.length
            && self.n_splice_sites == other.n_splice_sites
    }
}

impl PartialOrd for MergedIsoformOffsetPtr {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // Some(self.offset.cmp(&other.offset))

        // if offsets are equal, compare length, if length equal, compare n_splice_sites
        if self.offset == other.offset {
            if self.length == other.length {
                Some(self.n_splice_sites.cmp(&other.n_splice_sites))
            } else {
                Some(self.length.cmp(&other.length))
            }
        } else {
            Some(self.offset.cmp(&other.offset))
        }
    }
}

impl Eq for MergedIsoformOffsetPtr {}

impl Ord for MergedIsoformOffsetPtr {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // if offsets are equal, compare length, if length equal, compare n_splice_sites
        if self.offset == other.offset {
            if self.length == other.length {
                self.n_splice_sites.cmp(&other.n_splice_sites)
            } else {
                self.length.cmp(&other.length)
            }
        } else {
            self.offset.cmp(&other.offset)
        }
    }
}

impl Hash for MergedIsoformOffsetPtr {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // self.offset.hash(state);
        self.offset.hash(state);
        self.length.hash(state);
        self.n_splice_sites.hash(state);
    }
}

pub struct ChunkOffsetManager {
    chunk_idx: usize,
    base_path: PathBuf,
    bucket_bits: u32,
    buckets: LruCache<u32, FxHashMap<u64, u64>>, // 直接用 LruCache
    max_buckets: usize,
}

impl ChunkOffsetManager {
    pub fn new(chunk_idx: usize, base_path: PathBuf, bucket_bits: u32, max_buckets: usize) -> Self {
        Self {
            chunk_idx,
            base_path,
            bucket_bits,
            buckets: LruCache::new(
                NonZeroUsize::new(max_buckets).expect("max_buckets must be > 0"),
            ),
            max_buckets,
        }
    }

    #[inline]
    fn bucket_index(&self, offset: u64) -> u32 {
        (offset >> (64 - self.bucket_bits)) as u32
    }

    fn bucket_path(&self, bucket_idx: u32) -> PathBuf {
        // make a file under the base path with chunk idx and bucket idx
        self.base_path
            .with_extension(format!("chunk{}_bucket{}", self.chunk_idx, bucket_idx))
    }

    /// 保存 bucket 到磁盘
    fn save_bucket(&self, bucket_idx: u32, map: &FxHashMap<u64, u64>) {
        if map.is_empty() {
            return;
        }
        let mut writer = BufWriter::new(
            File::create(self.bucket_path(bucket_idx)).expect("Failed to create bucket file"),
        );
        for (&k, &v) in map.iter() {
            writer
                .write_all(&k.to_le_bytes())
                .expect("Failed to write key");
            writer
                .write_all(&v.to_le_bytes())
                .expect("Failed to write value");
        }
    }

    /// 从磁盘加载 bucket
    fn load_bucket(&self, bucket_idx: u32) -> FxHashMap<u64, u64> {
        let path = self.bucket_path(bucket_idx);
        if !path.exists() {
            return FxHashMap::default();
        }

        let mut map = FxHashMap::default();
        let mut reader = BufReader::new(File::open(&path).expect("Failed to open bucket file"));
        let mut buf = [0u8; 16];

        while reader.read_exact(&mut buf).is_ok() {
            let k = u64::from_le_bytes(
                buf[0..8]
                    .try_into()
                    .expect("Failed to convert bytes to u64"),
            );
            let v = u64::from_le_bytes(
                buf[8..16]
                    .try_into()
                    .expect("Failed to convert bytes to u64"),
            );
            map.insert(k, v);
        }

        std::fs::remove_file(&path).ok();
        map
    }

    /// 确保 bucket 在缓存中
    fn ensure_bucket(&mut self, bucket_idx: u32) {
        if self.buckets.contains(&bucket_idx) {
            return;
        }

        // 驱逐最老的 bucket
        if self.buckets.len() >= self.max_buckets {
            if let Some((old_idx, old_map)) = self.buckets.pop_lru() {
                self.save_bucket(old_idx, &old_map);
            }
        }

        // 加载新 bucket
        let map = self.load_bucket(bucket_idx);
        self.buckets.put(bucket_idx, map);
    }

    pub fn get(&mut self, old_offset: u64) -> Option<u64> {
        let idx = self.bucket_index(old_offset);
        self.ensure_bucket(idx);
        self.buckets
            .get(&idx)
            .and_then(|m| m.get(&old_offset).copied())
    }

    pub fn insert(&mut self, old_offset: u64, new_offset: u64) {
        let idx = self.bucket_index(old_offset);
        self.ensure_bucket(idx);
        self.buckets
            .get_mut(&idx)
            .unwrap()
            .insert(old_offset, new_offset);
    }

    pub fn cleanup(&mut self) {
        // 先保存所有内存中的 bucket
        while let Some((idx, map)) = self.buckets.pop_lru() {
            self.save_bucket(idx, &map);
        }
        // 删除所有磁盘文件
        for i in 0..(1 << self.bucket_bits) {
            std::fs::remove_file(self.bucket_path(i)).ok();
        }
    }
}

impl Drop for ChunkOffsetManager {
    fn drop(&mut self) {
        self.cleanup();
    }
}
