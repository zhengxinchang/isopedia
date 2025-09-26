use crate::constants::{BUF_SIZE_4M, BUF_SIZE_64M};
use crate::constants::{MAGIC, ORDER};
use indexmap::IndexMap;
use itertools::Itertools;
use log::error;
use log::info;
use memmap2::Mmap;
use num_format::{Locale, ToFormattedString};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::io::{self, BufRead};
use std::os::unix::fs::FileExt;
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Seek, Write},
    path::PathBuf,
};
use zerocopy::{FromBytes, IntoBytes};
use zerocopy_derive::{FromBytes, Immutable, IntoBytes};

use std::hash::Hash;
use std::hash::Hasher;

type ChromIdxStartart = u64;
type ChromIdxLength = u64;
type ChromId = u16;
use anyhow::Result;

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
    pub fn memory_usage_mb(&self) -> f64 {
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

        let mut total_idx_n = 0;

        // let mut processed_offsets = FxHashSet::default();
        let mut offset_mapping_vec: Vec<FxHashMap<u64, u64>> = Vec::new();

        for _ in 0..self.chunks {
            offset_mapping_vec.push(FxHashMap::default());
        }

        // mmap the input file

        let merged_data_file_list = (0..self.chunks)
            .map(|i| merged_data_name_base.with_extension(format!("chunk{}", i)))
            .collect::<Vec<_>>();

        let data_mmaps = merged_data_file_list
            .iter()
            .map(|path| {
                let file = fs::File::open(path).expect("Can not open chunk file");
                let mmap = unsafe { Mmap::map(&file).expect("mmap failed") };
                mmap.advise(memmap2::Advice::Sequential)
                    .expect("Can not set mmap advice");
                mmap
            })
            .collect::<Vec<_>>();

        let mut data_writer = BufWriter::with_capacity(
            BUF_SIZE_64M,
            fs::File::create(merged_data_name_base).expect("Can not create new record data file"),
        );

        let mut new_offset: u64 = 0;

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

        // process the heap
        while let Some(Reverse((interim_rec, idx))) = heap.pop() {
            chrom_count_map
                .entry(interim_rec.chrom_id)
                .and_modify(|c| *c += 1)
                .or_insert(1);

            total_idx_n += 1;

            if total_idx_n % self.chunk_size as u64 == 0 {
                info!(
                    "write {} offsets...",
                    total_idx_n.to_formatted_string(&Locale::en)
                );
            }

            if let Some(rec) = readers[idx]
                .by_ref()
                .bytes()
                .take(MergedIsoformOffsetPlusGenomeLoc::SIZE)
                .collect::<Result<Vec<u8>, _>>()
                .ok()
            {
                if rec.len() == MergedIsoformOffsetPlusGenomeLoc::SIZE {
                    let next_rec = MergedIsoformOffsetPlusGenomeLoc::from_bytes(&rec);
                    heap.push(Reverse((next_rec, idx)));
                }
            }

            let mut interim_rec = interim_rec.clone();
            let old_offset = interim_rec.record_ptr.offset;

            if let Some(&mapped_offset) = offset_mapping_vec[idx].get(&old_offset) {
                interim_rec.record_ptr.offset = mapped_offset;

                tmpidx_writer
                    .write_all(&interim_rec.to_bytes())
                    .expect("Can not write to tmpidx file");

                continue;
            }

            // first time seeing this record
            offset_mapping_vec[idx].insert(old_offset, new_offset);

            let length = interim_rec.record_ptr.length as usize;

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

        info!("Total {} index records written", total_idx_n);

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

        // remove chunk files
        for path in tmp_file_list.iter() {
            fs::remove_file(path).expect("Can not remove chunk file after merge");
        }

        // remove old record data file

        for path in merged_data_file_list.iter() {
            fs::remove_file(path).expect("Can not remove old record data chunk file");
        }

        // fs::remove_file(record_data_path).expect("Can not remove old record data file");
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

    pub fn get_blocks(&mut self, chrom_id: u16) -> Vec<Vec<MergedIsoformOffsetGroup>> {
        let (chrom_start_idx, chrom_length) = match self.meta.chrom_offsets.get(&chrom_id) {
            Some((start, length)) => (start, length),
            None => {
                // chromap.rs has all chromsome in the reference,
                // while the interim index file only has the chromsome that has records
                return Vec::new();
            }
        };

        let curr_offset: u64 =
            8u64 + (chrom_start_idx) * MergedIsoformOffsetPlusGenomeLoc::SIZE as u64;
        // dbg!(curr_offset);
        let mut reader: BufReader<&File> = BufReader::new(&self.file);

        let mut buffer = [0u8; MergedIsoformOffsetPlusGenomeLoc::SIZE];

        reader
            .seek(std::io::SeekFrom::Start(curr_offset))
            .expect("Can not move the cursor to start after read header of interim index file..");

        let mut interim_vec = Vec::new();

        for _ in 0..*chrom_length {
            reader
                .read_exact(&mut buffer)
                .expect("Can not read magic from interim index file..");
            let interim_rec = MergedIsoformOffsetPlusGenomeLoc::from_bytes(&buffer);
            interim_vec.push(interim_rec);
        }

        let groups: Vec<MergedIsoformOffsetGroup> = interim_vec
            .into_iter()
            .chunk_by(|x| (x.chrom_id, x.pos)) // group by chrom_id and pos
            .into_iter()
            .map(|(_, group)| {
                let mut aggr_intrim_rec: MergedIsoformOffsetGroup = MergedIsoformOffsetGroup::new();
                for rec in group {
                    if aggr_intrim_rec.is_emtpy() {
                        aggr_intrim_rec.update(rec.chrom_id, rec.pos, vec![rec.record_ptr]);
                    } else {
                        aggr_intrim_rec.record_ptr_vec.push(rec.record_ptr);
                    }
                }
                aggr_intrim_rec
            })
            .collect();

        groups.chunks(ORDER as usize).map(|x| x.to_vec()).collect()
    }

    pub fn groups(
        &self,
        chrom_id: u16,
    ) -> io::Result<impl Iterator<Item = MergedIsoformOffsetGroup>> {
        let (chrom_start_idx, chrom_length) = match self.meta.chrom_offsets.get(&chrom_id) {
            Some((s, l)) => (*s, *l),
            None => return Ok(Vec::new().into_iter()), // 空迭代器
        };

        let start_offset = 8 + chrom_start_idx * (MergedIsoformOffsetPlusGenomeLoc::SIZE as u64);

        let mut reader = BufReader::with_capacity(BUF_SIZE_4M, fs::File::open(&self.file_name)?);
        reader.seek(std::io::SeekFrom::Start(start_offset))?;

        // record iterator
        let rec_iter = (0..chrom_length).map(move |_| {
            let mut buf = [0u8; MergedIsoformOffsetPlusGenomeLoc::SIZE];
            let pos = reader.stream_position().unwrap();
            reader
                .get_ref()
                .read_exact_at(&mut buf, pos)
                .unwrap();
            reader.consume(MergedIsoformOffsetPlusGenomeLoc::SIZE);
            MergedIsoformOffsetPlusGenomeLoc::from_bytes(&buf)
        });

        // group_by (chrom_id,pos)
        let binding = rec_iter
                .chunk_by(|r| (r.chrom_id, r.pos));
        let group_iter =
            binding
                .into_iter()
                .map(|((chrom, pos), group)| {
                    let mut g = MergedIsoformOffsetGroup::new();
                    g.update(chrom, pos, group.map(|r| r.record_ptr).collect());
                    g
                });

        Ok(group_iter.collect::<Vec<_>>().into_iter())
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
    pub fn new() -> MergedIsoformOffsetGroup {
        MergedIsoformOffsetGroup {
            chrom_id: 0,
            pos: 0,
            record_ptr_vec: Vec::new(),
        }
    }

    pub fn is_emtpy(&self) -> bool {
        self.chrom_id == 0 && self.pos == 0
    }

    pub fn update(&mut self, chrom_id: u16, pos: u64, record_ptr_vec: Vec<MergedIsoformOffsetPtr>) {
        self.chrom_id = chrom_id;
        self.pos = pos;
        self.record_ptr_vec = record_ptr_vec;
    }
}

#[derive(Debug, Clone, IntoBytes, FromBytes, Immutable, Serialize, Deserialize)]
#[repr(C)]
pub struct MergedIsoformOffsetPtr {
    pub offset: u64,
    pub length: u32,
    pub n_splice_sites: u32,
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
    }
}

impl PartialOrd for MergedIsoformOffsetPtr {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.offset.cmp(&other.offset))
    }
}

impl Eq for MergedIsoformOffsetPtr {}

impl Ord for MergedIsoformOffsetPtr {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.offset.cmp(&other.offset)
    }
}

impl Hash for MergedIsoformOffsetPtr {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.offset.hash(state);
    }
}
