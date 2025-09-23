use crate::constants::{MAGIC, ORDER};
use clap::error;
use indexmap::IndexMap;
use itertools::Itertools;
use log::info;
use log::error;
use memmap2::Mmap;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
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

    pub fn create(file_name: &PathBuf) -> Tmpindex {
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
        };

        interim_index
    }

    pub fn add_one(&mut self, interim_record: MergedIsoformOffsetPlusGenomeLoc) {
        self.meta.data_size += 1;
        self.offsets.push(interim_record);
    }

    pub fn dump_chunk(&mut self, chunk_idx: usize) {

        if chunk_idx != self.chunks {
            error!("The chunk index of isoform archive does not match the current tmpidx chunk index");
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

    /// finalize the tmpidx, k-way merge the chunk files and rewrite the record data file, update the offsets at the same time.
    ///
    /// the final output would be a sorted tmpidx file , a rewritten record data file
    ///
    /// this function is the replace for the previous dump_to_disk function and rewrite_sorted_records function
    pub fn finalize(&mut self, new_record_data_path: &PathBuf) {
        // 同步写入新的 tmpdix和 record data file 
        // 明天只需要work on这个函数就可以了
        


        info!(
            "Merging {} chunk files to final tmpidx file: {}",
            self.chunks,
            self.file_name.display()
        );

        // tmpidx writer
        let mut tmpidx_writer =
            fs::File::create(&self.file_name).expect("Can not create tmpidx file");

        tmpidx_writer
            .seek(std::io::SeekFrom::Start(8))
            .expect("Can not move the cursor to start after write interim index file..");

        let mut chrom_count_map = IndexMap::new();

        // k-way merge the chunk files
        let file_list = (0..self.chunks)
            .map(|i| self.file_name.with_extension(format!("chunk{}", i)))
            .collect::<Vec<_>>();
        let mut readers = file_list
            .iter()
            .map(|path| {
                let file = fs::File::open(path).expect("Can not open chunk file");
                let reader = BufReader::new(file);
                reader
            })
            .collect::<Vec<_>>();

        let mut total_idx_n = 0;

        // let mut processed_offsets = FxHashSet::default();
        let mut offset_mapping: FxHashMap<u64, u64> = FxHashMap::default();

        // mmap the input file
        let data_file = fs::File::open(record_data_path).expect("Can not open record data file");
        let mmap = unsafe { Mmap::map(&data_file).expect("mmap failed") };

        // prepare writer
        let mut data_writer = BufWriter::new(
            fs::File::create(new_record_data_path).expect("Can not create new record data file"),
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

            if let Some(&mapped_offset) = offset_mapping.get(&old_offset) {
                interim_rec.record_ptr.offset = mapped_offset;

                tmpidx_writer
                    .write_all(&interim_rec.to_bytes())
                    .expect("Can not write to tmpidx file");

                continue;
            }

            // first time seeing this record
            offset_mapping.insert(old_offset, new_offset);

            let length = interim_rec.record_ptr.length as usize;

            let slice = &mmap[old_offset as usize..old_offset as usize + length];

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

        info!("Total {} index records written, total {} unique records in the record data file", total_idx_n, offset_mapping.len());

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

        // // remove chunk files
        // for path in file_list.iter() {
        //     fs::remove_file(path).expect("Can not remove chunk file after merge");
        // }

        // // remove old record data file
        // fs::remove_file(record_data_path).expect("Can not remove old record data file");
    }


    pub fn load(file_name: &PathBuf) -> Tmpindex {
        let file = fs::File::open(file_name).expect("Can not open file");
        let mut reader = BufReader::new(fs::File::open(file_name).expect("Can not open file"));
        let mut interim_index = Tmpindex {
            meta_start: 0,
            meta: TmpindexMeta {
                magic: 0,
                chrom_offsets: FxHashMap::default(),
                data_size: 0,
            },
            offsets: Vec::new(),
            file: file,
            file_name: file_name.clone(),
            chunks: 0,
        };
        let mut buffer = [0u8; 8];
        reader
            .read_exact(&mut buffer)
            .expect("Can not read meata offset from interim index file..");
        interim_index.meta_start = u64::from_le_bytes(buffer);
        // dbg!(&interim_index.meta_start);
        let mut buffer = Vec::new();
        reader
            .seek(std::io::SeekFrom::Start(interim_index.meta_start))
            .expect("Can not move the cursor to start after read header of interim index file..");
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
