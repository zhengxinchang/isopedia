/// B+ tree implementation for genomic data indexing
use crate::chromosome::ChromMapping;
use crate::constants::*;
pub type RangeSearchHits = (Option<(KeyType, u64, u64)>, Vec<ValueType>);
use crate::tmpidx::MergedIsoformOffsetGroup;
use crate::tmpidx::MergedIsoformOffsetPtr;
use crate::tmpidx::TmpIdxChunker;
use crate::tmpidx::Tmpindex;
use ahash::HashSet;
use anyhow::Result;
use itertools::Itertools;
use lru::LruCache;
use memmap2::Mmap;
use rustc_hash::FxHashMap;
use std;
use std::fmt::Debug;
use std::fs::OpenOptions;
use std::num::NonZeroUsize;
use std::os::unix::fs::FileExt;
use std::path::PathBuf;
use std::sync::Arc;
use std::{
    fs::File,
    io::{self, Read},
};
use zerocopy::FromBytes;
use zerocopy::IntoBytes;
use zerocopy_derive::FromBytes;
use zerocopy_derive::Immutable;
use zerocopy_derive::IntoBytes;

/// within the node header, used to store the offset and length of the isoform offsets in
/// node payload. the reason why use it beacuse node header is fixed size.
#[derive(Debug, Clone)]
pub struct LeafNodeUtils {}

impl LeafNodeUtils {
    pub fn to_combined(offset_in_payload: u32, size: u32) -> u64 {
        ((offset_in_payload as u64) << 32) | (size as u64)
    }

    pub fn from_combined(combined: u64) -> (u32, u32) {
        ((combined >> 32) as u32, combined as u32)
    }
}

#[derive(Clone, Copy, FromBytes, IntoBytes, Immutable)]
#[repr(C)]

pub struct NodeHeader {
    pub node_id: NodeIDType, //start with 1
    pub tree_level: u64,     // leaf is 0, increase by 1 for each level
    pub num_keys: u64,       // size of the node
    pub keys: [KeyType; ORDER as usize],
    pub childs: [ValueType; ORDER as usize],
    pub payload_offset: u64, // offset of the payload
    pub payload_size: u64,   // size of the payload
    pub is_leaf_flag: u8,    // 1 for leaf, 0 for internal
    _padding: [u8; 4096
        - (std::mem::size_of::<NodeIDType>()
            + std::mem::size_of::<u64>()
            + std::mem::size_of::<u64>()
            + std::mem::size_of::<[KeyType; ORDER as usize]>()
            + std::mem::size_of::<[ValueType; ORDER as usize]>()
            + std::mem::size_of::<u8>()
            + std::mem::size_of::<u64>()
            + std::mem::size_of::<u64>())],
}

impl NodeHeader {
    pub fn init(id: NodeIDType, is_leaf: bool, level: u64) -> Self {
        assert!(std::mem::size_of::<NodeHeader>() == 4096);
        Self {
            node_id: id,
            keys: [0; ORDER as usize],
            childs: [0; ORDER as usize],
            is_leaf_flag: is_leaf as u8,
            tree_level: level,
            num_keys: 0,
            // record_data_offset: NodeHeaderPayloadOffset::new(0, 0),
            payload_offset: 0, // offset of the payload
            payload_size: 0,   // size of the payload
            _padding: [0u8; 4096
                - (std::mem::size_of::<NodeIDType>()
                    + std::mem::size_of::<[KeyType; ORDER as usize]>()
                    + std::mem::size_of::<[ValueType; ORDER as usize]>()
                    + std::mem::size_of::<u8>()
                    + std::mem::size_of::<u64>()
                    + std::mem::size_of::<u64>()
                    + std::mem::size_of::<u64>()
                    + std::mem::size_of::<u64>())],
        }
    }

    pub fn is_full(&self) -> bool {
        self.num_keys == ORDER
    }

    pub fn is_leaf(&self) -> bool {
        self.is_leaf_flag == 1
    }

    pub fn set_keys(&mut self, keys: &[KeyType]) -> usize {
        let len = keys.len();
        if len as u64 > ORDER {
            panic!("keys length is larger than ORDER");
        }
        self.keys[..len].copy_from_slice(&keys);
        return len;
    }
    pub fn set_non_leaf_node_childrens(&mut self, node_ids: &[KeyType]) -> usize {
        let len = node_ids.len();
        if len as u64 > ORDER {
            panic!("keys length is larger than ORDER");
        }
        self.childs[..len].copy_from_slice(&node_ids);
        return len;
    }
}

#[derive(Debug, Clone)]
pub struct NodeData {
    num_records: usize,
    pub merge_isoform_offset_vec: Vec<MergedIsoformOffsetPtr>,
}

impl NodeData {
    pub fn new() -> Self {
        Self {
            num_records: 0,
            merge_isoform_offset_vec: Vec::new(),
        }
    }

    pub fn add_record(&mut self, archive_offset: MergedIsoformOffsetPtr) {
        self.merge_isoform_offset_vec.push(archive_offset);
        self.num_records += 1;
    }

    pub fn add_records(&mut self, merge_isoform_offsets: Vec<MergedIsoformOffsetPtr>) -> usize {
        let len = merge_isoform_offsets.len();
        self.merge_isoform_offset_vec.extend(merge_isoform_offsets);
        self.num_records += len.clone();
        return len;
    }

    pub fn get_records(&self, start: usize, length: usize) -> Vec<MergedIsoformOffsetPtr> {
        self.merge_isoform_offset_vec[start..start + length].to_vec()
    }
}

#[derive(Clone)]
pub struct Node {
    pub header: NodeHeader,
    pub data: NodeData,
}

impl Node {
    pub fn new(id: NodeIDType, is_leaf: bool, level: u64) -> Self {
        Self {
            header: NodeHeader::init(id, is_leaf, level),
            data: NodeData::new(),
        }
    }

    pub fn is_full(&self) -> bool {
        self.header.is_full()
    }

    pub fn is_leaf(&self) -> bool {
        self.header.is_leaf()
    }

    pub fn new_leaf_node_from_batch(id: NodeIDType, block: &Vec<MergedIsoformOffsetGroup>) -> Self {
        // init a leaf node
        // dbg!(&block.len());
        let mut node = Node::new(id, true, 0);
        // add RecordPtrs and update keys and childrens
        let mut curr_idx: u32 = 0;
        for (idx, interim_record) in block.into_iter().enumerate() {
            node.header.keys[idx] = interim_record.pos;
            node.header.childs[idx] =
                LeafNodeUtils::to_combined(curr_idx, interim_record.record_ptr_vec.len() as u32);
            node.add_header_size(1u64);
            node.data
                .add_records(interim_record.record_ptr_vec.to_owned()); // length updated in node.data.length;
            curr_idx += interim_record.record_ptr_vec.len() as u32;
        }
        node
    }

    pub fn get_header_bytes(&self) -> [u8; 4096] {
        let bytes: &[u8] = self.header.as_bytes();
        let mut array = [0u8; 4096];
        array.copy_from_slice(bytes);
        array
    }

    pub fn load_header_from_bytes(b: &[u8]) -> Self {
        let header: NodeHeader = NodeHeader::read_from_bytes(b).unwrap();
        Self {
            header: header,
            data: NodeData::new(),
        }
    }

    /// serialize the record pointers in the node to bytes
    pub fn serialize_payload2bytes(&self) -> Vec<u8> {
        let mut bytes: Vec<u8> = Vec::new();
        for record in &self.data.merge_isoform_offset_vec {
            bytes.extend(record.to_bytes());
        }
        bytes
    }

    /// get max_key of the node
    pub fn get_max_key(&self) -> KeyType {
        let max_key = self.header.keys[(self.header.num_keys - 1) as usize];
        max_key
    }

    pub fn get_min_key(&self) -> KeyType {
        self.header.keys[0]
    }

    pub fn set_inner_data_range(&mut self, offset: &u64, length: &u64) {
        self.header.payload_offset = *offset;
        self.header.payload_size = *length;
    }

    pub fn load_record_pointers_from_bytes(&mut self, b: &[u8]) {
        let mut offset = 0;
        let mut length = 0;
        let mut record_addr_list: Vec<MergedIsoformOffsetPtr> = Vec::new();
        while offset < b.len() {
            let record = MergedIsoformOffsetPtr::from_bytes(&b[offset..offset + 16]);
            record_addr_list.push(record);
            offset += 16;
            length += 1;
        }
        self.data.merge_isoform_offset_vec = record_addr_list;
        self.data.num_records = length;
    }

    pub fn add_header_size(&mut self, size: u64) {
        self.header.num_keys += size;
    }

    /// search one key in the node, return the value if found
    pub fn search(&self, key: KeyType) -> Option<ValueType> {
        match self.header.keys[..self.header.num_keys as usize].binary_search(&key) {
            Ok(idx) => Some(self.header.childs[idx]),
            Err(idx) => {
                // 对于 search，Err(idx) 表示 key 应该插入在 idx 处
                if idx < self.header.num_keys as usize {
                    Some(self.header.childs[idx])
                } else {
                    None
                }
            }
        }
    }

    /// search one key in the node
    pub fn exact_search(&self, key: KeyType) -> Option<ValueType> {
        match self.header.keys[..self.header.num_keys as usize].binary_search(&key) {
            Ok(idx) => Some(self.header.childs[idx]),
            Err(_) => None,
        }
    }

    // type RangeRes=(u64,u64,Vec<ValueType>);
    // seach a range of keys in the node, if the key's pos larger than the max_key,
    // return `(Option<(NodeIDType, u64, u64)>, Vec<ValueType>)
    pub fn range_search(&self, start: KeyType, end: KeyType) -> RangeSearchHits {
        let mut ret: RangeSearchHits = (None, Vec::new());

        /*
           a-----node----b
               start---------end
                         |-----| <- search range for next node

        */

        if self.header.keys[(self.header.num_keys - 1) as usize] < end {
            ret.0 = Some((
                self.header.node_id + 1,
                self.header.keys[self.header.num_keys as usize - 1],
                end,
            ));
        }
        let mut i = 0;
        while i < self.header.num_keys {
            if self.header.keys[i as usize] >= start && self.header.keys[i as usize] <= end {
                ret.1.push(self.header.childs[i as usize]);
            }
            i += 1;
        }
        ret
    }

    pub fn get_addresses_by_children_value(
        &self,
        value: &ValueType,
    ) -> Vec<MergedIsoformOffsetPtr> {
        if !self.is_leaf() {
            panic!("get_addresses_by_children_value is only for leaf node");
        }
        let (offset, size) = LeafNodeUtils::from_combined(*value);
        self.data.get_records(offset as usize, size as usize)
    }
}

#[derive(Debug, Clone, Copy, IntoBytes, FromBytes, Immutable)]
// The IntoBytes trait is affected by the rank of the fields in the struct. Large files should be placed at the start of the struct.
#[repr(C)]
pub struct CacheHeader {
    pub magic_number: u64,

    pub first_node_offset: u64,      // 8*u8 // always 4096
    pub root_node_id: u64,           // 8*u8 // this id also indicates the number of nodes
    pub payload_start_offset: usize, // the start offset of the payload area, always 4096 * (total_nodes + 1)
    pub curr_payload_offset: u64, // 8*u8 // point to the end of the last node, dont use it when load from disk
    pub total_nodes: u64,         // 8*u8 // this is the total number of nodes
    pub leaf_nodes: u64,          // 8*u8 // this is the number of leaf nodes
    pub num_of_genome_pos: u64,   // 8*u8 // total number of keys in the tree
    pub bptree_order: u64,        // 8*u8. // equals to ORDER
    pub height: u64,              // 8*u8 // height of the tree
    pub chromosome_id: u16,       // 2*u8 // chromosome id
    _padding: [u8; 4096 - (std::mem::size_of::<u64>() * 10 + std::mem::size_of::<u16>())],
}

impl CacheHeader {
    pub fn new() -> CacheHeader {
        CacheHeader {
            magic_number: MAGIC,

            first_node_offset: 0,
            chromosome_id: 0,
            root_node_id: 0,
            curr_payload_offset: 0,
            payload_start_offset: 0,
            total_nodes: 0,
            leaf_nodes: 0,
            num_of_genome_pos: 0,
            bptree_order: ORDER,
            height: 0,
            _padding: [0; (4096 - (std::mem::size_of::<u64>() * 10 + std::mem::size_of::<u16>()))],
        }
    }
    pub fn from_bytes(b: &[u8]) -> Result<CacheHeader> {
        let header: CacheHeader = CacheHeader::read_from_bytes(b).unwrap();
        Ok(header)
    }
    pub fn to_bytes(&self) -> [u8; 4096] {
        let header_bytes: [u8; 4096] = self
            .as_bytes()
            .try_into()
            .expect("slice with incorrect length");
        return header_bytes;
    }
}

// type NodeIDType = u64;
// type KeyType = u64;

#[derive(Debug)]
#[repr(C)]
pub struct Cache {
    pub header: CacheHeader,
    pub file: File,
    pub file_path: PathBuf,
    pub payload_file: Option<File>,
    pub payload_file_path: Option<PathBuf>,
    pub mmap: Option<memmap2::Mmap>,
    pub lru: LruCache<u64, Arc<Node>>,
}

impl Cache {
    pub fn create(file_path: &PathBuf) -> Self {
        let header = CacheHeader::new();
        let file = File::create(file_path).unwrap();
        let payload_file_path = file_path.with_extension("payload");
        // let payload_file = File::create(&payload_file_path).unwrap();

        let payload_file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&payload_file_path)
            .unwrap();

        Cache {
            header,
            file,
            file_path: file_path.clone(),
            payload_file: Some(payload_file),
            payload_file_path: Some(payload_file_path),
            mmap: None,
            lru: LruCache::new(NonZeroUsize::new(10).unwrap()),
        }
    }

    pub fn finish(&mut self) -> Result<()> {
        // update header and write to disk

        let header_bytes = self.header.to_bytes();
        self.file.write_all_at(&header_bytes, 0).unwrap();
        self.close();

        // copy the payload file to the main file
        let payload_file = self.payload_file.as_ref().expect("Payload file not found");
        let payload_size = payload_file.metadata().unwrap().len() as usize; // e.g., 200
        let mut payload_reader = io::BufReader::new(payload_file);
        let mut buffer = vec![0; BUF_SIZE_64M];
        let mut offset = self.header.payload_start_offset; // e.g., 1000
        let mut write_len = 0;

        while write_len < payload_size {
            let read_size = std::cmp::min(BUF_SIZE_64M, payload_size - write_len);
            payload_reader.read_exact(&mut buffer[..read_size])?;
            self.file
                .write_all_at(&buffer[..read_size], offset as u64)?;
            offset += read_size;
            write_len += read_size;
        }
        // remove the payload file
        self.payload_file = None;
        if let Some(payload_file_path) = self.payload_file_path.clone() {
            std::fs::remove_file(payload_file_path)?;
        }
        Ok(())
    }

    // pub fn update_header(&mut self, header: CacheHeader) {
    //     self.header = header;
    // }

    pub fn dump_one_node(&mut self, node: &mut Node) -> Result<()> {
        // dump the data first if the node is leaf
        assert!(std::mem::size_of_val(&node.header) == 4096); // the size of the node header should be 4096
        if node.is_leaf() {
            let payload_bytes: Vec<u8> = node.serialize_payload2bytes();
            let payload_bytes: &[u8] = payload_bytes.as_slice();
            let payload_bytes_len: usize = payload_bytes.len();
            // self.file
            //     .write_all_at(payload_bytes, self.header.curr_payload_offset)
            //     .unwrap();
            self.payload_file
                .as_ref()
                .unwrap()
                .write_all_at(payload_bytes, self.header.curr_payload_offset)?;
            // update the record_range
            node.set_inner_data_range(
                // need to be updated to a tmp offset starts from 0
                &self.header.curr_payload_offset,
                &(payload_bytes_len as u64),
            );
            // update max_offset
            self.header.curr_payload_offset += payload_bytes_len as u64;
        }
        // dump header
        let node_header_bytes: [u8; 4096] = node.get_header_bytes();
        self.file
            .write_all_at(&node_header_bytes, node.header.node_id * 4096)?;

        Ok(())
    }

    pub fn from_disk(file_path: &str, lru_size: usize) -> Cache {
        let mut file: File = File::open(file_path).unwrap();
        let mut header_bytes = [0; 4096];
        file.read_exact(&mut header_bytes).unwrap();

        let magic = u64::from_ne_bytes(header_bytes[0..8].try_into().unwrap());
        if magic != MAGIC {
            panic!(
                "{}",
                format!("The file {:?}is not a valid data cache file", file_path)
            );
        }

        let header = CacheHeader::from_bytes(&header_bytes).unwrap();
        Cache {
            header,
            file,
            file_path: PathBuf::from(file_path),
            payload_file: None,
            payload_file_path: None,
            mmap: None,
            lru: LruCache::new(NonZeroUsize::new(lru_size).unwrap()),
        }
    }

    pub fn from_disk_mmap(idx_path: &str, lru_size: usize) -> Result<Cache> {
        let file = File::options().read(true).write(true).open(idx_path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let header = CacheHeader::from_bytes(&mmap[0..4096])?;
        Ok(Cache {
            header,
            file,
            file_path: PathBuf::from(idx_path),
            payload_file: None,
            payload_file_path: None,
            mmap: Some(mmap),
            lru: LruCache::new(NonZeroUsize::new(lru_size).unwrap()),
        })
    }

    pub fn get_node2(&mut self, node_id: u64) -> Option<Arc<Node>> {
        if let Some(n) = self.lru.get(&node_id) {
            // 命中缓存，克隆 Arc 返回
            return Some(Arc::clone(n));
        }

        let mmap = self.mmap.as_ref().expect("Cache not mmap’d");
        if node_id == 0 || node_id > self.header.total_nodes {
            return None;
        }

        let off = (node_id * 4096) as usize;
        let hdr_slice = &mmap[off..off + 4096];
        let mut node = Node::load_header_from_bytes(hdr_slice);

        let po = node.header.payload_offset as usize + self.header.payload_start_offset;
        let ps = node.header.payload_size as usize;
        let data_slice = &mmap[po..po + ps];
        node.load_record_pointers_from_bytes(data_slice);

        // 放入缓存
        let arc_node = Arc::new(node);
        self.lru.put(node_id, Arc::clone(&arc_node));

        Some(arc_node)
    }

    pub fn get_root_node(&mut self) -> Arc<Node> {
        self.get_node2(self.header.root_node_id).unwrap()
    }

    pub fn get_max_key(&mut self) -> u64 {
        let max_leaf = self.get_node2(self.header.leaf_nodes);
        // dbg!(&max_leaf);
        let max_key = max_leaf.unwrap().get_max_key();
        max_key
    }

    pub fn get_min_key(&mut self) -> u64 {
        let min_leaf = self.get_node2(1);
        let min_key = min_leaf.unwrap().get_min_key();
        min_key
    }

    pub fn close(&mut self) {
        self.file.sync_all().unwrap();
        self.payload_file.as_ref().map(|f| {
            f.sync_all().unwrap();
        });
    }
}

pub struct BPTree {
    pub root: KeyType,
    pub idxdir: PathBuf,
    pub cache: Option<Cache>,
    pub size: u32,  // number of nodes
    pub hight: u16, // height of the tree;
    pub chrom_id: u16,
    pub order: u64,
}

impl BPTree {
    pub fn init(idx_path: &PathBuf, chrom_id: u16) -> BPTree {
        // let file: File = File::open(file_path).expect("can not open file");
        BPTree {
            root: 0,
            idxdir: idx_path.clone(),
            cache: None,
            size: 0,
            hight: 0,
            chrom_id: chrom_id,
            order: ORDER,
        }
    }

    pub fn from_disk(idx_path: &PathBuf, chrom_id: u16, lru_size: usize) -> BPTree {
        let cache = Cache::from_disk_mmap(
            idx_path
                .join(format!("bptree_{}.idx", chrom_id))
                .to_str()
                .unwrap(),
            lru_size,
        )
        .expect("Can not open cache file");
        BPTree {
            root: 0,
            idxdir: PathBuf::from(idx_path),
            cache: Some(cache),
            size: 0,
            hight: 0,
            chrom_id: chrom_id,
            order: ORDER,
        }
    }

    /// build the B+ tree from the interim records
    pub fn build_tree(
        // aggr_intrim_rec_vec: Vec<Vec<MergedIsoformOffsetGroup>>,
        tmpidx_chunker: TmpIdxChunker,
        idx_path: &PathBuf,
        chrom_id: u16,
    ) -> Result<()> {
        let mut tree: BPTree = BPTree {
            root: 0,
            idxdir: idx_path.clone(),
            cache: None,
            size: 0,
            hight: 0,
            chrom_id: 0,
            order: ORDER,
        };

        let mut cache = Cache::create(&idx_path.join(format!("bptree_{}.idx", chrom_id)));
        // let mut cache_header = CacheHeader::new();
        cache.header.chromosome_id = chrom_id;

        let mut node_list: Vec<(KeyType, NodeIDType)> = Vec::new();
        let mut node_id: u64 = 0;

        for block in tmpidx_chunker {
            // let block: Vec<MergedIsoformOffsetGroup> = block.collect();
            node_id += 1;
            // num_of_ptrs += block.len() as u64;
            cache.header.leaf_nodes += 1;
            cache.header.num_of_genome_pos += block.len() as u64;
            cache.header.total_nodes += 1;

            let mut node: Node = Node::new_leaf_node_from_batch(node_id.clone(), &block);
            node_list.push((node.get_max_key(), node_id));
            cache.dump_one_node(&mut node)?;
            // tree_size += 1;
        }

        let mut start_idx = 0;
        let mut level = 1;
        // let num_of_leaf_nodes = node_list.len();
        // start to construct internal nodes
        loop {
            let mut internal_node_list: Vec<(KeyType, NodeIDType)> = node_list[start_idx..]
                .chunks(ORDER as usize)
                .map(|group| {
                    node_id += 1;
                    cache.header.total_nodes += 1;

                    let mut node = Node::new(node_id, false, level);
                    let keys: Vec<u64> = group.iter().map(|(key, _)| *key).collect::<Vec<u64>>();
                    let childerns = group
                        .into_iter()
                        .map(|(_, node_id)| *node_id)
                        .collect::<Vec<NodeIDType>>();
                    let key_arr = keys.as_slice();
                    node.header.set_keys(key_arr);
                    let childerns_arr = childerns.as_slice();
                    node.header.set_non_leaf_node_childrens(childerns_arr);
                    node.header.num_keys = key_arr.len() as u64;
                    cache
                        .dump_one_node(&mut node)
                        .expect(&format!("Can not dump node {}", node_id));
                    (node.get_max_key(), node_id)
                })
                .collect();

            start_idx = node_list.len();
            node_list.extend(internal_node_list.clone());
            level += 1;
            // tree_size += internal_node_list.len();
            if internal_node_list.len() == 1 {
                tree.root = internal_node_list[0].1;
                break;
            }
            internal_node_list.clear();
        }
        tree.hight = level as u16;
        tree.size = node_list.len() as u32;
        tree.chrom_id = chrom_id;
        tree.order = ORDER;
        cache.header.root_node_id = tree.root.clone();
        cache.header.height = level;
        cache.header.first_node_offset = 4096; // always 4096
        cache.header.payload_start_offset = (4096 * (cache.header.total_nodes + 1)) as usize; // the start offset of the payload area, always 4096 * (total_nodes + 1)
        cache.finish()?;
        Ok(())
    }

    /// search the key in the B+ tree
    /// return the record pointer list
    /// if the key is not found, return None
    /// this function does not need chrom id since one bptree for one chromsome
    pub fn single_pos_search(&mut self, key: KeyType) -> Vec<MergedIsoformOffsetPtr> {
        // println!("search key: {}", &key);
        // let mut cache = cache::DataCache::from_disk("./test/block0.dat");
        let cache = self.cache.as_mut().expect("Can not get cache");
        let root = cache.get_root_node();

        let mut hit = match root.search(key) {
            None => {
                return vec![];
            }
            Some(hit) => {
                // println!("search node: {}", hit);
                hit
            }
        };

        loop {
            let node = cache.get_node2(hit).expect("Can not get node");
            // println!("search node: {}", node.header.id);
            if !node.is_leaf() {
                hit = match node.search(key) {
                    Some(hit) => hit,
                    None => return vec![],
                };
            } else {
                hit = match node.exact_search(key) {
                    Some(hit) => hit,
                    None => return vec![],
                };
                let addrs = node.get_addresses_by_children_value(&hit);
                return addrs;
            }
        }
    }

    pub fn range_search2(&mut self, pos: KeyType, flank: KeyType) -> Vec<MergedIsoformOffsetPtr> {
        let start = pos.saturating_sub(flank);
        let end = pos.saturating_add(flank);
        let cache = self.cache.as_mut().expect("Can not get cache");

        let mut node = {
            let mut n = cache.get_root_node();
            loop {
                if n.is_leaf() {
                    break n;
                }

                let slice = &n.header.keys[..n.header.num_keys as usize];
                let idx = match slice.binary_search(&start) {
                    Ok(idx) => idx,
                    Err(idx) => idx,
                };

                let child_id = n.header.childs[idx];
                // println!("searching node: {}, ", child_id);

                if child_id == 0 || child_id > cache.header.total_nodes {
                    // eprintln!("No such node: {}", child_id);
                    return vec![]; // no such node
                }

                n = cache.get_node2(child_id).expect("Can not get next node");
            }
        };

        // span the range in the leaf node

        let mut results = Vec::new();
        // results.reserve(estimated_size);

        loop {
            let keys = &node.header.keys[..node.header.num_keys as usize];
            for (nk, k) in keys.iter().enumerate() {
                if *k < start {
                    continue;
                }

                if *k > end {
                    return results;
                }

                let combined = node.header.childs[nk];
                results.extend(node.get_addresses_by_children_value(&combined));
            }

            let next_node_id = node.header.node_id + 1;
            if next_node_id > cache.header.leaf_nodes {
                break;
            }
            node = cache
                .get_node2(next_node_id)
                .expect("Can not get next leaf node");
        }
        results
    }
}

pub struct BPForest {
    pub index_dir: PathBuf,
    pub chrom_mapping: ChromMapping,
    // pub samples: AggrSampleMeta,
    pub trees_by_chrom: FxHashMap<u16, BPTree>,
}

impl BPForest {
    pub fn init(idx_path: &PathBuf) -> BPForest {
        let chrom_path = idx_path.join("chrom.map");
        let bytes = std::fs::read(chrom_path).unwrap();
        let chromamp = ChromMapping::decode(&bytes);

        BPForest {
            index_dir: idx_path.clone(),
            chrom_mapping: chromamp,
            // samples: AggrSampleMeta::load(&idx_path.join("sample.meta")),
            trees_by_chrom: FxHashMap::default(),
        }
    }

    pub fn build_forest(idx_path: &PathBuf) -> Result<()> {
        let chrom_path = idx_path.join("chrom.map");
        let bytes = std::fs::read(chrom_path).unwrap();
        let chromamp = ChromMapping::decode(&bytes);
        let intrim_path = idx_path.join("interim.idx");
        let mut count = 0;
        for chrom_id in chromamp.get_chrom_idxs() {
            let idx = Tmpindex::load(&intrim_path);
            // let blocks = idx.get_blocks(chrom_id);
            let blocks = idx.groups(chrom_id).unwrap();
            // if blocks.len() == 0 {
            //     continue;
            // }
            count += 1;
            BPTree::build_tree(blocks, &mut idx_path.clone(), chrom_id)?;
        }
        println!("Build {} trees", count);
        Ok(())
    }

    pub fn search_one_pos(
        &mut self,
        chrom_name: &str,
        pos: u64,
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        let chrom_id = match self.chrom_mapping.get_chrom_idx(chrom_name) {
            Some(id) => id,
            None => {
                return vec![];
            }
        };
        let tree: &mut BPTree = self
            .trees_by_chrom
            .entry(chrom_id)
            .or_insert_with(|| BPTree::from_disk(&self.index_dir, chrom_id, lru_size));

        // dbg!("aaa");
        tree.single_pos_search(pos)
    }

    pub fn search_one_range(
        &mut self,
        chrom_name: &str,
        pos: u64,
        flank: u64,
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        let chrom_id = match self.chrom_mapping.get_chrom_idx(chrom_name) {
            Some(id) => id,
            None => {
                return vec![];
            }
        };

        let tree: &mut BPTree = self
            .trees_by_chrom
            .entry(chrom_id)
            .or_insert_with(|| BPTree::from_disk(&self.index_dir, chrom_id, lru_size));

        tree.range_search2(pos, flank)
    }

    fn search_multi_exact(
        &mut self,
        positions: &Vec<(String, u64)>,
        min_match: usize,
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        let res_vec: Vec<Vec<MergedIsoformOffsetPtr>> = positions
            .iter()
            .map(|(chrom_name, pos)| {
                self.search_one_pos(&chrom_name.to_ascii_uppercase(), pos.clone(), lru_size)
            })
            .collect();

        if min_match == 0 || min_match >= positions.len() {
            return find_common(res_vec);
        } else {
            return find_partial_common(&res_vec, min_match);
        }
    }

    fn search_multi_range(
        &mut self,
        positions: &Vec<(String, u64)>,
        flank: u64,
        min_match: usize, // must larger than 0 and less than positions.len()
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        let res_vec: Vec<Vec<MergedIsoformOffsetPtr>> = positions
            .iter()
            .map(|(chrom_name, pos)| {
                self.search_one_range(
                    &chrom_name.to_ascii_uppercase(),
                    pos.clone(),
                    flank,
                    lru_size,
                )
            })
            .collect();

        // dbg!(&res_vec);

        if min_match == 0 || min_match >= positions.len() {
            find_common(res_vec)
        } else {
            find_partial_common(&res_vec, min_match)
        }
    }

    pub fn search_all_match(
        &mut self,
        positions: &Vec<(String, u64)>,
        flank: u64,
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        if flank == 0 {
            self.search_multi_exact(positions, 0, lru_size)
        } else {
            self.search_multi_range(positions, flank, 0, lru_size)
        }
    }

    pub fn search_partial_match(
        &mut self,
        positions: &Vec<(String, u64)>,
        flank: u64,
        min_match: usize, // must larger than 0 and less than positions.len()
        lru_size: usize,
    ) -> Vec<MergedIsoformOffsetPtr> {
        if flank == 0 {
            self.search_multi_exact(positions, min_match, lru_size)
        } else {
            self.search_multi_range(positions, flank, min_match, lru_size)
        }
    }
}

/// find the common elements in the vecs, if one element appears in at least min_match vecs, it is considered as common
/// min_match is set to 2 incase the mono exon isoforms
pub fn find_partial_common(
    vecs: &[Vec<MergedIsoformOffsetPtr>],
    min_match: usize,
) -> Vec<MergedIsoformOffsetPtr> {
    let mut count_map = FxHashMap::default();

    for vec in vecs {
        let unique_elements: HashSet<_> = vec.iter().cloned().collect();
        for elem in unique_elements {
            *count_map.entry(elem).or_insert(0) += 1;
        }
    }

    count_map
        .into_iter()
        .filter(|&(_, count)| count >= min_match)
        .map(|(item, _)| item)
        .collect()
}

pub fn find_common(mut vecs: Vec<Vec<MergedIsoformOffsetPtr>>) -> Vec<MergedIsoformOffsetPtr> {
    if vecs.is_empty() {
        return vec![];
    }

    // 先对每个 vec 排序
    for v in &mut vecs {
        v.sort();
        v.dedup(); // 避免重复
    }

    // 从第一个 vec 开始交集
    let mut result = vecs[0].clone();
    for v in vecs.iter().skip(1) {
        result = intersect_sorted(&result, v);
        if result.is_empty() {
            break;
        }
    }
    result
}

fn intersect_sorted<T: Ord + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let mut i = 0;
    let mut j = 0;
    let mut result = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                result.push(a[i].clone());
                i += 1;
                j += 1;
            }
        }
    }
    result
}
