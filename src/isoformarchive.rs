use std::{
    collections::HashMap,
    fs::File,
    hash::Hash,
    io::{BufWriter, Write},
    os::unix::fs::FileExt,
    path::Path,
    sync::Arc,
};

use anyhow::Result;

use crate::{isoform::MergedIsoform, tmpidx::MergedIsoformOffsetPtr};

pub struct IsoformArchiveWriter {
    writer: BufWriter<File>,
    buf: Vec<u8>,
}

impl IsoformArchiveWriter {
    pub fn create(path: &Path) -> IsoformArchiveWriter {
        let file = File::create(path).expect("Failed to open file");
        let writer = BufWriter::with_capacity(10 * 1024 * 1024, file);
        let buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
        IsoformArchiveWriter { writer, buf }
    }

    pub fn dump_to_disk(&mut self, record: &MergedIsoform) -> u32 {
        self.buf.clear();
        let byte_len = record.gz_encode(&mut self.buf);
        self.writer
            .write_all(&self.buf)
            .expect("Failed to write record");
        byte_len
    }

    pub fn close_file(&mut self) -> std::io::Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

pub fn read_record_from_mmap(
    mmap: &[u8],
    offset: &MergedIsoformOffsetPtr,
    buf: &mut Vec<u8>,
) -> MergedIsoform {
    let start = offset.offset as usize;
    let end = start + offset.length as usize;
    buf.clear();
    buf.extend_from_slice(&mmap[start..end]);
    match MergedIsoform::gz_decode(buf) {
        Ok(record) => record,
        Err(_) => {
            eprintln!("Failed to decode gzipped record");
            eprintln!("Offset: {:?}", offset);
            std::process::exit(1);
        }
    }
}

pub struct ArchiveCache {
    file: File,
    chunk_size: usize,
    lru_order: Vec<u64>,
    cache: HashMap<u64, Arc<Vec<u8>>>,
    max_chunks: usize,
    buf: Vec<u8>,
}

impl ArchiveCache {
    pub fn new<P: AsRef<Path>>(file: P, chunk_size: usize, max_chunks: usize) -> ArchiveCache {
        let f = File::open(&file).expect(&format!(
            "Failed to open file {:?}",
            file.as_ref().display()
        ));
        ArchiveCache {
            file: f,
            chunk_size,
            lru_order: Vec::new(),
            cache: HashMap::new(),
            max_chunks,
            buf: Vec::new(),
        }
    }

    pub fn read_chunk(&mut self, offset: u64) -> Result<Arc<Vec<u8>>> {
        let aligned = (offset / self.chunk_size as u64) * self.chunk_size as u64;

        if let Some(buf) = self.cache.get(&aligned) {
            return Ok(Arc::clone(buf));
        }

        // 新块
        let mut buf = vec![0u8; self.chunk_size];
        let bytes_read = self.file.read_at(&mut buf, aligned)?;
        buf.truncate(bytes_read);

        let arc_buf = Arc::new(buf);
        self.cache.insert(aligned, Arc::clone(&arc_buf));
        self.lru_order.push(aligned);

        // 控制 LRU 缓存数量
        if self.lru_order.len() > self.max_chunks {
            let old = self.lru_order.remove(0);
            self.cache.remove(&old);
        }

        Ok(arc_buf)
    }

    pub fn read_bytes(&mut self, offset: &MergedIsoformOffsetPtr) -> MergedIsoform {
        let chunk = self
            .read_chunk(offset.offset)
            .expect(&format!("can not read record from {:?}", offset));
        let start = (offset.offset % self.chunk_size as u64) as usize;
        let end = std::cmp::min(start + (offset.length as usize), chunk.len());
        self.buf.clear();
        self.buf.extend_from_slice(&chunk[start..end]);

        match MergedIsoform::gz_decode(&self.buf) {
            Ok(record) => record,
            Err(_) => {
                eprintln!("Failed to decode gzipped record");
                eprintln!("Offset: {:?}", offset);
                std::process::exit(1);
            }
        }
    }

    // pub fn clear(&mut self) {
    //     self.cache.clear();
    //     self.lru_order.clear();
    // }
}
