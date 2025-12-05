use std::{
    fs::File,
    io::{BufWriter, Write},
    num::NonZeroUsize,
    os::unix::fs::FileExt,
    path::Path,
    sync::Arc,
};

use anyhow::Result;
use lru::LruCache;

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

// pub fn read_record_from_mmap(
//     mmap: &[u8],
//     offset: &MergedIsoformOffsetPtr,
//     buf: &mut Vec<u8>,
// ) -> MergedIsoform {
//     let start = offset.offset as usize;
//     let end = start + offset.length as usize;
//     buf.clear();
//     buf.extend_from_slice(&mmap[start..end]);
//     match MergedIsoform::gz_decode(buf) {
//         Ok(record) => record,
//         Err(_) => {
//             eprintln!("Failed to decode gzipped record");
//             eprintln!("Offset: {:?}", offset);
//             std::process::exit(1);
//         }
//     }
// }

pub struct ArchiveCache {
    file: File,
    chunk_size: u64,
    lru_order: Vec<u64>,
    cache: LruCache<u64, Arc<Vec<u8>>>,
    buf: Vec<u8>,
}

impl ArchiveCache {
    pub fn new<P: AsRef<Path>>(file: P, chunk_size: u64, max_chunks: usize) -> ArchiveCache {
        let f = File::open(&file).expect(&format!(
            "Failed to open file {:?}",
            file.as_ref().display()
        ));
        ArchiveCache {
            file: f,
            chunk_size,
            lru_order: Vec::new(),
            cache: LruCache::new(NonZeroUsize::new(max_chunks).unwrap()),
            buf: Vec::new(),
        }
    }

    pub fn align_offset(&self, offset: u64) -> u64 {
        (offset / self.chunk_size) * self.chunk_size
    }

    pub fn read_chunk(&mut self, offset: u64) -> Result<Arc<Vec<u8>>> {
        let aligned = self.align_offset(offset);

        // if let Some(buf) = self.cache.get(&aligned) {
        //     return Ok(Arc::clone(buf));
        // }

        if let Some(buf) = self.cache.get(&aligned) {
            return Ok(Arc::clone(buf));
        }
        // dbg!(aligned);

        let mut buf = vec![0u8; self.chunk_size as usize];
        let bytes_read = self.file.read_at(&mut buf, aligned)?;
        buf.truncate(bytes_read);

        let arc_buf = Arc::new(buf);

        self.cache.put(aligned, Arc::clone(&arc_buf));

        Ok(arc_buf)
    }

    pub fn read_bytes(&mut self, offset: &MergedIsoformOffsetPtr) -> MergedIsoform {
        // info!("Reading record at offset: {:?}", offset);
        if self.buf.len() > 100 * 1024 * 1024 {
            // if buffer larger than 100MB, reset it
            self.buf.clear();
            self.buf.shrink_to_fit();
        } else {
            self.buf.clear();
        }

        let chunk = self
            .read_chunk(offset.offset)
            .expect(&format!("can not read record from {:?}", offset));

        let local_start = (offset.offset % self.chunk_size) as usize;
        let local_end = local_start + (offset.length as usize);

        // if end larger than chunk.len(), then we have problem, we need to load next chunk
        if local_end > chunk.len() {
            // load the partial data from current chunk
            self.buf.extend_from_slice(&chunk[local_start..]);

            let mut new_offset = self.align_offset(offset.offset);
            // iteratively load next chunks until we have enough data
            loop {
                if self.buf.len() == offset.length as usize {
                    break;
                }
                new_offset += self.chunk_size;

                let remaining = (offset.length as usize) - self.buf.len();

                let next_chunk = self.read_chunk(new_offset).expect(&format!(
                    "can not read record from next chunk at {:?}",
                    offset
                ));

                self.buf
                    .extend_from_slice(&next_chunk[0..remaining.min(next_chunk.len())]);
            }
        } else {
            self.buf.extend_from_slice(&chunk[local_start..local_end]);
        }

        match MergedIsoform::gz_decode(&self.buf) {
            Ok(record) => record,
            Err(_) => {
                eprintln!("Failed to decode gzipped record");
                eprintln!("Offset: {:?}", offset);
                std::process::exit(1);
            }
        }
    }

    pub fn clear_cache(&mut self) {
        self.cache.clear();
        self.lru_order.clear();
        self.buf.clear();

        self.lru_order.shrink_to_fit();
        self.buf.shrink_to_fit();
    }

    pub fn clear_buffer(&mut self) {
        self.buf.clear();
        if self.buf.capacity() > 10 * 1024 * 1024 {
            self.buf.shrink_to_fit();
        }
    }
}
