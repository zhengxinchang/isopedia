use std::{
    fs::File,
    io::{BufReader, BufWriter, Read, Seek, Write}, path::Path,
};

use crate::{isoform::MergedIsoform, tmpidx::MergedIsoformOffset};

/// extract the code about the data storage and place in this file.
///
///
///
///
///
///

pub struct IsoformArchive {
    writer: BufWriter<File>,
    buf: Vec<u8>,
}

impl IsoformArchive {
    pub fn create(path: &Path) -> IsoformArchive {
        let file = File::create(path).expect("Failed to open file");
        let writer = BufWriter::new(file);
        let buf = Vec::with_capacity(1024 * 1024); // 1MB buffer
        IsoformArchive { writer, buf }
    }

    pub fn dump_to_disk(&mut self, record: &MergedIsoform) -> u32 {
        self.buf.clear();
        let byte_len = record.gz_encode(&mut self.buf);
        self.writer
            .write_all(&self.buf)
            .expect("Failed to write record");
        byte_len
    }
}

pub fn read_record_from_aggr_file(
    reader: &mut BufReader<File>,
    record_ptr: &MergedIsoformOffset,
) -> MergedIsoform {
    let mut buf = Vec::new();
    reader
        .seek(std::io::SeekFrom::Start(record_ptr.offset))
        .unwrap();
    reader
        .take(record_ptr.length as u64)
        .read_to_end(&mut buf)
        .unwrap();

    MergedIsoform::gz_decode(&buf).unwrap()
}
