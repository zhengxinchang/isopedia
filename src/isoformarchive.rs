use std::{
    fs::File,
    io::{BufReader, BufWriter, Read, Seek, Write}, path::Path,
};

use crate::{isoform::MergedIsoform, tmpidx::MergedIsoformOffset};

/// extract the code about the data storage and place in this file.

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

// pub fn read_record_from_archive(
//     reader: &mut BufReader<File>,
//     record_ptr: &MergedIsoformOffset,
// ) -> MergedIsoform {
//     let mut buf = Vec::new();
//     reader
//         .seek(std::io::SeekFrom::Start(record_ptr.offset))
//         .expect("Failed to seek in archive file");
//     reader
//         .take(record_ptr.length as u64)
//         .read_to_end(&mut buf)
//         .expect("Failed to read record from archive file");

//     MergedIsoform::gz_decode(&buf).expect("Failed to decode gzipped record")
// }


pub fn read_record_from_archive(
    reader: &mut BufReader<File>,
    record_ptr: &MergedIsoformOffset,
    buf: &mut Vec<u8>, // ✅ 外部提供可复用 buffer
) -> MergedIsoform {
    buf.clear(); // ✅ 避免重新分配

    // ✅ 移动光标到记录起始
    reader
        .seek(std::io::SeekFrom::Start(record_ptr.offset))
        .expect("Failed to seek in archive file");

    // ✅ 只读取指定字节数
    reader
        .take(record_ptr.length as u64)
        .read_to_end(buf)
        .expect("Failed to read record from archive file");

    // ✅ 解压 gzipped 数据
    MergedIsoform::gz_decode(buf).expect("Failed to decode gzipped record")
}
