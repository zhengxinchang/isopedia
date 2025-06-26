use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use flate2::{write::GzEncoder, Compression};

pub struct MyGzWriter {
    file_name: String,
    inner: Option<GzEncoder<BufWriter<File>>>,
}

impl MyGzWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(&path)?;

        let fname = path.as_ref().to_string_lossy();
        let writer = GzEncoder::new(BufWriter::new(file), Compression::default());

        Ok(MyGzWriter {
            file_name: fname.to_string(),
            inner: Some(writer),
        })
    }

    pub fn write_all_bytes(&mut self, bytes: &[u8]) -> io::Result<()> {
        match &mut self.inner {
            Some(writer) => writer.write_all(bytes),
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been finished",
            )),
        }
    }

    pub fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            Some(writer) => writer.flush(),
            None => Ok(()), // 已完成，无需 flush
        }
    }

    pub fn finish(&mut self) -> io::Result<()> {
        if let Some(writer) = self.inner.take() {
            // take() 移出所有权
            // let mut buf_writer = 
            writer.finish()?;
            // buf_writer.flush()?;
        }
        Ok(())
    }
}

impl Drop for MyGzWriter {
    fn drop(&mut self) {
        self.finish().expect("Can not drop gz writer...");
    }
}
