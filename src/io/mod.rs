use anyhow::Result;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

pub struct InputReader {
    inner: InputReaderInner,
}

enum InputReaderInner {
    Stdin(BufReader<std::io::Stdin>),
    File(BufReader<File>),
    Gzip(BufReader<GzDecoder<File>>),
}

impl InputReader {
    pub fn open(path: &str) -> Result<Self> {
        let inner = if path == "-" || path.is_empty() {
            InputReaderInner::Stdin(BufReader::new(std::io::stdin()))
        } else if path.ends_with(".gz") {
            let file = File::open(path)?;
            InputReaderInner::Gzip(BufReader::new(GzDecoder::new(file)))
        } else {
            let file = File::open(path)?;
            InputReaderInner::File(BufReader::new(file))
        };
        Ok(Self { inner })
    }

    pub fn read_line(&mut self, buf: &mut String) -> Result<Option<usize>> {
        buf.clear();
        let n = match &mut self.inner {
            InputReaderInner::Stdin(r) => r.read_line(buf)?,
            InputReaderInner::File(r) => r.read_line(buf)?,
            InputReaderInner::Gzip(r) => r.read_line(buf)?,
        };
        Ok(if n == 0 { None } else { Some(n) })
    }

    pub fn read_header(&mut self) -> Result<(Vec<String>, Option<String>)> {
        let mut header_lines = Vec::new();
        let mut buf = String::new();
        let mut first_data_line: Option<String> = None;

        loop {
            buf.clear();
            let n = match &mut self.inner {
                InputReaderInner::Stdin(r) => r.read_line(&mut buf)?,
                InputReaderInner::File(r) => r.read_line(&mut buf)?,
                InputReaderInner::Gzip(r) => r.read_line(&mut buf)?,
            };

            if n == 0 {
                break;
            }

            let line = buf.trim_end().to_string();
            if line.starts_with('#') {
                header_lines.push(line);
            } else {
                first_data_line = Some(line);
                break;
            }
        }

        Ok((header_lines, first_data_line))
    }
}

pub struct OutputWriter {
    inner: OutputWriterInner,
}

enum OutputWriterInner {
    Stdout(BufWriter<std::io::Stdout>),
    File(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
}

impl OutputWriter {
    pub fn open(path: &str) -> Result<Self> {
        let inner = if path == "-" || path.is_empty() {
            OutputWriterInner::Stdout(BufWriter::new(std::io::stdout()))
        } else if path.ends_with(".gz") {
            let file = File::create(path)?;
            OutputWriterInner::Gzip(BufWriter::new(GzEncoder::new(file, Compression::default())))
        } else {
            let file = File::create(path)?;
            OutputWriterInner::File(BufWriter::new(file))
        };
        Ok(Self { inner })
    }

    pub fn write_line(&mut self, line: &str) -> Result<()> {
        let mut line_with_newline = String::with_capacity(line.len() + 1);
        line_with_newline.push_str(line);
        line_with_newline.push('\n');
        self.write_all(line_with_newline.as_bytes())?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        match &mut self.inner {
            OutputWriterInner::Stdout(w) => w.flush()?,
            OutputWriterInner::File(w) => w.flush()?,
            OutputWriterInner::Gzip(w) => w.flush()?,
        }
        Ok(())
    }
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.inner {
            OutputWriterInner::Stdout(w) => w.write(buf),
            OutputWriterInner::File(w) => w.write(buf),
            OutputWriterInner::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match &mut self.inner {
            OutputWriterInner::Stdout(w) => w.flush(),
            OutputWriterInner::File(w) => w.flush(),
            OutputWriterInner::Gzip(w) => w.flush(),
        }
    }
}
