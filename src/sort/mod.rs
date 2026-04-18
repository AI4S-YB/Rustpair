use crate::format::PairsamHeader;
use crate::header;
use anyhow::Result;
use log::debug;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::io::Write;

#[derive(Debug, Clone)]
pub struct SortConfig {
    pub nproc: usize,
    pub tmpdir: Option<String>,
    pub memory: usize,
    pub c1_col: String,
    pub c2_col: String,
    pub p1_col: String,
    pub p2_col: String,
    pub pt_col: String,
    pub extra_cols: Vec<String>,
}

impl Default for SortConfig {
    fn default() -> Self {
        Self {
            nproc: 8,
            tmpdir: None,
            memory: 2 * 1024 * 1024 * 1024,
            c1_col: "chrom1".to_string(),
            c2_col: "chrom2".to_string(),
            p1_col: "pos1".to_string(),
            p2_col: "pos2".to_string(),
            pt_col: "pair_type".to_string(),
            extra_cols: Vec::new(),
        }
    }
}

pub fn sort_pairs(input_path: &str, output_path: &str, config: &SortConfig) -> Result<()> {
    let mut reader = crate::io::InputReader::open(input_path)?;
    let (header_lines, first_data_line) = reader.read_header()?;

    let (mut header, _) = PairsamHeader::parse_header_lines(&header_lines);
    header = header::mark_header_as_sorted(header);

    let column_names = header::extract_column_names(&header);

    let c1_idx = resolve_column_index(&config.c1_col, &column_names);
    let c2_idx = resolve_column_index(&config.c2_col, &column_names);
    let p1_idx = resolve_column_index(&config.p1_col, &column_names);
    let p2_idx = resolve_column_index(&config.p2_col, &column_names);
    let pt_idx = resolve_column_index(&config.pt_col, &column_names);

    debug!(
        "Sorting by columns: c1={}, c2={}, p1={}, p2={}, pt={}",
        c1_idx, c2_idx, p1_idx, p2_idx, pt_idx
    );

    let mut all_lines: Vec<String> = Vec::new();
    let mut buf = String::new();

    if let Some(line) = first_data_line {
        if !line.trim().is_empty() {
            all_lines.push(line.trim().to_string());
        }
    }

    loop {
        match reader.read_line(&mut buf)? {
            None => break,
            Some(_) => {
                if !buf.trim().is_empty() && !buf.starts_with('#') {
                    all_lines.push(buf.trim().to_string());
                }
            }
        }
    }

    debug!("Read {} lines to sort", all_lines.len());

    if all_lines.is_empty() {
        let mut writer = crate::io::OutputWriter::open(output_path)?;
        for line in header.to_lines() {
            writer.write_line(&line)?;
        }
        writer.flush()?;
        return Ok(());
    }

    let chunk_size = (all_lines.len() / config.nproc).max(1000);

    let sorted_chunks: Vec<Vec<String>> = all_lines
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_chunk: Vec<&str> = chunk.iter().map(|s| s.as_str()).collect();
            local_chunk.sort_by(|a, b| {
                let fa: Vec<&str> = a.split('\t').collect();
                let fb: Vec<&str> = b.split('\t').collect();
                let mut ord = Ordering::Equal;
                if c1_idx < fa.len() && c1_idx < fb.len() {
                    ord = fa[c1_idx].cmp(&fb[c1_idx]);
                }
                if ord == Ordering::Equal && c2_idx < fa.len() && c2_idx < fb.len() {
                    ord = fa[c2_idx].cmp(&fb[c2_idx]);
                }
                if ord == Ordering::Equal {
                    let p1_a: u64 = fa.get(p1_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let p1_b: u64 = fb.get(p1_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
                    ord = p1_a.cmp(&p1_b);
                }
                if ord == Ordering::Equal {
                    let p2_a: u64 = fa.get(p2_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
                    let p2_b: u64 = fb.get(p2_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
                    ord = p2_a.cmp(&p2_b);
                }
                if ord == Ordering::Equal && pt_idx < fa.len() && pt_idx < fb.len() {
                    ord = fa[pt_idx].cmp(&fb[pt_idx]);
                }
                ord
            });
            local_chunk.into_iter().map(|s| s.to_string()).collect()
        })
        .collect();

    let merged = kway_merge(sorted_chunks, c1_idx, c2_idx, p1_idx, p2_idx, pt_idx);

    let mut writer = crate::io::OutputWriter::open(output_path)?;

    for line in header.to_lines() {
        writer.write_line(&line)?;
    }

    for line in merged {
        writeln!(writer, "{}", line)?;
    }

    writer.flush()?;
    debug!("Sorting complete");
    Ok(())
}

fn resolve_column_index(col_name: &str, columns: &[String]) -> usize {
    if let Ok(idx) = col_name.parse::<usize>() {
        idx
    } else if let Some(pos) = columns.iter().position(|c| c == col_name) {
        pos
    } else {
        panic!("Column '{}' not found in header", col_name)
    }
}

fn kway_merge(
    chunks: Vec<Vec<String>>,
    c1_idx: usize,
    c2_idx: usize,
    p1_idx: usize,
    p2_idx: usize,
    pt_idx: usize,
) -> Vec<String> {
    use std::collections::BinaryHeap;

    #[derive(Clone)]
    struct HeapEntry {
        line: String,
        chunk_idx: usize,
        line_idx: usize,
        c1_idx: usize,
        c2_idx: usize,
        p1_idx: usize,
        p2_idx: usize,
        pt_idx: usize,
    }

    impl PartialEq for HeapEntry {
        fn eq(&self, _other: &Self) -> bool {
            false
        }
    }

    impl Eq for HeapEntry {}

    impl PartialOrd for HeapEntry {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for HeapEntry {
        fn cmp(&self, other: &Self) -> Ordering {
            let fields_a: Vec<&str> = other.line.split('\t').collect();
            let fields_b: Vec<&str> = self.line.split('\t').collect();

            let mut ord = Ordering::Equal;

            if self.c1_idx < fields_b.len() && self.c1_idx < fields_a.len() {
                ord = fields_a[self.c1_idx].cmp(fields_b[self.c1_idx]);
            }
            if ord == Ordering::Equal
                && self.c2_idx < fields_b.len()
                && self.c2_idx < fields_a.len()
            {
                ord = fields_a[self.c2_idx].cmp(fields_b[self.c2_idx]);
            }
            if ord == Ordering::Equal {
                let pos1_a: u64 = fields_a
                    .get(self.p1_idx)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                let pos1_b: u64 = fields_b
                    .get(self.p1_idx)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                ord = pos1_a.cmp(&pos1_b);
            }
            if ord == Ordering::Equal {
                let pos2_a: u64 = fields_a
                    .get(self.p2_idx)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                let pos2_b: u64 = fields_b
                    .get(self.p2_idx)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                ord = pos2_a.cmp(&pos2_b);
            }
            if ord == Ordering::Equal
                && self.pt_idx < fields_a.len()
                && self.pt_idx < fields_b.len()
            {
                ord = fields_a[self.pt_idx].cmp(fields_b[self.pt_idx]);
            }

            ord
        }
    }

    let mut heap = BinaryHeap::new();
    let mut chunk_iters: Vec<std::vec::IntoIter<String>> =
        chunks.into_iter().map(|c| c.into_iter()).collect();

    for (i, iter) in chunk_iters.iter_mut().enumerate() {
        if let Some(line) = iter.next() {
            heap.push(HeapEntry {
                line,
                chunk_idx: i,
                line_idx: 0,
                c1_idx,
                c2_idx,
                p1_idx,
                p2_idx,
                pt_idx,
            });
        }
    }

    let mut result = Vec::new();
    while let Some(entry) = heap.pop() {
        result.push(entry.line);

        if let Some(next_line) = chunk_iters[entry.chunk_idx].next() {
            heap.push(HeapEntry {
                line: next_line,
                chunk_idx: entry.chunk_idx,
                line_idx: entry.line_idx + 1,
                c1_idx: entry.c1_idx,
                c2_idx: entry.c2_idx,
                p1_idx: entry.p1_idx,
                p2_idx: entry.p2_idx,
                pt_idx: entry.pt_idx,
            });
        }
    }

    result
}
