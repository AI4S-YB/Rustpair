use crate::format::{PairType, UNMAPPED_CHROM};
use anyhow::Result;
use log::debug;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct DedupConfig {
    pub max_mismatch: u32,
    pub method: DedupMethod,
    pub mark_dups: bool,
    pub keep_parent_id: bool,
    pub backend: DedupBackend,
    pub chunksize: usize,
    pub carryover: usize,
    pub n_proc: usize,
    pub c1_col: String,
    pub c2_col: String,
    pub p1_col: String,
    pub p2_col: String,
    pub s1_col: String,
    pub s2_col: String,
}

#[derive(Debug, Clone, Copy)]
pub enum DedupMethod {
    Max,
    Sum,
}

impl Default for DedupMethod {
    fn default() -> Self {
        Self::Max
    }
}

#[derive(Debug, Clone, Copy)]
pub enum DedupBackend {
    Scipy,
    Cython,
}

impl Default for DedupBackend {
    fn default() -> Self {
        Self::Scipy
    }
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            max_mismatch: 3,
            method: DedupMethod::Max,
            mark_dups: true,
            keep_parent_id: false,
            backend: DedupBackend::Scipy,
            chunksize: 10_000,
            carryover: 100,
            n_proc: 1,
            c1_col: "chrom1".to_string(),
            c2_col: "chrom2".to_string(),
            p1_col: "pos1".to_string(),
            p2_col: "pos2".to_string(),
            s1_col: "strand1".to_string(),
            s2_col: "strand2".to_string(),
        }
    }
}

pub struct OnlineDuplicateDetector {
    c1: Vec<i32>,
    c2: Vec<i32>,
    p1: Vec<i32>,
    p2: Vec<i32>,
    s1: Vec<i32>,
    s2: Vec<i32>,
    rm: Vec<bool>,
    parent_idxs: Vec<i32>,
    method_max: bool,
    max_mismatch: i32,
    low: usize,
    high: usize,
    n: usize,
    _return_data: bool,
    keep_parent_id: bool,
}

impl OnlineDuplicateDetector {
    pub fn new(method: &DedupMethod, max_mismatch: u32, keep_parent_id: bool) -> Self {
        let method_max = matches!(method, DedupMethod::Max);
        Self {
            c1: Vec::new(),
            c2: Vec::new(),
            p1: Vec::new(),
            p2: Vec::new(),
            s1: Vec::new(),
            s2: Vec::new(),
            rm: Vec::new(),
            parent_idxs: Vec::new(),
            method_max,
            max_mismatch: max_mismatch as i32,
            low: 0,
            high: 1,
            n: 0,
            _return_data: false,
            keep_parent_id,
        }
    }

    pub fn push(
        &mut self,
        c1: &[i32],
        c2: &[i32],
        p1: &[i32],
        p2: &[i32],
        s1: &[i32],
        s2: &[i32],
    ) -> (Vec<bool>, Option<Vec<i32>>) {
        self.c1.extend_from_slice(c1);
        self.c2.extend_from_slice(c2);
        self.p1.extend_from_slice(p1);
        self.p2.extend_from_slice(p2);
        self.s1.extend_from_slice(s1);
        self.s2.extend_from_slice(s2);
        self.rm.resize(self.n + c1.len(), false);

        if self.keep_parent_id {
            self.parent_idxs.extend(c1.iter().map(|_| 0i32));
        }

        self.n += c1.len();
        self.run(false)
    }

    pub fn finish(&mut self) -> (Vec<bool>, Option<Vec<i32>>) {
        self.run(true)
    }

    fn run(&mut self, finishing: bool) -> (Vec<bool>, Option<Vec<i32>>) {
        loop {
            if self.low == self.n {
                break;
            }

            if self.high == self.n {
                if finishing {
                    self.low += 1;
                    self.high = self.low + 1;
                    continue;
                } else {
                    break;
                }
            }

            if self.rm[self.low] {
                self.low += 1;
                self.high = self.low + 1;
                continue;
            }

            if self.rm[self.high] {
                self.high += 1;
                continue;
            }

            if self.c1[self.high] != self.c1[self.low]
                || (self.p1[self.high] - self.p1[self.low]) > self.max_mismatch
                || (self.p1[self.high] - self.p1[self.low]) < 0
            {
                self.low += 1;
                self.high = self.low + 1;
                continue;
            }

            let extra_condition = if self.method_max {
                (self.p1[self.low] - self.p1[self.high])
                    .abs()
                    .max((self.p2[self.low] - self.p2[self.high]).abs())
                    <= self.max_mismatch
            } else {
                ((self.p1[self.low] - self.p1[self.high]).abs()
                    + (self.p2[self.low] - self.p2[self.high]).abs())
                    <= self.max_mismatch
            };

            if self.c2[self.low] == self.c2[self.high]
                && self.s1[self.low] == self.s1[self.high]
                && self.s2[self.low] == self.s2[self.high]
                && extra_condition
            {
                self.rm[self.high] = true;
                if self.keep_parent_id {
                    self.parent_idxs[self.high] = self.low as i32;
                }
                self.high += 1;
                continue;
            }

            self.high += 1;
        }

        self.shrink()
    }

    fn shrink(&mut self) -> (Vec<bool>, Option<Vec<i32>>) {
        let result_rm = self.rm[..self.low].to_vec();

        let result_parents = if self.keep_parent_id {
            Some(self.parent_idxs[..self.low].to_vec())
        } else {
            None
        };

        let drain_n = self.low;

        self.c1.drain(..drain_n);
        self.c2.drain(..drain_n);
        self.p1.drain(..drain_n);
        self.p2.drain(..drain_n);
        self.s1.drain(..drain_n);
        self.s2.drain(..drain_n);
        self.rm.drain(..drain_n);

        if self.keep_parent_id {
            self.parent_idxs.drain(..drain_n);
        }

        self.n -= drain_n;
        self.high = self.high.saturating_sub(drain_n);
        self.low = 0;

        (result_rm, result_parents)
    }

    pub fn len(&self) -> usize {
        self.n
    }
}

pub fn streaming_dedup(
    input_path: &str,
    output_path: &str,
    output_dups_path: Option<&str>,
    config: &DedupConfig,
) -> Result<()> {
    debug!("Starting dedup with max_mismatch={}", config.max_mismatch);

    let mut reader = crate::io::InputReader::open(input_path)?;
    let (header_lines, first_data_line) = reader.read_header()?;

    let mut writer = crate::io::OutputWriter::open(output_path)?;
    let mut dups_writer = if let Some(dups_path) = output_dups_path {
        if dups_path == output_path || dups_path == "-" {
            None
        } else {
            Some(crate::io::OutputWriter::open(dups_path)?)
        }
    } else {
        None
    };

    for line in &header_lines {
        writer.write_line(line)?;
        if let Some(ref mut dw) = dups_writer {
            dw.write_line(line)?;
        }
    }

    let header = crate::format::PairsamHeader::parse_header_lines(&header_lines).0;
    let column_names = crate::header::extract_column_names(&header);

    let c1_idx = resolve_column_idx(&config.c1_col, &column_names);
    let c2_idx = resolve_column_idx(&config.c2_col, &column_names);
    let p1_idx = resolve_column_idx(&config.p1_col, &column_names);
    let p2_idx = resolve_column_idx(&config.p2_col, &column_names);
    let s1_idx = resolve_column_idx(&config.s1_col, &column_names);
    let s2_idx = resolve_column_idx(&config.s2_col, &column_names);
    let pt_idx = 7;

    let max_idx = *[c1_idx, c2_idx, p1_idx, p2_idx, s1_idx, s2_idx, pt_idx]
        .iter()
        .max()
        .unwrap();

    let mut dd =
        OnlineDuplicateDetector::new(&config.method, config.max_mismatch, config.keep_parent_id);

    let mut chrom_dict: HashMap<String, i32> = HashMap::new();
    let mut strand_dict: HashMap<String, i32> = HashMap::new();

    let mut line_buffer: Vec<String> = Vec::new();
    let mut cols_buffer: Vec<Vec<String>> = Vec::new();

    let mut c1_buf: Vec<i32> = Vec::new();
    let mut c2_buf: Vec<i32> = Vec::new();
    let mut p1_buf: Vec<i32> = Vec::new();
    let mut p2_buf: Vec<i32> = Vec::new();
    let mut s1_buf: Vec<i32> = Vec::new();
    let mut s2_buf: Vec<i32> = Vec::new();

    const MAX_LEN: usize = 10_000;

    let mut prefetched = first_data_line;
    let mut buf = String::new();

    loop {
        let stripline = if let Some(line) = prefetched.take() {
            line
        } else {
            match reader.read_line(&mut buf)? {
                None => break,
                Some(_) => buf.trim().to_string(),
            }
        };
        if stripline.is_empty() {
            continue;
        }

        let cols: Vec<&str> = stripline.split('\t').collect();
        if cols.len() <= max_idx {
            anyhow::bail!(
                "Error parsing line: expected {} columns, got {}",
                max_idx + 1,
                cols.len()
            );
        }

        let is_unmapped = cols[c1_idx] == UNMAPPED_CHROM || cols[c2_idx] == UNMAPPED_CHROM;

        if is_unmapped {
            writer.write_line(&stripline)?;
        } else {
            line_buffer.push(stripline.to_string());
            cols_buffer.push(cols.iter().map(|s| s.to_string()).collect());

            c1_buf.push(fetch_or_add(cols[c1_idx], &mut chrom_dict));
            c2_buf.push(fetch_or_add(cols[c2_idx], &mut chrom_dict));
            p1_buf.push(cols[p1_idx].parse::<i32>().unwrap_or(0));
            p2_buf.push(cols[p2_idx].parse::<i32>().unwrap_or(0));
            s1_buf.push(fetch_or_add(cols[s1_idx], &mut strand_dict));
            s2_buf.push(fetch_or_add(cols[s2_idx], &mut strand_dict));
        }

        if line_buffer.len() >= MAX_LEN {
            flush_dedup_buffer(
                &mut dd,
                &mut line_buffer,
                &mut cols_buffer,
                &mut c1_buf,
                &mut c2_buf,
                &mut p1_buf,
                &mut p2_buf,
                &mut s1_buf,
                &mut s2_buf,
                &mut writer,
                dups_writer.as_mut(),
                config,
                pt_idx,
                false,
            )?;
        }
    }

    if !line_buffer.is_empty() {
        flush_dedup_buffer(
            &mut dd,
            &mut line_buffer,
            &mut cols_buffer,
            &mut c1_buf,
            &mut c2_buf,
            &mut p1_buf,
            &mut p2_buf,
            &mut s1_buf,
            &mut s2_buf,
            &mut writer,
            dups_writer.as_mut(),
            config,
            pt_idx,
            true,
        )?;
    }

    writer.flush()?;
    if let Some(ref mut dw) = dups_writer {
        dw.flush()?;
    }

    debug!("Dedup complete");
    Ok(())
}

fn flush_dedup_buffer(
    dd: &mut OnlineDuplicateDetector,
    line_buffer: &mut Vec<String>,
    cols_buffer: &mut Vec<Vec<String>>,
    c1_buf: &mut Vec<i32>,
    c2_buf: &mut Vec<i32>,
    p1_buf: &mut Vec<i32>,
    p2_buf: &mut Vec<i32>,
    s1_buf: &mut Vec<i32>,
    s2_buf: &mut Vec<i32>,
    writer: &mut crate::io::OutputWriter,
    mut dups_writer: Option<&mut crate::io::OutputWriter>,
    config: &DedupConfig,
    pt_idx: usize,
    finishing: bool,
) -> Result<()> {
    let (res, parents_opt) = dd.push(c1_buf, c2_buf, p1_buf, p2_buf, s1_buf, s2_buf);

    let (res_final, parents): (Vec<bool>, Vec<i32>) = if finishing {
        let (r, p_finish) = dd.finish();
        let mut combined = res;
        combined.extend_from_slice(&r);
        match (parents_opt, p_finish) {
            (Some(mut po), Some(pf)) => {
                po.extend_from_slice(&pf);
                (combined, po)
            }
            (Some(po), None) => (combined, po),
            (None, Some(pf)) => (combined, pf),
            (None, None) => (combined, vec![]),
        }
    } else {
        (res, parents_opt.unwrap_or_default())
    };

    let n: usize = res_final.len();
    for i in 0..n {
        if !res_final[i] {
            writer.write_line(&line_buffer[i])?;
        } else {
            if let Some(ref mut dw) = dups_writer {
                if config.mark_dups && i < cols_buffer.len() && pt_idx < cols_buffer[i].len() {
                    cols_buffer[i][pt_idx] = PairType::DD.to_string();
                }
                let output = cols_buffer[i].join("\t");
                if config.keep_parent_id {
                    let parent_readid = if (parents[i] as usize) < line_buffer.len() {
                        line_buffer[parents[i] as usize]
                            .split('\t')
                            .next()
                            .unwrap_or("")
                            .to_string()
                    } else {
                        String::new()
                    };
                    dw.write_line(&format!("{}\t{}", output, parent_readid))?;
                } else {
                    dw.write_line(&output)?;
                }
            }
        }
    }

    line_buffer.drain(..n);
    cols_buffer.drain(..n);
    c1_buf.clear();
    c2_buf.clear();
    p1_buf.clear();
    p2_buf.clear();
    s1_buf.clear();
    s2_buf.clear();

    Ok(())
}

fn fetch_or_add(key: &str, dict: &mut HashMap<String, i32>) -> i32 {
    let key = key.trim();
    if let Some(&idx) = dict.get(key) {
        idx
    } else {
        let idx = dict.len() as i32;
        dict.insert(key.to_string(), idx);
        idx
    }
}

fn resolve_column_idx(col_name: &str, columns: &[String]) -> usize {
    if let Ok(idx) = col_name.parse::<usize>() {
        idx
    } else if let Some(pos) = columns.iter().position(|c| c == col_name) {
        pos
    } else {
        panic!("Column '{}' not found in header", col_name)
    }
}
