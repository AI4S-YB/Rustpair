use crate::format::PairsamHeader;
use anyhow::{Context, Result};
use std::collections::HashMap;

pub fn get_chrom_order(
    chroms_path: &str,
    sam_chroms: Option<&[String]>,
) -> Result<HashMap<String, usize>> {
    let mut chrom_enum = HashMap::new();
    chrom_enum.insert(crate::format::UNMAPPED_CHROM.to_string(), 0usize);

    let content = std::fs::read_to_string(chroms_path)
        .with_context(|| format!("Failed to read chromosome file: {}", chroms_path))?;

    let mut i = 1usize;
    for line in content.lines() {
        let chrom = line.split_whitespace().next().unwrap_or("");
        if !chrom.is_empty() {
            if let Some(sam_chroms) = sam_chroms {
                if sam_chroms.iter().any(|c| c == chrom) {
                    chrom_enum.insert(chrom.to_string(), i);
                    i += 1;
                }
            } else {
                chrom_enum.insert(chrom.to_string(), i);
                i += 1;
            }
        }
    }

    if let Some(sam_chroms) = sam_chroms {
        let mut remaining: Vec<&String> = sam_chroms
            .iter()
            .filter(|c| !chrom_enum.contains_key(*c))
            .collect();
        remaining.sort();

        for chrom in remaining {
            chrom_enum.insert(chrom.clone(), i);
            i += 1;
        }
    }

    Ok(chrom_enum)
}

pub fn make_standard_pairsheader(
    assembly: &str,
    chromsizes: &[(String, u64)],
    columns: &[String],
    shape: &str,
) -> PairsamHeader {
    let mut header = PairsamHeader::new();
    header.genome_assembly = assembly.to_string();
    header.chromsizes = chromsizes.to_vec();
    header.columns = columns.to_vec();
    header.shape = shape.to_string();
    header
}

pub fn mark_header_as_sorted(mut header: PairsamHeader) -> PairsamHeader {
    if header.sorted.is_none() {
        header.sorted = Some("chr1-chr2-pos1-pos2".to_string());
    }
    header
}

pub fn extract_column_names(header: &PairsamHeader) -> Vec<String> {
    header.columns.clone()
}
