use crate::format::{Alignment, PairType, UNMAPPED_CHROM, UNMAPPED_POS, UNMAPPED_STRAND};
use crate::header;
use anyhow::{Context, Result};
use log::debug;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct ParseConfig {
    pub min_mapq: u8,
    pub max_molecule_size: u32,
    pub max_inter_align_gap: u32,
    pub walks_policy: WalksPolicy,
    pub drop_readid: bool,
    pub drop_seq: bool,
    pub drop_sam: bool,
    pub add_pair_index: bool,
    pub add_columns: Vec<String>,
    pub report_alignment_end: char,
    pub flip: bool,
}

#[derive(Debug, Clone, Copy)]
pub enum WalksPolicy {
    Mask,
    FiveAny,
    FiveUnique,
    ThreeAny,
    ThreeUnique,
    All,
}

impl Default for ParseConfig {
    fn default() -> Self {
        Self {
            min_mapq: 1,
            max_molecule_size: 750,
            max_inter_align_gap: 20,
            walks_policy: WalksPolicy::FiveUnique,
            drop_readid: false,
            drop_seq: false,
            drop_sam: false,
            add_pair_index: false,
            add_columns: Vec::new(),
            report_alignment_end: '5',
            flip: true,
        }
    }
}

const BAM_FPAIRED: u16 = 1;
const _BAM_FPROPER_PAIR: u16 = 2;
const FUNMAP: u16 = 4;
const FREVERSE: u16 = 16;
const FSUPPLEMENTARY: u16 = 2048;
const FFIRST_MATE: u16 = 64;

pub fn streaming_classify(
    sam_path: &str,
    chroms_path: &str,
    output_path: &str,
    config: &ParseConfig,
) -> Result<()> {
    let mut reader = bam::Reader::from_path(sam_path)
        .with_context(|| format!("Failed to open SAM/BAM file: {}", sam_path))?;

    let header_view = reader.header().clone();

    let mut sam_chromsizes: HashMap<String, i64> = HashMap::new();
    for tid in 0..header_view.target_count() {
        if let (Some(name), Some(len)) = (
            std::str::from_utf8(header_view.tid2name(tid)).ok(),
            header_view.target_len(tid),
        ) {
            sam_chromsizes.insert(name.to_string(), len as i64);
        }
    }

    let chrom_names: Vec<String> = sam_chromsizes.keys().cloned().collect();
    let chromosomes = header::get_chrom_order(chroms_path, Some(&chrom_names))?;

    let mut chrom_list: Vec<(String, usize)> = chromosomes
        .iter()
        .filter(|(k, _)| **k != UNMAPPED_CHROM)
        .map(|(k, v)| (k.clone(), *v))
        .collect();
    chrom_list.sort_by_key(|(_, v)| *v);
    let chrom_list: Vec<String> = chrom_list.into_iter().map(|(k, _)| k).collect();

    let chromsizes: Vec<(String, u64)> = chrom_list
        .iter()
        .filter_map(|c| sam_chromsizes.get(c).map(|&size| (c.clone(), size as u64)))
        .collect();

    let columns = if config.drop_sam {
        vec![
            "readID".into(),
            "chrom1".into(),
            "pos1".into(),
            "chrom2".into(),
            "pos2".into(),
            "strand1".into(),
            "strand2".into(),
            "pair_type".into(),
        ]
    } else {
        vec![
            "readID".into(),
            "chrom1".into(),
            "pos1".into(),
            "chrom2".into(),
            "pos2".into(),
            "strand1".into(),
            "strand2".into(),
            "pair_type".into(),
            "sam1".into(),
            "sam2".into(),
        ]
    };

    let pairsam_header = header::make_standard_pairsheader(
        "",
        &chromsizes,
        &columns,
        if config.flip {
            "upper triangle"
        } else {
            "whole matrix"
        },
    );

    let mut writer = crate::io::OutputWriter::open(output_path)?;

    for line in pairsam_header.to_lines() {
        writer.write_line(&line)?;
    }

    debug!("Starting to parse {}", sam_path);

    let mut prev_read_id: Option<String> = None;
    let mut current_sams1: Vec<bam::Record> = Vec::new();
    let mut current_sams2: Vec<bam::Record> = Vec::new();
    let mut record = bam::Record::new();

    loop {
        let read_result = reader.read(&mut record);
        match read_result {
            Some(Ok(())) => {}
            Some(Err(e)) => return Err(e).context("Error reading BAM record"),
            None => break,
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();

        if let Some(ref prev_id) = prev_read_id {
            if read_id != *prev_id {
                process_read(
                    &current_sams1,
                    &current_sams2,
                    prev_id,
                    config,
                    &chromosomes,
                    &header_view,
                    &mut writer,
                )?;
                current_sams1.clear();
                current_sams2.clear();
            }
        }

        prev_read_id = Some(read_id.clone());
        let flag = record.flags();

        if (flag & FFIRST_MATE != 0) || ((flag & BAM_FPAIRED == 0) && (flag & FSUPPLEMENTARY == 0))
        {
            current_sams1.push(record.clone());
        } else {
            current_sams2.push(record.clone());
        }
    }

    if let Some(read_id) = prev_read_id {
        process_read(
            &current_sams1,
            &current_sams2,
            &read_id,
            config,
            &chromosomes,
            &header_view,
            &mut writer,
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn process_read(
    sams1: &[bam::Record],
    sams2: &[bam::Record],
    read_id: &str,
    config: &ParseConfig,
    chrom_enum: &HashMap<String, usize>,
    header_view: &bam::HeaderView,
    writer: &mut crate::io::OutputWriter,
) -> Result<()> {
    let is_empty = match config.walks_policy {
        WalksPolicy::All => {
            (sams1.is_empty() && sams2.len() < 2) || (sams2.is_empty() && sams1.len() < 2)
        }
        _ => sams1.is_empty() || sams2.is_empty(),
    };

    if is_empty {
        let mut algn1 = Alignment::empty();
        let mut algn2 = Alignment::empty();
        algn1.pair_type_char = 'X';
        algn2.pair_type_char = 'X';
        write_pair(&algn1, &algn2, read_id, config, writer)?;
        return Ok(());
    }

    let mut algns1: Vec<Alignment> = sams1
        .iter()
        .map(|s| parse_pysam_entry(s, config.min_mapq, &header_view))
        .collect();

    let mut algns2: Vec<Alignment> = sams2
        .iter()
        .map(|s| parse_pysam_entry(s, config.min_mapq, &header_view))
        .collect();

    normalize_alignment_list(&mut algns1, 1, config.max_inter_align_gap);
    normalize_alignment_list(&mut algns2, 2, config.max_inter_align_gap);

    let is_chimeric_1 = algns1.len() > 1;
    let is_chimeric_2 = algns2.len() > 1;

    let pairs: Vec<(Alignment, Alignment, (u32, &'static str))> = if is_chimeric_1 || is_chimeric_2
    {
        match config.walks_policy {
            WalksPolicy::All => {
                let walk_pairs: Vec<(Alignment, Alignment, (u32, &'static str))> =
                    parse_complex_walk(
                        &algns1,
                        &algns2,
                        config.max_molecule_size,
                        "outer",
                        "pair",
                    )?
                    .collect();
                walk_pairs
            }
            _ => {
                let rescued = rescue_walk(&algns1, &algns2, config.max_molecule_size);
                let pair_index = match rescued {
                    Some(1) => (1, "R1"),
                    Some(2) => (1, "R2"),
                    Some(_) => (1, "R1-2"),
                    None => (1, "R1-2"),
                };

                let (mut hic_algn1, mut hic_algn2) = (algns1[0].clone(), algns2[0].clone());

                if rescued.is_none() {
                    match config.walks_policy {
                        WalksPolicy::Mask => {
                            hic_algn1.mask();
                            hic_algn2.mask();
                            hic_algn1.pair_type_char = 'W';
                            hic_algn2.pair_type_char = 'W';
                        }
                        WalksPolicy::FiveAny => {}
                        WalksPolicy::FiveUnique => {
                            for a in &algns1 {
                                if a.is_mapped && a.is_unique {
                                    hic_algn1 = a.clone();
                                    break;
                                }
                            }
                            for a in &algns2 {
                                if a.is_mapped && a.is_unique {
                                    hic_algn2 = a.clone();
                                    break;
                                }
                            }
                        }
                        WalksPolicy::ThreeAny => {
                            if !algns1.is_empty() {
                                hic_algn1 = algns1.last().unwrap().clone();
                            }
                            if !algns2.is_empty() {
                                hic_algn2 = algns2.last().unwrap().clone();
                            }
                        }
                        WalksPolicy::ThreeUnique => {
                            for a in algns1.iter().rev() {
                                if a.is_mapped && a.is_unique {
                                    hic_algn1 = a.clone();
                                    break;
                                }
                            }
                            for a in algns2.iter().rev() {
                                if a.is_mapped && a.is_unique {
                                    hic_algn2 = a.clone();
                                    break;
                                }
                            }
                        }
                        _ => unreachable!(),
                    }
                } else {
                    if rescued == Some(1) {
                        if algns2.len() > 1 {
                            algns2[1].pair_type_char = 'X';
                        }
                        hic_algn2.pair_type_char = 'R';
                    } else {
                        if algns1.len() > 1 {
                            algns1[1].pair_type_char = 'X';
                        }
                        hic_algn1.pair_type_char = 'R';
                    }
                }

                vec![(hic_algn1, hic_algn2, pair_index)]
            }
        }
    } else {
        vec![(
            algns1.get(0).cloned().unwrap_or_default(),
            algns2.get(0).cloned().unwrap_or_default(),
            (1, "R1-2"),
        )]
    };

    for (mut algn1, mut algn2, _pair_index) in pairs {
        if config.report_alignment_end == '5' {
            algn1.pos = algn1.pos5;
            algn2.pos = algn2.pos5;
        } else {
            algn1.pos = algn1.pos3;
            algn2.pos = algn2.pos3;
        }

        if config.flip {
            if !check_pair_order(&algn1, &algn2, chrom_enum) {
                std::mem::swap(&mut algn1, &mut algn2);
            }
        }

        write_pair(&algn1, &algn2, read_id, config, writer)?;
    }

    Ok(())
}

fn write_pair(
    algn1: &Alignment,
    algn2: &Alignment,
    read_id: &str,
    config: &ParseConfig,
    writer: &mut crate::io::OutputWriter,
) -> Result<()> {
    let type1 = algn1.type_str();
    let type2 = algn2.type_str();
    let pair_type = PairType::from_chars(type1, type2);

    let mut fields = Vec::new();

    if !config.drop_readid {
        fields.push(read_id.to_string());
    }
    fields.push(algn1.chrom.clone());
    fields.push(algn1.pos.to_string());
    fields.push(algn2.chrom.clone());
    fields.push(algn2.pos.to_string());
    fields.push(algn1.strand.to_string());
    fields.push(algn2.strand.to_string());
    fields.push(pair_type.to_string());

    writer.write_line(&fields.join("\t"))?;

    Ok(())
}

fn parse_pysam_entry(sam: &bam::Record, min_mapq: u8, hv: &bam::HeaderView) -> Alignment {
    let mut algn = Alignment::empty();

    let flag = sam.flags();
    algn.is_mapped = (flag & FUNMAP) == 0;
    algn.mapq = sam.mapq();
    algn.is_unique = algn.mapq >= min_mapq;

    let cigar = sam.cigar();
    let mut matched_bp = 0u32;
    let mut algn_ref_span = 0u32;
    let mut algn_read_span = 0u32;
    let mut read_len = 0u32;
    let mut clip5_ref = 0u32;
    let mut clip3_ref = 0u32;

    for c in cigar.iter() {
        match c.char() {
            'M' | '=' | 'X' => {
                matched_bp += c.len();
                algn_ref_span += c.len();
                algn_read_span += c.len();
                read_len += c.len();
            }
            'I' => {
                algn_read_span += c.len();
                read_len += c.len();
            }
            'D' | 'N' => {
                algn_ref_span += c.len();
            }
            'S' | 'H' => {
                read_len += c.len();
                if matched_bp == 0 {
                    clip5_ref = c.len();
                } else {
                    clip3_ref = c.len();
                }
            }
            'P' => {}
            _ => {}
        }
    }

    algn.matched_bp = matched_bp;
    algn.algn_ref_span = algn_ref_span;
    algn.algn_read_span = algn_read_span;
    algn.read_len = read_len;
    algn.clip5_ref = clip5_ref;
    algn.clip3_ref = clip3_ref;

    if algn.is_mapped {
        algn.is_linear = aux_tag_exists(sam, b"SA");

        if (flag & FREVERSE) == 0 {
            algn.strand = '+';
            algn.dist_to_5 = clip5_ref;
            algn.dist_to_3 = clip3_ref;
        } else {
            algn.strand = '-';
            algn.dist_to_5 = clip3_ref;
            algn.dist_to_3 = clip5_ref;
        }

        if algn.is_unique {
            let tid = sam.tid();
            if tid >= 0 {
                if let Some(name_bytes) = std::str::from_utf8(hv.tid2name(tid as u32)).ok() {
                    algn.chrom = name_bytes.to_string();
                }

                if algn.strand == '+' {
                    algn.pos5 = (sam.pos() + 1) as u32;
                    algn.pos3 = sam.pos() as u32 + algn.algn_ref_span;
                } else {
                    algn.pos5 = sam.pos() as u32 + algn.algn_ref_span;
                    algn.pos3 = (sam.pos() + 1) as u32;
                }
            }
        } else {
            algn.chrom = UNMAPPED_CHROM.to_string();
            algn.strand = UNMAPPED_STRAND.chars().next().unwrap_or('-');
            algn.pos5 = UNMAPPED_POS;
            algn.pos3 = UNMAPPED_POS;
        }
    }

    algn.pair_type_char = algn.type_str();
    algn
}

fn aux_tag_exists(record: &bam::Record, tag: &[u8]) -> bool {
    record.aux(tag).is_ok()
}

fn normalize_alignment_list(algns: &mut Vec<Alignment>, _side: u32, max_inter_align_gap: u32) {
    if algns.is_empty() {
        algns.push(Alignment::empty());
        return;
    }

    algns.sort_by_key(|a| a.dist_to_5);

    if max_inter_align_gap > 0 && algns.len() > 1 {
        convert_gaps_into_alignments(algns, max_inter_align_gap);
    }

    for i in 0..algns.len() {
        algns[i].algn_idx = i as i32;
        algns[i].same_side_algn_count = algns.len() as i32;
    }
}

fn convert_gaps_into_alignments(algns: &mut Vec<Alignment>, max_inter_align_gap: u32) {
    if algns.len() == 1 && !algns[0].is_mapped {
        return;
    }

    let mut last_5_pos: u32 = 0;
    let mut i = 0;
    while i < algns.len() {
        let algn = &algns[i];
        if algn.dist_to_5.saturating_sub(last_5_pos) > max_inter_align_gap && algn.is_mapped {
            let mut new_algn = Alignment::empty();
            new_algn.dist_to_5 = last_5_pos;
            new_algn.algn_read_span = algn.dist_to_5.saturating_sub(last_5_pos);
            new_algn.read_len = algn.read_len;
            new_algn.dist_to_3 = new_algn.read_len.saturating_sub(algn.dist_to_5);
            algns.insert(i, new_algn);
            i += 2;
        } else {
            last_5_pos = last_5_pos.max(algn.dist_to_5 + algn.algn_read_span);
            i += 1;
        }
    }
}

fn rescue_walk(algns1: &[Alignment], algns2: &[Alignment], max_molecule_size: u32) -> Option<u32> {
    let n1 = algns1.len();
    let n2 = algns2.len();

    if n1 <= 1 && n2 <= 1 {
        return None;
    }

    if !((n1 == 1 && n2 == 2) || (n1 == 2 && n2 == 1)) {
        return None;
    }

    let first_read_is_chimeric = n1 > 1;
    let chim5 = if first_read_is_chimeric {
        &algns1[0]
    } else {
        &algns2[0]
    };
    let chim3 = if first_read_is_chimeric {
        &algns1[1]
    } else {
        &algns2[1]
    };
    let linear_algn = if first_read_is_chimeric {
        &algns2[0]
    } else {
        &algns1[0]
    };

    if !(linear_algn.is_mapped && linear_algn.is_unique) {
        return None;
    }

    let can_rescue = if chim3.is_mapped && chim5.is_unique {
        let mut ok = chim3.chrom == linear_algn.chrom;
        ok &= chim3.strand != linear_algn.strand;

        if ok {
            if linear_algn.strand == '+' {
                ok &= linear_algn.pos5 < chim3.pos5;
            } else {
                ok &= linear_algn.pos5 > chim3.pos5;
            }
        }

        if ok {
            let molecule_size = if linear_algn.strand == '+' {
                chim3.pos5.saturating_sub(linear_algn.pos5)
                    + chim3.dist_to_5
                    + linear_algn.dist_to_5
            } else {
                linear_algn.pos5.saturating_sub(chim3.pos5)
                    + chim3.dist_to_5
                    + linear_algn.dist_to_5
            };
            ok &= molecule_size <= max_molecule_size;
        }

        ok
    } else {
        true
    };

    if can_rescue {
        Some(if first_read_is_chimeric { 1 } else { 2 })
    } else {
        None
    }
}

fn check_pair_order(
    algn1: &Alignment,
    algn2: &Alignment,
    chrom_enum: &HashMap<String, usize>,
) -> bool {
    let has_correct_order = (algn1.is_mapped as u8, algn1.is_unique as u8)
        <= (algn2.is_mapped as u8, algn2.is_unique as u8);

    if algn1.chrom != UNMAPPED_CHROM && algn2.chrom != UNMAPPED_CHROM {
        let order1 = chrom_enum.get(&algn1.chrom).copied().unwrap_or(usize::MAX);
        let order2 = chrom_enum.get(&algn2.chrom).copied().unwrap_or(usize::MAX);
        return (order1, algn1.pos) <= (order2, algn2.pos);
    }

    has_correct_order
}

fn parse_complex_walk<'a>(
    _algns1: &'a [Alignment],
    _algns2: &'a [Alignment],
    _max_insert_size: u32,
    _report_position: &'a str,
    _report_orientation: &'a str,
) -> Result<Box<dyn Iterator<Item = (Alignment, Alignment, (u32, &'static str))> + 'a>> {
    Ok(Box::new(std::iter::empty()))
}
