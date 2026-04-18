pub const PAIRSAM_FORMAT_VERSION: &str = "1.0.0";

pub const PAIRSAM_SEP: char = '\t';
pub const PAIRSAM_SEP_STR: &str = "\t";
pub const SAM_SEP: char = '\x1d';
pub const SAM_SEP_STR: &str = "\x1d";
pub const INTER_SAM_SEP: &str = "\x1dNEXT_SAM\x1d";

pub const COL_READID: usize = 0;
pub const COL_C1: usize = 1;
pub const COL_P1: usize = 2;
pub const COL_C2: usize = 3;
pub const COL_P2: usize = 4;
pub const COL_S1: usize = 5;
pub const COL_S2: usize = 6;
pub const COL_PTYPE: usize = 7;
pub const COL_SAM1: usize = 8;
pub const COL_SAM2: usize = 9;

pub const UNMAPPED_CHROM: &str = "!";
pub const UNMAPPED_POS: u32 = 0;
pub const UNMAPPED_STRAND: &str = "-";

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PairType {
    UU,
    UM,
    MU,
    MM,
    UR,
    RU,
    NU,
    UN,
    MN,
    NM,
    NN,
    UW,
    WU,
    WW,
    XX,
    DD,
}

impl PairType {
    pub fn from_chars(side1: char, side2: char) -> Self {
        match (side1, side2) {
            ('U', 'U') => Self::UU,
            ('U', 'M') => Self::UM,
            ('M', 'U') => Self::MU,
            ('M', 'M') => Self::MM,
            ('U', 'R') => Self::UR,
            ('R', 'U') => Self::RU,
            ('U', 'N') => Self::NU,
            ('N', 'U') => Self::UN,
            ('M', 'N') => Self::MN,
            ('N', 'M') => Self::NM,
            ('N', 'N') => Self::NN,
            ('U', 'W') => Self::UW,
            ('W', 'U') => Self::WU,
            ('W', 'W') => Self::WW,
            ('X', 'X') => Self::XX,
            ('D', 'D') => Self::DD,
            _ => panic!("Invalid pair type chars: {}{}", side1, side2),
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            Self::UU => "UU".to_string(),
            Self::UM => "UM".to_string(),
            Self::MU => "MU".to_string(),
            Self::MM => "MM".to_string(),
            Self::UR => "UR".to_string(),
            Self::RU => "RU".to_string(),
            Self::NU => "NU".to_string(),
            Self::UN => "UN".to_string(),
            Self::MN => "MN".to_string(),
            Self::NM => "NM".to_string(),
            Self::NN => "NN".to_string(),
            Self::UW => "UW".to_string(),
            Self::WU => "WU".to_string(),
            Self::WW => "WW".to_string(),
            Self::XX => "XX".to_string(),
            Self::DD => "DD".to_string(),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::UU => "UU",
            Self::UM => "UM",
            Self::MU => "MU",
            Self::MM => "MM",
            Self::UR => "UR",
            Self::RU => "RU",
            Self::NU => "NU",
            Self::UN => "UN",
            Self::MN => "MN",
            Self::NM => "NM",
            Self::NN => "NN",
            Self::UW => "UW",
            Self::WU => "WU",
            Self::WW => "WW",
            Self::XX => "XX",
            Self::DD => "DD",
        }
    }

    pub fn is_mapped(&self) -> bool {
        !matches!(self, Self::NN | Self::XX)
    }
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub chrom: String,
    pub pos5: u32,
    pub pos3: u32,
    pub pos: u32,
    pub strand: char,
    pub mapq: u8,
    pub is_mapped: bool,
    pub is_unique: bool,
    pub is_linear: bool,
    pub dist_to_5: u32,
    pub dist_to_3: u32,
    pub cigar: String,
    pub algn_ref_span: u32,
    pub algn_read_span: u32,
    pub matched_bp: u32,
    pub clip5_ref: u32,
    pub clip3_ref: u32,
    pub read_len: u32,
    pub pair_type_char: char,
    pub mismatches: String,
    pub sam_tags: Vec<(String, String)>,
    pub read_side: i32,
    pub algn_idx: i32,
    pub same_side_algn_count: i32,
}

impl Default for Alignment {
    fn default() -> Self {
        Self {
            chrom: UNMAPPED_CHROM.to_string(),
            pos5: UNMAPPED_POS,
            pos3: UNMAPPED_POS,
            pos: UNMAPPED_POS,
            strand: UNMAPPED_STRAND.chars().next().unwrap_or('-'),
            mapq: 0,
            is_mapped: false,
            is_unique: false,
            is_linear: true,
            dist_to_5: 0,
            dist_to_3: 0,
            cigar: "*".to_string(),
            algn_ref_span: 0,
            algn_read_span: 0,
            matched_bp: 0,
            clip5_ref: 0,
            clip3_ref: 0,
            read_len: 0,
            pair_type_char: 'N',
            mismatches: String::new(),
            sam_tags: Vec::new(),
            read_side: 0,
            algn_idx: 0,
            same_side_algn_count: 0,
        }
    }
}

impl Alignment {
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn mask(&mut self) {
        self.chrom = UNMAPPED_CHROM.to_string();
        self.pos5 = UNMAPPED_POS;
        self.pos3 = UNMAPPED_POS;
        self.pos = UNMAPPED_POS;
        self.strand = UNMAPPED_STRAND.chars().next().unwrap_or('-');
    }

    pub fn flip(&mut self) {
        std::mem::swap(&mut self.pos5, &mut self.pos3);
        self.strand = if self.strand == '+' { '-' } else { '+' };
    }

    pub fn flip_orientation(&mut self) {
        self.strand = if self.strand == '+' { '-' } else { '+' };
    }

    pub fn flip_position(&mut self) {
        std::mem::swap(&mut self.pos5, &mut self.pos3);
    }

    pub fn type_str(&self) -> char {
        if !self.is_mapped {
            'N'
        } else if !self.is_unique {
            'M'
        } else {
            'U'
        }
    }
}

#[derive(Debug, Clone)]
pub struct PairRecord {
    pub read_id: String,
    pub chrom1: String,
    pub pos1: u32,
    pub chrom2: String,
    pub pos2: u32,
    pub strand1: char,
    pub strand2: char,
    pub pair_type: PairType,
    pub sam1: Option<String>,
    pub sam2: Option<String>,
    pub pair_index: Option<u32>,
    pub walk_pair_type: Option<String>,
    pub extra_columns: Vec<(String, String)>,
}

impl PairRecord {
    pub fn new() -> Self {
        Self {
            read_id: String::new(),
            chrom1: UNMAPPED_CHROM.to_string(),
            pos1: 0,
            chrom2: UNMAPPED_CHROM.to_string(),
            pos2: 0,
            strand1: '-',
            strand2: '-',
            pair_type: PairType::NN,
            sam1: None,
            sam2: None,
            pair_index: None,
            walk_pair_type: None,
            extra_columns: Vec::new(),
        }
    }

    pub fn to_pairsam_line(
        &self,
        drop_readid: bool,
        drop_sam: bool,
        add_pair_index: bool,
    ) -> String {
        let mut fields = Vec::new();

        if !drop_readid {
            fields.push(self.read_id.clone());
        }
        fields.push(self.chrom1.clone());
        fields.push(self.pos1.to_string());
        fields.push(self.chrom2.clone());
        fields.push(self.pos2.to_string());
        fields.push(self.strand1.to_string());
        fields.push(self.strand2.to_string());
        fields.push(self.pair_type.to_string());

        if !drop_sam {
            if let Some(ref sam1) = self.sam1 {
                fields.push(sam1.clone());
            } else {
                fields.push(String::new());
            }
            if let Some(ref sam2) = self.sam2 {
                fields.push(sam2.clone());
            } else {
                fields.push(String::new());
            }
        }

        if add_pair_index {
            if let Some(idx) = self.pair_index {
                fields.push(idx.to_string());
            }
            if let Some(ref wpt) = self.walk_pair_type {
                fields.push(wpt.clone());
            }
        }

        for (k, v) in &self.extra_columns {
            fields.push(format!("{}={}", k, v));
        }

        fields.join(PAIRSAM_SEP_STR)
    }
}

#[derive(Debug, Clone)]
pub struct PairsamHeader {
    pub format_version: String,
    pub shape: String,
    pub genome_assembly: String,
    pub chromsizes: Vec<(String, u64)>,
    pub columns: Vec<String>,
    pub samheader_lines: Vec<String>,
    pub sorted: Option<String>,
    pub extra_fields: Vec<String>,
}

impl PairsamHeader {
    pub fn new() -> Self {
        Self {
            format_version: PAIRSAM_FORMAT_VERSION.to_string(),
            shape: "upper triangle".to_string(),
            genome_assembly: "unknown".to_string(),
            chromsizes: Vec::new(),
            columns: vec![
                "readID".to_string(),
                "chrom1".to_string(),
                "pos1".to_string(),
                "chrom2".to_string(),
                "pos2".to_string(),
                "strand1".to_string(),
                "strand2".to_string(),
                "pair_type".to_string(),
                "sam1".to_string(),
                "sam2".to_string(),
            ],
            samheader_lines: Vec::new(),
            sorted: None,
            extra_fields: Vec::new(),
        }
    }

    pub fn to_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();

        lines.push(format!("## pairs format v{}", self.format_version));
        lines.push(format!("#shape: {}", self.shape));
        lines.push(format!("#genome_assembly: {}", self.genome_assembly));

        for (chrom, size) in &self.chromsizes {
            lines.push(format!("#chromsize: {} {}", chrom, size));
        }

        if let Some(ref sorted) = self.sorted {
            lines.push(format!("#sorted: {}", sorted));
        }

        for line in &self.samheader_lines {
            lines.push(format!("#samheader: {}", line));
        }

        for field in &self.extra_fields {
            lines.push(field.clone());
        }

        lines.push(format!("#columns: {}", self.columns.join(" ")));

        lines
    }

    pub fn parse_header_lines(lines: &[String]) -> (Self, usize) {
        let mut header = Self::new();
        let mut header_end = 0;

        for (i, line) in lines.iter().enumerate() {
            if !line.starts_with('#') {
                break;
            }
            header_end = i + 1;

            if line.starts_with("##") {
                if let Some(ver) = line.strip_prefix("## pairs format v") {
                    header.format_version = ver.trim().to_string();
                }
                continue;
            }

            if let Some(rest) = line.strip_prefix("#shape: ") {
                header.shape = rest.trim().to_string();
            } else if let Some(rest) = line.strip_prefix("#genome_assembly: ") {
                header.genome_assembly = rest.trim().to_string();
            } else if let Some(rest) = line.strip_prefix("#chromsize: ") {
                let parts: Vec<&str> = rest.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(size) = parts[1].parse::<u64>() {
                        header.chromsizes.push((parts[0].to_string(), size));
                    }
                }
            } else if let Some(rest) = line.strip_prefix("#columns: ") {
                header.columns = rest.split_whitespace().map(|s| s.to_string()).collect();
            } else if let Some(rest) = line.strip_prefix("#sorted: ") {
                header.sorted = Some(rest.trim().to_string());
            } else if let Some(rest) = line.strip_prefix("#samheader: ") {
                header.samheader_lines.push(rest.trim().to_string());
            } else {
                header.extra_fields.push(line.clone());
            }
        }

        (header, header_end)
    }
}
