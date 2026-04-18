use crate::dedup::{DedupBackend, DedupConfig, DedupMethod};
use crate::parse::{ParseConfig, WalksPolicy};
use crate::sort::SortConfig;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "pairtools")]
#[command(about = "Fast Rust implementation of pairtools for Hi-C data processing", long_about = None)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
}

#[derive(Subcommand)]
enum Commands {
    Parse {
        #[arg(help = "Input SAM/BAM file path")]
        sam_path: String,

        #[arg(help = "Chromosome order file (chrom.sizes)")]
        chroms_path: String,

        #[arg(short = 'o', long, default_value = "", help = "Output file")]
        output: String,

        #[arg(long, default_value = "1", help = "Min MAPQ for unique mapping")]
        min_mapq: u8,

        #[arg(
            long,
            default_value = "750",
            help = "Max molecule size for walk rescue"
        )]
        max_molecule_size: u32,

        #[arg(long, default_value = "20", help = "Max inter-alignment gap")]
        max_inter_align_gap: u32,

        #[arg(
            long,
            default_value = "5unique",
            help = "Walks policy (mask/5any/5unique/3any/3unique/all)"
        )]
        walks_policy: String,

        #[arg(long, help = "Drop read IDs from output")]
        drop_readid: bool,

        #[arg(long, help = "Drop sequences from SAM fields")]
        drop_seq: bool,

        #[arg(long, help = "Drop SAM fields from output")]
        drop_sam: bool,

        #[arg(long, help = "Add pair index to output")]
        add_pair_index: bool,

        #[arg(long, default_value = "5", help = "Report 5' or 3' alignment end")]
        report_alignment_end: String,

        #[arg(long, default_value_t = true, help = "Flip to upper triangular order")]
        flip: bool,
    },
    Sort {
        #[arg(help = "Input pairs/pairsam file")]
        pairs_path: Option<String>,

        #[arg(short = 'o', long, default_value = "", help = "Output file")]
        output: String,

        #[arg(long, default_value = "8", help = "Number of threads for sorting")]
        nproc: usize,

        #[arg(long, default_value = "2G", help = "Memory limit (e.g. 2G, 4G)")]
        memory: String,

        #[arg(long, default_value = "chrom1", help = "Chrom1 column name/index")]
        c1: String,

        #[arg(long, default_value = "chrom2", help = "Chrom2 column name/index")]
        c2: String,

        #[arg(long, default_value = "pos1", help = "Pos1 column name/index")]
        p1: String,

        #[arg(long, default_value = "pos2", help = "Pos2 column name/index")]
        p2: String,

        #[arg(
            long,
            default_value = "pair_type",
            help = "Pair type column name/index"
        )]
        pt: String,
    },
    Dedup {
        #[arg(help = "Input sorted pairs/pairsam file")]
        pairs_path: Option<String>,

        #[arg(
            short = 'o',
            long,
            default_value = "",
            help = "Output deduplicated file"
        )]
        output: String,

        #[arg(long, help = "Output duplicates file")]
        output_dups: Option<String>,

        #[arg(long, default_value = "3", help = "Max mismatch distance (bp)")]
        max_mismatch: u32,

        #[arg(long, default_value = "max", help = "Dedup method (max or sum)")]
        method: String,

        #[arg(long, default_value = "scipy", help = "Backend (scipy/cython)")]
        backend: String,

        #[arg(
            long,
            default_value = "10000",
            help = "Chunk size for scipy/sklearn backend"
        )]
        chunksize: usize,

        #[arg(long, default_value = "100", help = "Carryover between chunks")]
        carryover: usize,

        #[arg(long, default_value_t = true, help = "Mark duplicates as DD")]
        mark_dups: bool,

        #[arg(long, help = "Keep parent read ID for dups")]
        keep_parent_id: bool,
    },
}

pub fn run() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.verbose {
        0 => {}
        1 => std::env::set_var("RUST_LOG", "info"),
        _ => std::env::set_var("RUST_LOG", "debug"),
    }
    env_logger::init();

    match cli.command {
        Commands::Parse {
            sam_path,
            chroms_path,
            output,
            min_mapq,
            max_molecule_size,
            max_inter_align_gap,
            walks_policy,
            drop_readid,
            drop_seq,
            drop_sam,
            add_pair_index,
            report_alignment_end,
            flip,
        } => {
            let wp = match walks_policy.as_str() {
                "mask" => WalksPolicy::Mask,
                "5any" => WalksPolicy::FiveAny,
                "5unique" => WalksPolicy::FiveUnique,
                "3any" => WalksPolicy::ThreeAny,
                "3unique" => WalksPolicy::ThreeUnique,
                "all" => WalksPolicy::All,
                _ => anyhow::bail!("Unknown walks policy: {}", walks_policy),
            };

            let config = ParseConfig {
                min_mapq,
                max_molecule_size,
                max_inter_align_gap,
                walks_policy: wp,
                drop_readid,
                drop_seq,
                drop_sam,
                add_pair_index,
                report_alignment_end: report_alignment_end.chars().next().unwrap_or('5'),
                flip,
                ..Default::default()
            };
            crate::parse::streaming_classify(&sam_path, &chroms_path, &output, &config)?;
        }
        Commands::Sort {
            pairs_path,
            output,
            nproc,
            memory,
            c1,
            c2,
            p1,
            p2,
            pt,
        } => {
            let input = pairs_path.unwrap_or_else(|| "-".to_string());
            let mem_bytes = parse_memory_size(&memory);
            let config = SortConfig {
                nproc,
                memory: mem_bytes,
                c1_col: c1,
                c2_col: c2,
                p1_col: p1,
                p2_col: p2,
                pt_col: pt,
                ..Default::default()
            };
            crate::sort::sort_pairs(&input, &output, &config)?;
        }
        Commands::Dedup {
            pairs_path,
            output,
            output_dups,
            max_mismatch,
            method,
            backend,
            chunksize,
            carryover,
            mark_dups,
            keep_parent_id,
        } => {
            let input = pairs_path.unwrap_or_else(|| "-".to_string());
            let method = if method == "sum" {
                DedupMethod::Sum
            } else {
                DedupMethod::Max
            };
            let backend = if backend == "cython" {
                DedupBackend::Cython
            } else {
                DedupBackend::Scipy
            };
            let config = DedupConfig {
                max_mismatch,
                method,
                mark_dups,
                keep_parent_id,
                backend,
                chunksize,
                carryover,
                ..Default::default()
            };
            crate::dedup::streaming_dedup(&input, &output, output_dups.as_deref(), &config)?;
        }
    }

    Ok(())
}

fn parse_memory_size(s: &str) -> usize {
    let s = s.trim().to_uppercase();
    let (num_str, suffix) = if s.ends_with('G') {
        (&s[..s.len() - 1], 1024u64 * 1024 * 1024)
    } else if s.ends_with('M') {
        (&s[..s.len() - 1], 1024u64 * 1024)
    } else if s.ends_with('K') {
        (&s[..s.len() - 1], 1024)
    } else {
        (s.as_str(), 1)
    };
    num_str.parse::<f64>().unwrap_or(2.0) as usize * suffix as usize
}
