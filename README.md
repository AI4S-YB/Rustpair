# Rustpair
A Rust rewrite of pairtools for processing paired-end Hi-C sequencing data files.

https://img.shields.io/badge/rust-1.75+-orange.svg
https://img.shields.io/badge/license-MIT-blue.svg

📊 Comparison with Python Version
Feature	Python pairtools	Rustpair
Language	Python 3 / Cython	Rust
SAM/BAM Parsing	pysam (libhts)	rust-htslib
Sorting	Calls Unix sort	rayon parallel + k-way merge
Deduplication	scipy KD-tree / Cython	Online sliding window
Multithreading	Limited	Fully parallel with rayon
Memory Safety	GC-managed	Compile-time guarantees
🚀 Installation
📦 Building from Source
bash
# Requires Rust toolchain (>= 1.75)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

git clone <repo-url>
cd pairtools-rs
cargo build --release

# Binary will be at target/release/pairtools
🔧 System Dependencies
<details> <summary><b>Ubuntu / Debian</b></summary>
bash
sudo apt install libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
</details><details> <summary><b>macOS</b></summary>
bash
brew install zlib xz
</details>
📖 Usage

text
pairtools [OPTIONS] <COMMAND>

Commands:
  parse    Convert SAM/BAM files to pairs/pairsam format
  sort     Sort pairs/pairsam files
  dedup    Remove PCR duplicates
  help     Print help information

Options:
  -v, --verbose...  Increase logging verbosity (-v info, -vv debug)
  -h, --help        Print help
  -V, --version     Print version
🔄 parse — SAM/BAM to pairsam Conversion
Convert aligned SAM/BAM files to standard pairs format, automatically detecting pairing relationships, chimeric reads, and duplicates.

Syntax
bash
pairtools parse <SAM_PATH> <CHROMS_PATH> [OPTIONS]
Required Arguments
Argument	Description
SAM_PATH	Input SAM or BAM file path (- for stdin)
CHROMS_PATH	Chromosome sizes file (chrom.sizes), format: chr_name length per line
Common Options
Option	Default	Description
-o, --output	"" (stdout)	Output file path
--min-mapq	1	Minimum MAPQ threshold; alignments below this are considered non-unique (M)
--walks-policy	5unique	Chimeric walks handling policy (see below)
--max-molecule-size	750	Maximum molecule length for walk rescue (bp)
--max-inter-align-gap	20	Maximum gap between adjacent alignments (bp)
--flip	enabled	Flip to upper-triangular coordinates (chrom1 <= chrom2, pos1 <= pos2)
--report-alignment-end	5	Report 5' or 3' end position
--drop-readid	disabled	Drop read ID column from output
--drop-sam	disabled	Drop SAM field columns from output
--drop-seq	disabled	Drop sequence information
Walks Policy Details
Policy	Behavior
mask	Mark chimeric pairs as WW (masked)
5any	Take first alignment (5' end)
5unique	Take first unique alignment at 5' end (recommended)
3any	Take last alignment (3' end)
3unique	Take last unique alignment at 3' end
all	Expand all walk combinations into multiple rows
Examples
bash
# Basic usage
pairtools parse aligned.sam chrom.sizes -o output.pairsam

# Read BAM from stdin
samtools view -h input.bam | pairtools parse - chrom.sizes -o output.pairsam

# Drop SAM fields, output pure pairs format
pairtools parse input.bam chrom.sizes --drop-sam -o output.pairs

# Use 3' unique policy for chimeric reads
pairtools parse input.bam chrom.sizes --walks-policy 3unique -o output.pairsam
Output Format (pairsam)
text
## pairs format v1.0.0
#shape: upper triangle
#chromsize: chr1 248956422
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1 sam2
read001	chr1	1000	chr5	50000	+	-	UU	read001\x1d...	read001\x1d...
Pair Type Codes
Code	Meaning
UU	Both ends uniquely mapped
UM / MU	One unique, one multi-mapped
MM	Both ends multi-mapped
UR / RU	One unique, one rescued
NU / UN	One unmapped, one unique
NM / MN	One unmapped, one multi-mapped
NN	Both ends unmapped
WW	Masked chimeric pair
XX	Empty/invalid pair
DD	PCR duplicate (marked after dedup)
📊 sort — Sorting
Sort pairs/pairsam files by chromosome and position, with multi-threaded parallel support.

Syntax
bash
pairtools sort [PAIRS_PATH] [OPTIONS]
Options
Option	Default	Description
-o, --output	"" (stdout)	Output file path
--nproc	8	Number of sorting threads
--memory	2G	Memory limit (e.g., 4G, 500M)
--c1	chrom1	Column name/index for first chromosome
--c2	chrom2	Column name/index for second chromosome
--p1	pos1	Column name/index for first position
--p2	pos2	Column name/index for second position
--pt	pair_type	Column name/index for pair type
Examples
bash
# Basic sorting
pairtools sort input.pairsam -o sorted.pairsam

# Use 16 threads, 4GB memory
pairtools sort input.pairsam -o sorted.pairsam --nproc 16 --memory 4G

# Pipeline usage
pairtools parse ... | pairtools sort - -o sorted.pairsam
Algorithm
Data is split into chunks and sorted in parallel using rayon

Sorted chunks are merged via k-way merge (min-heap)

Sort key: (chrom1, chrom2, pos1, pos2, pair_type)

🔍 dedup — Deduplication
Remove PCR/optical duplicates using an online sliding window algorithm, without loading all data into memory.

Syntax
bash
pairtools dedup [PAIRS_PATH] [OPTIONS]
⚠️ Note: Input must be pre-sorted! Typically used after sort.

Options
Option	Default	Description
-o, --output	"" (stdout)	Deduplicated output file
--output-dups	none	Separate output file for reads marked as duplicates
--max-mismatch	3	Maximum positional tolerance for duplicate detection (bp)
--method	max	Distance calculation: max (max difference) or sum (sum of differences)
--mark-dups	enabled	Mark duplicates as DD instead of removing them
--keep-parent-id	disabled	Append parent read ID to duplicate rows
Examples
bash
# Basic deduplication
pairtools sort input.pairsam -o sorted.pairsam && \
pairtools dedup sorted.pairsam -o deduped.pairsam

# Separate duplicate reads to a different file
pairtools dedup sorted.pairsam -o deduped.pairsam --output-dups dups.pairsam

# Relax tolerance to 10bp
pairtools dedup sorted.pairsam -o deduped.pairsam --max-mismatch 10
Algorithm
Maintains a sliding window, comparing only nearby candidate pairs

Time complexity ≈ O(n), memory usage O(window_size)

Pairs on the same chromosome, same strand, and within max_mismatch position are considered duplicates

🔗 Typical Workflow
Step-by-step
bash
# Complete Hi-C data processing pipeline
pairtools parse aligned.bam chrom.sizes \
  --drop-sam \
  -o raw.pairs && \
pairtools sort raw.pairs \
  --nproc 16 \
  -o sorted.pairs && \
pairtools dedup sorted.pairs \
  --output-dups dups.pairs \
  -o final.pairs
One-liner Pipeline
bash
pairtools parse aligned.bam chrom.sizes --drop-sam -o - | \
pairtools sort - --nproc 16 -o - | \
pairtools dedup - -o final.pairs --output-dups dups.pairs
⚡ Performance Recommendations
Command	Bottleneck	Recommendation
parse	I/O	Use SSD/NVMe for significantly improved throughput
sort	CPU + Memory	Set --nproc to number of cores; adjust --memory based on available RAM
dedup	Memory	Larger --max-mismatch increases window size; 3-5bp is usually sufficient
💡 Tip: Release builds enable LTO and opt-level=3. Always use cargo build --release for production.
