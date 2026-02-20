# About Isopedia <img src="./img/logo3.png" align="right" alt="" width=120 />

**Isopedia** is a scalable tool for analyzing hundreds to thousands of long-read transcriptome datasets using a read-level indexing approach.

It provides two key capabilities:

- Population-level transcript quantification and frequency profiling.
- Isoform diversity exploration, including fusion and splice-junction analysis.

## Table of Contents

- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Download Pre-built Index](#download-pre-built-index)
- [Build Your Own Index](#build-your-own-index)
- [How to Use Isopedia](#how-to-use-isopedia)
  - [Isoform Query](#isoform-query)
  - [Fusion Query and Discovery](#fusion-query-and-discovery)
  - [Splice Query and Visualization](#splice-query-and-visualization)
- [Command Parameters (Latest)](#command-parameters-latest)
- [How to Install Isopedia](#how-to-install-isopedia)
- [Annotate ORF with ORFannotate](#annotate-orf-with-orfannotate)
- [Computational Resource Usage](#computational-resource-usage)
- [Contact](#contact)

## Quick Start

Isopedia has two binaries:

- `isopedia`: main workflows (query + index building)
- `isopedia-tools`: helper utilities

```bash
# install isopedia from conda
conda install -c bioconda isopedia

# clone the repo and enter toy example directory
git clone https://github.com/zhengxinchang/isopedia && cd isopedia/toy_ex/

# query transcripts
isopedia isoform -i index/ -g query.gtf -o out.profile.tsv.gz

# query one fusion event (two breakpoints)
isopedia fusion -i index/ -p chr1:181130,chr1:201853853 -o out.fusion.tsv.gz

# query multiple fusion events
isopedia fusion -i index/ -P fusion_query.bed -o out.fusion.tsv.gz

# discover candidate fusions from gene regions
isopedia fusion -i index/ -G gene.gtf -o out.fusion.discovery.tsv.gz

# query one splice junction
isopedia splice -i index/ -s 22:17744013,22:17750104 -o out.splice.tsv.gz

# visualize splice query output
isopedia-splice-viz.py -i out.splice.tsv.gz -g gencode.v47.basic.chr22.gtf -o isopedia-splice-view
```

For GTF indexing details, see [Indexing GTF Files](doc/indexing_gtf.md).

## How It Works

![how-it-works](./img/how-it-works.png)

Typical workflow:

1. `isopedia profile`: extract isoform signals per sample.
2. `isopedia merge`: aggregate all profiles into one merged archive.
3. `isopedia index`: build tree index for fast query.
4. Query with `isopedia isoform`, `isopedia fusion`, or `isopedia splice`.

<details>
<summary><strong>How Isopedia determines a positive hit for a query transcript</strong></summary>

<img src="./img/how-it-works2.png" align="center" alt="" />

</details>

## Download Pre-built Index

Isopedia provides pre-built indexes from public long-read RNA-seq datasets.

| Name | Organism | Reference | Link | Sample Size | Index Size (Compressed) | Minimum Memory for Query | Description |
|---|---|---|---|---:|---|---|---|
| `isopedia_index_hs_v1.0` | *Homo sapiens* | GRCh38 | `ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz` | 1,007 | 305G (107G) | 16 GB | 107 ENCODE samples + 900 SRA samples |

Checksum file:

- `ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz.md5sum`

### Download with `isopedia download`

```bash
# list indexes from the default remote manifest
isopedia download --list

# download a named index from the default manifest
isopedia download --name isopedia_index_hs_v1.0.tar.xz --outdir ./downloads

# download from custom URL
isopedia download --url ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz --outdir ./downloads

# use custom manifest TOML
isopedia download --manifest /path/to/manifest.toml --name isopedia_index_hs_v1.0.tar.xz --outdir ./downloads

```

Example custom manifest entry:

```toml
[[index]]
name = "isopedia_index_hs_v1.0.tar.xz"
url = "ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz"
size = 114596438692
md5 = "eb922ef27257d363969a835d4175da26"
source = "ftp"
protocol = "ftp"
description = "Human pre-built index (GRCh38)"
```

Unpack:

```bash
tar -xvf isopedia_index_hs_v1.0.tar.xz
```

## Build Your Own Index

Prerequisites:

1. Long-read alignments (`.bam`/`.cram`) or transcript annotations (`.gtf`).
2. A manifest TSV with at least two columns: sample name and profile path.

Example manifest:

```tsv
name	path	platform
HG002_pb_chr22	/path/to/hg002_pb_chr22.isoform.gz	PacBio
HG002_ont_chr22	/path/to/hg002_ont_chr22.isoform.gz	ONT
```

Workflow:

```bash
# 1) profile each sample
isopedia profile -i ./chr22.pb.grch38.bam -o ./hg002_pb_chr22.isoform.gz
isopedia profile -i ./chr22.ont.grch38.bam -o ./hg002_ont_chr22.isoform.gz

# optional profile fields
# --rname (BAM/CRAM), --tid and --gid (GTF)

# 2) merge profiles
isopedia merge -i manifest.tsv -o index/

# 3) build tree index
isopedia index -i index/ -m manifest.tsv

# 4) test query
isopedia isoform -i index/ -g gencode.v47.basic.chr22.gtf -o out.isoform.tsv.gz
```

## How to Use Isopedia

### Isoform Query

Purpose:

- Annotate transcripts from an input GTF against the index.

Example:

```bash
# input GTF must be sorted
gffread -T -o- origin.gtf | sort -k1,1 -k4,4n | gffread - -o query.sorted.gtf

isopedia isoform -i index/ -g query.sorted.gtf -o isoform.out.tsv.gz
```

#### Output columns (latest code)

| Column | Description |
|---|---|
| `chrom` | Chromosome |
| `start` | Transcript start |
| `end` | Transcript end |
| `length` | Transcript length |
| `exon_count` | Exon count |
| `trans_id` | Transcript ID |
| `gene_id` | Gene ID |
| `confidence` | Confidence score across cohort |
| `detected(total:fsm:em)` | Detection status for total/FSM/EM |
| `min_read` | Minimum read threshold used |
| `n_pos_samples(total:fsm:em/sample_size)` | Positive sample counts and sample size |
| `attributes` | Original GTF attributes |
| `FORMAT` | Per-sample field format |
| `sample_*` | Per-sample values |

Current `FORMAT` string in output header:

- `CPM:COUNT:FSM_CPM:FSM_COUNT:EM_CPM:EM_COUNT:INFO`

Per-sample values include abundance/count components for total, FSM, and EM estimates.

### Fusion Query and Discovery

`fusion` supports three modes:

1. Breakpoint query (`--pos` or `--pos-bed`)
2. Gene-region discovery (`--gene-gtf`)
3. Region-pair detailed read output (`--region`)

Examples:

```bash
# single breakpoint pair
isopedia fusion -i index/ -p chr1:1000,chr2:2000 -o fusion.single.tsv.gz

# batch breakpoint pairs from bed-like file
isopedia fusion -i index/ -P fusion_breakpoints.bed -o fusion.batch.tsv.gz

# discover candidate fusions from gene GTF
isopedia fusion -i index/ -G gene.gtf -o fusion.discovery.tsv.gz

# detailed read-level output for two regions
isopedia fusion -i index/ -r chr8:92017266-92017466,chr21:34859374-34859574 -o fusion.region.tsv.gz
```

Breakpoint-query output columns:

- `chr1`, `pos1`, `chr2`, `pos2`, `id`, `min_read`, `positive/sample_size`, `left_isoforms`, `right_isoforms`, then per-sample counts.

Gene-discovery output columns:

- `geneA_name`, `geneB_name`, `total_evidences`, `total_samples`, `is_two_strand`,
  `AtoB_primary_start`, `AtoB_primary_end`, `AtoB_supp_start`, `AtoB_supp_end`,
  `BtoA_primary_start`, `BtoA_primary_end`, `BtoA_supp_start`, `BtoA_supp_end`, then per-sample counts.

Region-pair detailed output columns:

- `chr1`, `start1`, `end1`, `chr2`, `start2`, `end2`, `main_exon_count1`,
  `supp_segment_count2`, `query_part`, `main_isoforms`, `supp_aln_regions`, `sample_name`.

### Splice Query and Visualization

Purpose:

- Query isoforms overlapping a splice junction and optionally visualize them.

Examples:

```bash
# single splice query
isopedia splice -i index/ -s chr22:41100500,chr22:41101500 -o splice.out.tsv.gz

# batch splice query
isopedia splice -i index/ -S splice_queries.bed -o splice.batch.tsv.gz

# visualize
python script/isopedia-splice-viz.py \
  -i splice.out.tsv.gz \
  -g gencode.v47.basic.annotation.gtf \
  -t script/isopedia-splice-viz-temp.html \
  -o isopedia-splice-view
```

`splice` output columns:

- `id`, `chr1`, `pos1`, `chr2`, `pos2`, `total_evidence`, `cpm`, `matched_sj_idx`,
  `dist_to_matched_sj`, `n_exons`, `start_pos_left`, `start_pos_right`,
  `end_pos_left`, `end_pos_right`, `splice_junctions`, then per-sample values.

`splice` per-sample `FORMAT`:

- `COUNT:CPM:START|END|STRAND`

## Command Parameters (Latest)

Parameter lists below are based on the current CLI in source.

### `isopedia isoform`

Core options:

- `-i, --idxdir <IDXDIR>`
- `-g, --gtf <GTF>`
- `-o, --output <OUTPUT>`
- `-f, --flank <FLANK>` (default: `10`)
- `-m, --min-read <MIN_READ>` (default: `1`)
- `--info`
- `-n, --num-threads <NUM_THREADS>` (default: `4`)
- EM options: `--em-max-iter`, `--em-conv-min-diff`, `--em-chunk-size`, `--em-effective-len-coef`, `--em-damping-factor`, `--min-em-abundance`
- TSS/TES options: `--no-check-tss-tes`, `--tss-degrad-bp`, `--tes-degrad-bp`, `--terminal-tolerance-bp`
- Cache options: `-c, --cached-nodes`, `--cached-chunk-num`, `--cached-chunk-size-mb`
- `--output-tmp-shard-counts`
- `--verbose`

### `isopedia fusion`

Core options:

- `-i, --idxdir <IDXDIR>`
- Query selectors: `-p, --pos`, `-P, --pos-bed`, `-r, --region`, `-G, --gene-gtf`
- `-o, --output <OUTPUT>`
- `-f, --flank <FLANK>` (default: `10`)
- `-m, --min-read <MIN_READ>` (default: `1`)
- Cache options: `-c, --cached_nodes`, `--cached-chunk-number`, `--cached-chunk-size-mb`
- `--verbose`

### `isopedia splice`

Core options:

- `-i, --idxdir <IDXDIR>`
- Query selectors: `-s, --splice` or `-S, --splice-bed`
- `-o, --output <OUTPUT>`
- `-f, --flank <FLANK>` (default: `10`)
- `-m, --min-read <MIN_READ>` (default: `1`)
- Cache options: `-c, --cached_nodes`, `--cached-chunk-number`, `--cached-chunk-size-mb`
- `--verbose`

### `isopedia profile`

Core options:

- Input selectors: `-i, --bam` or `-g, --gtf`
- `-r, --reference <REFERENCE>` for CRAM
- `-o, --output <OUTPUT>`
- `--mapq <MAPQ>` (default: `5`)
- `--use-secondary`, `--rname`, `--tid`, `--gid`, `--verbose`

### `isopedia merge`

Core options:

- `-i, --input <INPUT>`
- `-o, --outdir <OUTDIR>`
- `-c, --chunk-size <CHUNK_SIZE>` (default: `1000000`)

### `isopedia index`

Core options:

- `-i, --idxdir <IDXDIR>`
- `-m, --manifest <MANIFEST>`
- `-t, --threads <THREADS>` (default: `4`)

### `isopedia download`

Core options:

- `-l, --list`
- `-n, --name <NAME>`
- `-u, --url <URL>`
- `-m, --manifest <MANIFEST>`
- `-o, --outdir <OUTDIR>` (default: `.`)

## How to Install Isopedia

### Install from conda

```bash
conda install -c zhengxinchang isopedia
```

### Latest release

- https://github.com/zhengxinchang/isopedia/releases

`isopedia-<version>.linux.tar.gz` is compiled on Amazon Linux 2 (GCC 7.3, glibc 2.26, binutils 2.29.1).
If your distribution is not compatible, use `isopedia-<version>.musl.tar.gz`.

### Build from source

Rust and Cargo are required.

```bash
git clone https://github.com/zhengxinchang/isopedia.git
cd isopedia
cargo build --release
cargo build --release --target x86_64-unknown-linux-musl
```

## Annotate ORF with ORFannotate

Isopedia interoperates with [ORFannotate](https://github.com/egustavsson/ORFannotate), which can consume Isopedia outputs to predict ORFs/UTRs and annotate CDS features.

```text
<placeholder>
```

## Computational Resource Usage

Isopedia 1.4.0 benchmark on 1,007 long-read transcriptome datasets from SRA and ENCODE.

- Hardware: AMD Ryzen 9 7940HX, 64 GB RAM.

| Step | Peak Mem (GB) | Time (H:MM:SS) |
|---|---:|---|
| `isopedia merge` | 54.77 | 28:26:08 |
| `isopedia index` | 45.85 | 5:48:55 |
| `isopedia isoform` (280K GENCODE v49 basic) | 9.19 | 4:52:18 |

## Contact

- zhengxc93@gmail.com
- fritz.sedlazeck@bcm.edu
