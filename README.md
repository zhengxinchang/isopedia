<img src="./img/logo2.png" alt="isopedia" width="150" align="right" style="vertical-align: middle;z-index: 1000;"> 


# About the Isopedia

Isopedia is a scalable tool that evaluates novel isoforms by leveraging population-scale long-read transcriptome data to distinguish true biological variants from technical artifacts.

# Quick Q&A:

**Q1: What gap does Isopedia aim to fill?**
A: Long-read RNA sequencing often reveals a large number of novel isoforms. However, evaluating whether these novel isoforms are biologically meaningful—or simply artifacts caused by RNA degradation or sequencing errors—is challenging. Understanding both the existence and population frequency of a novel isoform is essential for downstream analyses, yet current approaches are limited in scalability and robustness.

**Q2: How does Isopedia address this problem?**
A: Rather than relying on model-based classification or paired sequencing datasets from the same sample, Isopedia takes a different approach. It searches for supporting evidence of novel isoforms across large-scale long-read transcriptome datasets, leveraging population-level data to distinguish true biological isoforms from noise.

**Q3: How does Isopedia manage large datasets and provide a practical solution for isoform assessment?**
A: Isopedia introduces several innovations for scalable and efficient isoform evaluation. It uses a B+ tree–based data structure to rapidly index and compare isoform-related genomic positions. In addition, it employs a read-level signal extraction algorithm to build a compact yet informative index of transcriptome data. When assessing a query isoform, Isopedia integrates evidence from splicing junctions and alignment quality to ensure a robust evaluation. The entire tool is implemented in Rust, offering high performance and a user-friendly interface.

**Q4: What is the best use case for Isopedia?**
A: Isopedia is designed to be integrated into standard long-read transcriptome analysis pipelines. After isoforms are identified using tools like IsoQuant, FLAMES, TALON, or others, their output GTF files can be passed to Isopedia. It then evaluates the isoforms against a background index constructed from hundreds of publicly available long-read transcriptome datasets, providing a scalable and population-aware assessment.


# Quick Start

isopeida consists of multiple binaries that have prefix isopedia-*. This naming strategy help isopedia update each command individually and easily to expand.


**Download prebuild index and run**
```
# download index and uncompress it
wget xxx -O index.tar.gz
tar xzvf index.tar.gz

sopedia-anno-isoform -i lr_idx/ -g query.gtf -o isoform.anno.tsv

isopedia-anno-fusion -i lr_idx/ -p chr1:181130,chr1:201853853 -o fusion.anno.tsv
```

**Build your own index**

```
# extract isoform signals on each bam individually
isopedia-extr -b in.bam -o out.isoform.gz

# make a manifest.txt(tab-seprated) of isoform.gz files such as:
# path                               sample
# path/to/pb.isoform.gz              pb-bam
# path/to/ont.isoform.gz             ont-bam

# aggregate
isopedia-aggr  -i manifest.txt -o lr_idx/

# build index
isopedia-idx  -i lr_idx/

# annoate isoform by provide a GTF file
isopedia-anno-isoform -i lr_idx/ -g query.gtf -o isoform.anno.tsv

# annotate fusion by provide breapoint string
isopedia-anno-fusion -i lr_idx/ -p chr1:181130,chr1:201853853 -o fusion.anno.tsv
```


# How it works


![how-it-works](./img/how-it-works.png)



# Usage

```
Usage: isopedia <COMMAND>

Commands:
  extract  Extract isoform signals from Alignment file(BAM/CRAM)
  aggr     Aggregate and merge signals from multiple samples into a unified file
  idx      Create searchable indices
  search   Search and retrieve signals across indexed samples
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```


There are four main steps to use isopedia:

1. Extract isoform raw signal from long-read full-length transcriptome datasets
2. Aggregate isoform raw signal from multiple samples for building index
3. Build index on aggregated data
4. Annotate the population frequency of isoform in the population

## Extract singals

Extract isoform signals from Alignment file(BAM/CRAM)

```
Usage: isopedia extract [OPTIONS] --bam <BAM> --output <OUTPUT>

Options:
  -b, --bam <BAM>              Input file in BAM/CRAM format
  -r, --reference <REFERENCE>  Reference file for CRAMs
  -o, --output <OUTPUT>        Output
      --mapq <MAPQ>            Minimal mapping quality of reads [default: 5]
  -h, --help                   Print help

```
Please note that if the input file is in CRAM format, you need to provide the reference file for CRAMs.

**Example**

isopedia extract -b sampleA.mapped.bam -o sampleA.isoform.raw.txt


## Aggregate signals

Aggregate and merge signals from multiple samples into a unified file

```
Usage: isopedia aggr --input <INPUT> --outdir <OUTDIR>

Options:
  -i, --input <INPUT>
          Input manifest(tabs-separated) file.
          Frist 2 cols are required: sample path and sample name.
          The rest of cols are optional and will be used as sample meta.
          First row is header.

  -o, --outdir <OUTDIR>
          Output index directory

  -h, --help
          Print help (see a summary with '-h')

```
**Example**

The input file example(tab-separated)
(reqired)            (reqired)  (optional) (optional) ...
-------------------- ---------- ---------- ----------
path                 name       meta1      meta2     <- header
/path/to/sample1.bam sample1    value1     value2
/path/to/sample2.bam sample2    value1     value2

Run aggregation:

isopedia aggr --input input.tsv --outdir index_dir


## Build index

Create searchable indices

```
Usage: isopedia idx --idxdir <IDXDIR>

Options:
  -i, --idxdir <IDXDIR>  
  -h, --help             Print help
```

**Example**

isopedia idx -i idx/




## Search 

Search and retrieve signals across indexed samples

```
Usage: isopedia search [OPTIONS] --idxdir <IDXDIR>

Options:
  -i, --idxdir <IDXDIR>      index directory
  -p, --pos <POS>            positions to be search(-p chr:pos * N)
  -g, --gtf <GTF>            gtf file
  -f, --flank <FLANK>        flank size for search, before and after the position [default: 2]
  -m, --min-read <MIN_READ>  minimal reads to define a positive sample [default: 1]
  -o, --output <OUTPUT>      output file for search results
  -h, --help                 Print help
```



**Example**

isopedia search -i idx/ -g gencode.v47.basic.annotation.gtf -o gencode.report.txt


# Output

The output of the search command is a tab-separated file with the following columns:

| Column name                  | Description                              |
|------------------------------|------------------------------------------|
| chrom                        | Chromosome                               |
| start                        | Start position                           |
| end                          | End position                             |
| trans_id                     | Transcript ID                            |
| gene_id                      | Gene ID                                  |
| hit                          | How many samples support this transcript |
| min_read                     | Minimal reads to define a positive sample|
| positive_count/sample_size   | Positive count/sample size               |
| sample1                      | Evidence of sample1 for this transcript  |
| .......                      | ......                                   |
| sampleN                      | Evidence of sampleN for this transcript  |

# Installation

## Check out the latest release


https://github.com/zhengxinchang/isopedia/releases


## From source code

git clone https://github.com/zhengxinchang/isopedia.git
cd isopedia
cargo build --release



# Roadmap