<img src="./img/logo2.png" alt="isopedia" width="200"/>


# isopedia

isopedia is a tool that extends the concept of [STIX](https://github.com/ryanlayer/stix) to isoform analysis. It can extract isoform raw signal from long-read full-length transcriptome datasets and build index on it. isopedia provides a fast approach to annotate the population frequency of isoforms in the population.


# Quick Start




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