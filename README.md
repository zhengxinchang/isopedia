

# About the Isopedia <img src="./img/logo2.png" align="right" alt="" width=120 />

Isopedia is a scalable tool that evaluates novel isoforms by leveraging population-scale long-read transcriptome data to distinguish true biological variants from technical artifacts.

<details>
<summary>Table of Content</summary>


- [Quick Start](#quick-start)
- [How it works](#how-it-works)
- [Usage](#usage)
  - [Annotate (search) isoforms/transcripts](#annotate-search-isoformstranscripts)
  - [Annotate (search) fusion genes](#annotate-search-fusion-genes)
  - [Find potential fusion genes](#find-potential-fusion-genes)
- [Installation](#installation)
- [Quick Q&A](#quick-qa)

</details>

# Quick Start

isopeida consists of multiple binaries that have prefix isopedia-*. This naming strategy help isopedia update each command individually and easily to expand.


**Download prebuild index and run**
```bash
isopedia isoform -i index/ -g query.gtf -o isoform.anno.tsv

isopedia fusion  -i index/ -p chr1:181130,chr1:201853853 -o fusion.anno.tsv

isopedia fusion  -i index/ -P fusion_query.bed -o fusion.anno.tsv

isopedia fusion  -i index/ -g gene.gtf -o fusion.discovery.tsv
```

**Build your own index**


Isopedia supports building local index in your own datasets. prerequests are listed below:

1. latest isopedia binaries
2. A set of mapped bam files(sorted bam are not required)
3. A manifest file that describe the sample name, isoform file path, and other optional meta data in tabular(\t sperated and with a header line) format


<details>

<summary>
You can find example files and commands at here [click to expand]
</summary>

```bash
# make sure isopedia in your $PATH or use absolute path to the binaries.

# download the toy_ex 
git clone https://github.com/zhengxinchang/isopedia && cd isopedia/toy_ex/

# extract isoform signals on each bam individually
isopedia extr -b ./chr22.pb.grch38.bam -o ./hg002_pb_chr22.isoform.gz
isopedia extr -b ./chr22.ont.grch38.bam -o ./hg002_ont_chr22.isoform.gz

# make a manifest.tsv(tab-seprated) for *.isoform.gz files. example can be found at ./manifest.tsv

# aggregate, only first two column will be read in this step.
isopedia aggr -i manifest.tsv -o index/

# build index. provide the same manifest file, the rest of meta columns will be read.
isopedia-idx  -i index/ -m manifest.tsv 

# test your index by run a small annotation task.
isopedia isoform -i index/ -g gencode.v47.basic.chr22.gtf -o isoform.anno.tsv

```
</details>


# How it works


![how-it-works](./img/how-it-works.png)


# Usage

## Annotate (search) isoforms/transcripts

### Purpose:

search transcripts from input gtf file and return how many samples in the index have evidence. 

### Example:

```bash
isopedia isoform -i index/ -g query.gtf -f 15 -o out.tsv.gz
```

key parameters:

`--idxdir(-i)` path to index file

`--gtf(-g)` path to gtf that to be annotate

`--min-read(-m)` minimal support read in each sample to define a postive sample

`--flank(-f)` flank base pairs when searching splice sites. large value will slow down the run time but allow more wobble splice site.

<details>
<summary>
All parameters:
</summary>

```bash
Usage: isopedia isoform [OPTIONS] --idxdir <IDXDIR> --gtf <GTF>

Options:
  -i, --idxdir <IDXDIR>
          index directory

  -g, --gtf <GTF>
          gtf file

  -f, --flank <FLANK>
          flank size for search, before and after the position
          
          [default: 10]

  -m, --min-read <MIN_READ>
          minimal reads to define a positive sample
          
          [default: 1]

  -o, --output <OUTPUT>
          output file for search results

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

```

</details>


### Output


The output of the search command is a tab-separated file with the following columns:

| Column name               | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| chrom                      | Chromosome                                                                 |
| start                      | Start position of the query transcript                                      |
| end                        | End position of the query transcript                                        |
| length                     | Length of the query transcript                                              |
| exon_count                 | Number of exons in the query transcript                                     |
| trans_id                   | Transcript ID                                                              |
| gene_id                    | Gene ID                                                                    |
| confidence                  | Confidence value for detecting the query transcript in the index            |
| detected     | Whether at least one sample supports this transcript with ≥ `--min-read` reads |
| min_read                   | Minimum number of reads to define a positive sample                        |
| positive_count/sample_size  | Positive count / sample size                                                |
| attributes                 | Original attributes of the transcript from the input GTF file               |
| FORMAT                     | Format of the values in each sample column                                  |
| sample1                    | Values                                                                     |
| …                          | …                                                                           |
| sampleN                    | Values                                                                     |


There are a few columns can be used to filter the results.

`detected` this binary value indicates if at least one sample has evidence to support the query transcirpt. it can be used to quickly filterout transcirpts without evidence.

`positive_count/sample_size` this value is a combination of two values. it indicates how many samples have engouth evidence(defined by `--min-read`). it can be used to quckly filter the transcirpts that have at least several samples in the index.

`confidence` a value that summarize the confidence of observing a transcript in the entire index

<details>

$$C = \frac{k}{n}* (\prod_{i}^{n}CPM_{i})^{1/n} *G$$

where $n$ is the total number of samples in the index. $k$ is the sample number that found evidence(at least 1 support read) for a query. $CPM_{i}$ is the count per million value of the transcript in the sample $i$, which is defined as:
$$CPM_{i}=\frac{ \text{Number of support reads for the query transcript}} {\text{Total number of valid reads in the sample }i} * 1,000,000$$
$G$ is the GINI coefficient in positive samples$(i=0..k)$:

$$G = 2 \frac{\sum_{i=1}^{n} i*CPM_{i}}{n \sum_{i=1}^{n} CPM_{i} } - \frac{n+1}{n}$$


</details>


## Annotate (search) fusion genes

### Purpose:

search fusion from the index and report evidence.

### Example:

```bash
# query a single fusion
isopedia fusion -i index/ -f 10 -p chr1:pos1,chr2:pos2 -o fusion.anno.bed.gz

# query multiple fusions at the same time
isopedia fusion -i index/ -f 10 -P fusion_breakpoints.bed -o fusion_all.anno.bed.gz
```

<details>

<summary>
All parameters:
</summary>

```bash

Usage: isopedia fusion [OPTIONS] --idxdir <IDXDIR> --output <OUTPUT>

Options:
  -i, --idxdir <IDXDIR>
          index directory

  -p, --pos <POS>
          two breakpoints for gene fusion to be search(-p chr1:pos1,chr2:pos2)

  -P, --pos-bed <POS_BED>
          bed file that has the breakpoints for gene fusions. First four columns are chr1, pos1, chr2, pos2, and starts from the fifth column is the fusion id

  -G, --gene-gtf <GENE_GTF>
          bed file that has the start-end positions of the genes, used to find any possible gene fusions within the provided gene regions

  -f, --flank <FLANK>
          flank size for search, before and after the position
          
          [default: 10]

  -m, --min-read <MIN_READ>
          minimal reads to define a positive sample
          
          [default: 1]

  -o, --output <OUTPUT>
          output file for search results

      --debug
          debug mode

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```
</details>


### Output

| Column name              | Description                                |
|---------------------------|--------------------------------------------|
| chr1                      | Chromosome of the first breakpoint         |
| pos1                      | Position of the first breakpoint           |
| chr2                      | Chromosome of the second breakpoint        |
| pos2                      | Position of the second breakpoint          |
| id                        | Event or transcript ID                     |
| min_read                  | Minimum number of reads required           |
| sample_size               | Total number of samples considered         |
| positive_sample_count     | Number of positive samples                 |
| sample1                      | Value/status in sample1                |
|...|...|
| sampleN                       | Value/status in sampleN                 |



## Find potential fusion genes

### Purpose:

search fusion from the index and report evidence.

### Example:

```bash

isopedia fusion -i index/  -G gene.gtf -o fusion.discovery.out.gz
```

key parameters:

`--gene-gtf(-G)` a gtf file that has gene records. the rest of feature will be ignored.

### Output


| Column name               | Description                                             |
|----------------------------|---------------------------------------------------------|
| gene1_name                 | Name of gene 1                                         |
| gene1_id                   | ID of gene 1                                           |
| gene2_name                 | Name of gene 2                                         |
| gene2_id                   | ID of gene 2                                           |
| chr1                       | Chromosome of gene 1                                   |
| start1                     | Consensus start position for gene 1 mapped region      |
| end1                       | Consensus end position for gene 1 mapped region        |
| chr2                       | Chromosome of gene 2                                   |
| start2                     | Consensus start position for gene 2 mapped region      |
| end2                       | Consensus end position for gene 2 mapped region        |
| total_evidences            | Total number of supporting evidences                   |
| total_samples              | Total number of samples supporting the event           |
| splice_junctions_count1    | Number of splice junctions supporting gene 1           |
| splice_junctions_count2    | Number of splice junctions supporting gene 2           |
| Sample1                       | Number of supporting reads in sample1                           |
|...|...|
| SampleN                    | Number of supporting reads in sampleN                            |




## Annotate splice junctions and visualize isoforms

### Purpose:
This command is designed for cases where you have a specific splice junction of interest and want to explore its isoform context in detail. It provides both tabular output and visualization.

### Example:

```bash
isopedia splice \
  -i index/ \
  -g gencode.v47.basic.chr22.gtf \
  -p chr22:41100500-41101500 \
  -o splice.out.gz
```

key parameters:

| Option | Argument        | Description                                                                 | Default   |
|--------|-----------------|-----------------------------------------------------------------------------|-----------|
| `-i, --idxdir`       | `<IDXDIR>`      | Path to the index directory                                           | —         |
| `-s, --splice`       | `<SPLICE>`      | Splice junction in `chr1:pos1,chr2:pos2` format                      | —         |
| `-S, --splice-bed`   | `<SPLICE_BED>`  | Path to splice junction BED file                                     | —         |
| `-f, --flank`        | `<FLANK>`       | Flanking size (in bases) before and after the position               | `10`      |
| `-m, --min-read`     | `<MIN_READ>`    | Minimum number of reads required to define a positive sample         | `1`       |
| `-o, --output`       | `<OUTPUT>`      | Output file for search results (gzip-compressed)                     | —         |
| `-w, --warmup-mem`   | `<WARMUP_MEM>`  | Memory size for warming up (GB). Larger values improve performance   | `4`       |
| `-c, --cached_nodes` | `<LRU_SIZE>`    | Maximum number of cached nodes per tree                              | `100000`  |




### Output

The output is a gzip-compressed file containing detailed information about the splice junction and associated isoforms. Each isoform record includes:
| Field Name          | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| #id                 | Isoform identifier                                                          |
| chr1, pos1          | Chromosome and position of splice donor                                     |
| chr2, pos2          | Chromosome and position of splice acceptor                                  |
| total_evidence      | Total number of supporting reads                                            |
| cpm                 | Normalized counts (CPM)                                                     |
| matched_sj_idx      | Index of the matched splice junction                                        |
| dist_to_matched_sj  | Distance to the matched splice junction                                     |
| n_exons             | Number of exons in the isoform                                              |
| start_pos_left      | Leftmost starting position of isoform                                       |
| start_pos_right     | Rightmost starting position of isoform                                      |
| end_pos_left        | Leftmost ending position of isoform                                         |
| end_pos_right       | Rightmost ending position of isoform                                        |
| splice_junctions    | List of splice junctions in the isoform                                     |
| format              | Format of record                                                            |
| ENCSR***            | Per-sample evidence (columns for each dataset, e.g., ENCODE accessions)     |



Visualization output example:
https://zhengxinchang.github.io/isopedia/ 


# Memory Usage


## ENCODE long-read RNA-seq datasets(107 samples)
| Step                       | Peak Memory Usage (GB) |
|----------------------------|------------------------|
| isopedia aggr             | 7.12                   |
| isopedia-idx               | 3.84                   |
| isopedia isoform(158K transcripts from GENCODE)      | 15.82                  |



# Installation

## Check out the latest release


https://github.com/zhengxinchang/isopedia/releases


## From source code


Rust, cargo, and musl are required for building the project from source.

```bash
git clone https://github.com/zhengxinchang/isopedia.git
cd isopedia
cargo build --release
cargo build --release --target x86_64-unknown-linux-musl
```



# Quick Q&A:

**Q1: What gap does Isopedia aim to fill?**
A: Long-read RNA sequencing often reveals a large number of novel isoforms. However, evaluating whether these novel isoforms are biologically meaningful—or simply artifacts caused by RNA degradation or sequencing errors—is challenging. Understanding both the existence and population frequency of a novel isoform is essential for downstream analyses, yet current approaches are limited in scalability and robustness.

**Q2: How does Isopedia address this problem?**
A: Rather than relying on model-based classification or paired sequencing datasets from the same sample, Isopedia takes a different approach. It searches for supporting evidence of novel isoforms across large-scale long-read transcriptome datasets, leveraging population-level data to distinguish true biological isoforms from noise.

**Q3: How does Isopedia manage large datasets and provide a practical solution for isoform assessment?**
A: Isopedia introduces several innovations for scalable and efficient isoform evaluation. It uses a B+ tree–based data structure to rapidly index and compare isoform-related genomic positions. In addition, it employs a read-level signal extraction algorithm to build a compact yet informative index of transcriptome data. When assessing a query isoform, Isopedia integrates evidence from splicing junctions and alignment quality to ensure a robust evaluation. The entire tool is implemented in Rust, offering high performance and a user-friendly interface.

**Q4: What is the best use case for Isopedia?**
A: Isopedia is designed to be integrated into standard long-read transcriptome analysis pipelines. After isoforms are identified using tools like IsoQuant, FLAMES, TALON, or others, their output GTF files can be passed to Isopedia. It then evaluates the isoforms against a background index constructed from hundreds of publicly available long-read transcriptome datasets, providing a scalable and population-aware assessment.



# Contact

* zhengxc93@gmail.com
* fritz.sedlazeck@bcm.edu