# About the Isopedia <img src="./img/logo3.png" align="right" alt="" width=120 />

**Isopedia** is a scalable tool designed for the analysis of hundreds to thousands of long-read transcriptome datasets simultaneously by read-level indexing approach. It provides two key capabilities:

**Population-level transcript quantification and frequency profiling** — With a single command, Isopedia can quantify isoform expression and estimate their occurrence frequencies across large cohorts using minimal computational resources.

**Isoform diversity exploration and visualization** — Isopedia enables systematic analysis of fusion genes and specific splicing events across populations, offering insights into transcript diversity beyond individual samples.


# Table of Content

- [Quick Start](#quick-start)
- [How it works](#how-it-works)
- [Download pre-built index](#download-pre-built-index)
- [Build your own index](#build-your-own-index)
- [How to use Isopedia](#how-to-use-isopedia)
  - [Query(quantification) transcripts](#annotate-search-isoformstranscripts)
  - [Query(search) fusion genes](#annotate-search-fusion-genes)
  - [Find potential fusion genes](#find-potential-fusion-genes)
  - [Query(search)/visualization of splice junction](#query-splice-junctions-and-visualize-isoforms)
- [How to install Isopedia](#how-to-install-isopedia)
- [Annotate ORF with ORFannotate](#annotate-orf-with-orfannotate)
<!-- - [Quick Q&A](#quick-qa) -->


# Quick Start

Isopedia has two binaries: `isopedia` and `isopedia-tools`. The main binary `isopedia` is used for all the main functions, and `isopedia-tools` has some helper functions.

**Download prebuild index and run**

```bash
# install isopedia from conda
conda install -y -c zhengxinchang isopedia

# clone the repo and enter the toy_ex directory
git clone https://github.com/zhengxinchang/isopedia && cd isopedia/toy_ex/

# query transcripts
isopedia isoform -i index/ -g query.gtf -o ./out.profile.tsv.gz

# query one fusion gene (two breakpoints)
isopedia fusion -i index/ -p chr1:181130,chr1:201853853 -o ./out.fusion.tsv.gz

# query multiple fusion genes
isopedia fusion -i index/ -P fusion_query.bed -o ./out.fusion.tsv.gz

# query gene regions and discover potential fusion events
isopedia fusion -i index/ -G gene.gtf -o ./out.fusion.discovery.tsv.gz

# query a splice junction from the BID gene
isopedia splice -i index/ -s 22:17744013,22:17750104 -o ./out.splice.tsv.gz

# visualize the splice junction
python script/isopedia-splice-viz.py -i ./out.splice.tsv.gz -g gencode.v47.basic.chr22.gtf -o isopedia-splice-view
```

For indexing GTF files, please refer to [Indexing GTF Files](doc/indexing_gtf.md) section.


# How it works


![how-it-works](./img/how-it-works.png)

The workflow of Isopedia involves several key steps, including isoform profiling, merging, indexing, and querying. Users can start by profiling isoform signals from individual BAM files, then merge the results to build a comprehensive index. Once the index is ready, it can be used to query isoforms, fusion genes, and explore splice junctions across multiple samples. Users can also visualize specific splicing events using the provided visualization tools. 


<details>
<summary> 
<strong>How Isopedia determines a positive hit for a query transcript?</strong>
</summary>

<img src="./img/how-it-works2.png" align="center" alt=""  />

</details>



# Download pre-built index

Isopedia comes with pre-built indexes from hundreds of publicly available long-read RNA-seq datasets, which can be used directly for isoform and fusion gene annotation.

|Name|Organism|Reference|Link|Sample Size|Index size(compressed size)|Minimal required Memory for querying | Description |
|----|----|-----------|-----|------| -----|-----|-----|
|isopedia_index_hs_v1.0| Homo sapiens| GRCh38| Download ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz MD5SUM: ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz.md5|1,007|305G(107G)|64Gb|107 ENCODE samples + 900 SRA samples|

**Download and unpack the index:**

```bash
wget ftp://hgsc-sftp1.hgsc.bcm.tmc.edu//rt38520/isopedia_index_hs_v1.0.tar.xz
tar -xvf isopedia_index_hs_v1.0.tar.xz
```

**Unpack the index**

The index file was compressed by xz to achieve the highest compression rate, it can be decompressed by:

```bash
tar -xvf isopedia_index_hs_v1.0.tar.xz
```

# Build your own index

Isopedia supports building local index in your own datasets. Prerequisites are listed below:

1. A set of mapped bam files (sorted bam are not required)
2. A manifest file that describes the sample name, isoform file path, and other optional metadata in tabular (\t separated and with a header line) format. 

Example manifest file:

```
sample_name   path   platform
HG002_pb_chr22   /path/to/hg002_pb_chr22.isoform.gz   PacBio
HG002_ont_chr22   /path/to/hg002_ont_chr22.isoform.gz   ONT
```

## Example workflow

```bash
# make sure isopedia is in your $PATH or use absolute path to the binaries.

# download the toy_ex 
git clone https://github.com/zhengxinchang/isopedia && cd isopedia/toy_ex/

# profile isoform signals on each bam individually
isopedia profile -i ./chr22.pb.grch38.bam -o ./hg002_pb_chr22.isoform.gz
isopedia profile -i ./chr22.ont.grch38.bam -o ./hg002_ont_chr22.isoform.gz

# Optional flags for profile command:
# --tid : include transcript id in the profile file, only works when -g (GTF) is passed.
# --gid : include gene id in the profile file, only works when -g (GTF) is passed.
# --rname : include read name in the profile file, only works when -i (BAM/CRAM) is passed.

# make a manifest.tsv (tab-separated) for *.isoform.gz files. example can be found at ./manifest.tsv

# merging, only first two columns will be read in this step.
isopedia merge -i manifest.tsv -o index/

# build index. provide the same manifest file, the rest of meta columns will be read.
isopedia index -i index/ -m manifest.tsv 

# test your index by running a small annotation task.
isopedia isoform -i index/ -g gencode.v47.basic.chr22.gtf -o out.isoform.tsv.gz
```

# How to use Isopedia

## Query(quantification) transcripts

### Purpose:

Search transcripts from input gtf file and return how many samples in the index have evidence. 

### Example:

```bash
# Isopedia assumes the GTF file is sorted.
gffread -T -o- origin.gtf | sort -k1,1 -k4,4n | gffread - -o query.gtf

isopedia isoform -i index/ -g query.gtf -o out.tsv.gz
```

### Key parameters:

`--min-read(-m)` Minimum support reads in each sample to define a positive sample (default: 1)

`--flank(-f)` Flank base pairs when searching splice sites. Larger value will slow down the run time but allow more wobble splice sites (default: 10)

`--info` Include additional information such as read name (indexing BAM/CRAM), transcript IDs (indexing GTF), and/or gene IDs (indexing GTF) in the annotation output

`--num-threads(-n)` Number of threads to use (default: 4)

### EM Algorithm Parameters:

`--em-max-iter` Maximum EM iterations (default: 100)

`--em-conv-min-diff` EM convergence threshold (default: 0.01)

`--em-chunk-size` EM chunk size, reduce it if you have low memory (default: 4)

`--em-effective-len-coef` EM effective length coefficient, avoid divide by zero when transcript is very short (default: 2)

`--em-damping-factor` EM damping factor (default: 0.3)

`--min-em-abundance` Minimum EM abundance to report (default: 0.0001)

### TSS/TES Filtering Parameters:

`--no-check-tss-tes` Disable TSS and TES checking

`--tss-degrad-bp` Maximum allowed degradation bp for TSS (default: 2000)

`--tes-degrad-bp` Maximum allowed degradation bp for TES (default: 8000)

`--terminal-tolerance-bp` Maximum allowed deviation (bp) beyond annotated TSS and TES. Isoforms whose TSS and TES fall outside the annotation but within this tolerance are still classified as FSM (default: 10)

### Performance Tuning Parameters:

`--cached-nodes(-c)` Maximum number of cached tree nodes in memory (default: 10)

`--cached-chunk-num` Maximum number of cached isoform chunks in memory (default: 4)

`--cached-chunk-size-mb` Cached isoform chunk size in MB (default: 128)

`--output-tmp-shard-counts` Output temporary shard record size (default: 10000)

`--verbose` Enable verbose mode for detailed logging

<details>
<summary>
All parameters:
</summary>

```bash
Usage: isopedia isoform [OPTIONS] --idxdir <IDXDIR> --gtf <GTF> --output <OUTPUT>

Options:
  -i, --idxdir <IDXDIR>
          Path to the index directory

  -g, --gtf <GTF>
          Path to the GTF file

  -f, --flank <FLANK>
          Flanking size (in bases) before and after the position
          [default: 10]

  -m, --min-read <MIN_READ>
          Minimum number of reads required to define a positive sample
          [default: 1]

  -o, --output <OUTPUT>
          Output file for search results

      --info
          Whether to include additional information in the output

  -n, --num-threads <NUM_THREADS>
          Number of threads to use
          [default: 4]

      --em-max-iter <EM_MAX_ITER>
          Max EM iterations
          [default: 100]

      --em-conv-min-diff <EM_CONV_MIN_DIFF>
          EM convergence threshold
          [default: 0.01]

      --em-chunk-size <EM_CHUNK_SIZE>
          EM chunk size, reduce it if you have low memory
          [default: 4]

      --em-effective-len-coef <EM_EFFECTIVE_LEN_COEF>
          EM effective length coefficient, avoid divide by zero when transcript is very short
          [default: 2]

      --em-damping-factor <EM_DAMPING_FACTOR>
          EM damping factor
          [default: 0.3]

      --min-em-abundance <MIN_EM_ABUNDANCE>
          Minimum EM abundance to report
          [default: 0.0001]

      --no-check-tss-tes
          No check TSS and TES

      --tss-degrad-bp <TSS_DEGRAD_BP>
          Maximum allowed degradation bp for TSS
          [default: 2000]

      --tes-degrad-bp <TES_DEGRAD_BP>
          Maximum allowed degradation bp for TES
          [default: 8000]

      --terminal-tolerance-bp <TERMINAL_TOLERANCE_BP>
          Maximum allowed deviation (bp) beyond annotated TSS and TES
          [default: 10]

  -c, --cached-nodes <CACHED_NODES>
          Maximum number of cached tree nodes in memory
          [default: 10]

      --cached-chunk-num <CACHED_CHUNK_NUM>
          Maximum number of cached isoform chunks in memory
          [default: 4]

      --cached-chunk-size-mb <CACHED_CHUNK_SIZE_MB>
          Cached isoform chunk size in Mb
          [default: 128]

      --verbose
          Verbose mode

      --output-tmp-shard-counts <OUTPUT_TMP_SHARD_COUNTS>
          Output temporary shard record size
          [default: 10000]

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
| sample1                    | Values (CPM:COUNT:INFO)                                                                   |
| …                          | …                                                                           |
| sampleN                    | Values (CPM:COUNT:INFO)                                                                      |


There are a few columns that can be used to filter the results.

`detected` This binary value indicates if at least one sample has evidence to support the query transcript. It can be used to quickly filter out transcripts without evidence.

`positive_count/sample_size` This value is a combination of two values. It indicates how many samples have enough evidence (defined by `--min-read`). It can be used to quickly filter the transcripts that have at least several samples in the index.

`confidence` A value that summarizes the confidence of observing a transcript in the entire index

*CPM* values are provided in each sample column, which is defined as:

$$CPM=\frac{ \text{Number of support reads for the query transcript}} {\text{Total number of valid reads in the sample}} \times 1,000,000 $$ 

`INFO` field in each sample column contains additional information depending on the index content, such as read names (if BAM/CRAM files were indexed) or transcript/gene IDs (if GTF files were indexed). This information can be useful for further analysis or validation. It is recommended to include this field when building the index from GTF files.

**Note that** 

The `INFO` field will only be included when the index was created with `--rname` (BAM/CRAM) or `--tid`/`--gid` option and the `--info` flag is used during the query.


<details>
<summary>Confidence score calculation details</summary>

$$C = \frac{k}{n} \times (\prod_{i}^{n}CPM_{i})^{1/n} \times G$$

where $n$ is the total number of samples in the index. $k$ is the sample number that found evidence (at least 1 support read) for a query. $CPM_{i}$ is the count per million value of the transcript in the sample $i$, which is defined as:

$$CPM_{i}=\frac{ \text{Number of support reads for the query transcript}} {\text{Total number of valid reads in the sample }i} \times 1,000,000$$

$G$ is the GINI coefficient in positive samples $(i=0..k)$:

$$G = 2 \frac{\sum_{i=1}^{n} i \times CPM_{i}}{n \sum_{i=1}^{n} CPM_{i} } - \frac{n+1}{n}$$

</details>


## Query fusion gene breakpoints

### Purpose:

This command is used to search for evidence of specific gene fusion events in the index based on provided breakpoints.

### Example:

```bash
# query a single fusion
isopedia fusion -i index/ -f 10 -p chr1:pos1,chr2:pos2 -o fusion.anno.bed.gz

# query multiple fusions at the same time
isopedia fusion -i index/ -f 10 -P fusion_breakpoints.bed -o fusion_all.anno.bed.gz

```

### Key parameters:

`--min-read(-m)` Minimum support reads in each sample to define a positive sample (default: 1)

`--flank(-f)` Flank base pairs when searching splice sites. Larger value will slow down the run time but allow more wobble splice sites (default: 10)

`--pos(-p)` Two breakpoints for gene fusion to be searched (format: chr1:pos1,chr2:pos2)

`--pos-bed(-P)` BED file with breakpoints for gene fusions. First six columns are chr1, start, end, chr2, start, end, and from the seventh column is the fusion id

`--region(-r)` Query two regions for possible gene fusions, returns any reads that have breakpoints within the two regions (format: chr1:start-end,chr2:start-end)

`--gene-gtf(-G)` GTF file with gene annotations for discovering potential gene fusions

### Performance Tuning Parameters:

`--cached_nodes(-c)` Number of cached nodes for each tree in maximal (default: 1000000)

`--cached-chunk-number` Maximum number of cached isoform chunks in memory (default: 4)

`--cached-chunk-size-mb` Cached isoform chunk size in MB (default: 128)

`--verbose` Enable verbose mode for detailed logging

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
          bed file that has the breakpoints for gene fusions. First six columns are 
          chr1, start, end, chr2, start, end, and starts from the seventh column is 
          the fusion id

  -r, --region <REGION>
          query two regions for possible gene fusions, return any reads that has 
          breakpoints within the two regions. Format: chr1:start-end,chr2:start-end

  -G, --gene-gtf <GENE_GTF>
          bed file that has the start-end positions of the genes, used to find any 
          possible gene fusions within the provided gene regions

  -f, --flank <FLANK>
          flank size for search, before and after the position
          [default: 10]

  -m, --min-read <MIN_READ>
          minimal reads to define a positive sample
          [default: 1]

  -o, --output <OUTPUT>
          output file for search results

  -c, --cached_nodes <LRU_SIZE>
          number of cached nodes for each tree in maximal
          [default: 1000000]

      --cached-chunk-number <CACHED_CHUNK_NUMBER>
          Maximum number of cached isoform chunks in memory
          [default: 4]

      --cached-chunk-size-mb <CACHED_CHUNK_SIZE_MB>
          Cached isoform chunk size in Mb
          [default: 128]

      --verbose
          Verbose mode

  -h, --help
          Print help (see a summary with '-h')
```
</details>

### Example BED format for fusion breakpoints:

```
chr2    50795173    50795273    chr17    61368325    61368425    BCAS4:BCAS3,UHR(optional)
```

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

Query candidate fusion genes within specified gene regions. It identifies potential fusion events by examining all possible gene pairs within the provided regions and reporting those with supporting evidence in the index.

### Example:

```bash
isopedia fusion -i index/ -G gene.gtf -o fusion.discovery.out.gz
```

### Key parameters:

`--gene-gtf(-G)` A GTF file that has gene records. The rest of the features will be ignored.

`--min-read(-m)` Minimum support reads in each sample to define a positive sample (default: 1)

`--flank(-f)` Flank base pairs when searching splice sites. Larger value will slow down the run time but allow more wobble splice sites (default: 10)

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




## Query splice junctions and visualize isoforms

### Purpose:
This command is designed for cases where you have a specific splice junction of interest and want to explore its isoform context in detail. It provides both tabular output and visualization.

### Example:

```bash
isopedia splice -i index/ -s chr22:41100500,chr22:41101500 -o splice.out.gz
python script/isopedia-splice-viz.py -i splice.out.gz -g gencode.v47.basic.annotation.gtf -o isopedia-splice-view
```

### Key parameters (isopedia splice):

`--splice(-s)` Splice junction in 'chr1:pos1,chr2:pos2' format

`--splice-bed(-S)` Path to splice junction BED file

`--min-read(-m)` Minimum support reads in each sample to define a positive sample (default: 1)

`--flank(-f)` Flank base pairs when searching splice sites. Larger value will slow down the run time but allow more wobble splice sites (default: 10)

### Performance Tuning Parameters:

`--cached_nodes(-c)` Maximum number of cached nodes per tree (default: 1000)

`--cached-chunk-number` Maximum number of cached isoform chunks in memory (default: 4)

`--cached-chunk-size-mb` Cached isoform chunk size in MB (default: 128)

`--verbose` Enable verbose mode for detailed logging

### Key parameters (isopedia-splice-viz.py):

`-g/--gtf` A GTF file that has gene annotations. It will be used to annotate the splice junction and isoforms.

`-t/--template` A template HTML file that will be used to generate the interactive visualization.


<details>
<summary>
All parameters:
</summary>

```bash
Usage: isopedia splice [OPTIONS] --idxdir <IDXDIR> --output <OUTPUT>

Options:
  -i, --idxdir <IDXDIR>
          Path to the index directory

  -s, --splice <SPLICE>
          Splice junction in 'chr1:pos1,chr2:pos2' format

  -S, --splice-bed <SPLICE_BED>
          Path to splice junction bed file

  -f, --flank <FLANK>
          Flanking size (in bases) before and after the position
          [default: 10]

  -m, --min-read <MIN_READ>
          Minimum number of reads required to define a positive sample
          [default: 1]

  -o, --output <OUTPUT>
          Output file for search results

  -c, --cached_nodes <LRU_SIZE>
          Maximum number of cached nodes per tree
          [default: 1000]

      --cached-chunk-number <CACHED_CHUNK_NUMBER>
          Maximum number of cached isoform chunks in memory
          [default: 4]

      --cached-chunk-size-mb <CACHED_CHUNK_SIZE_MB>
          Cached isoform chunk size in Mb
          [default: 128]

      --verbose
          Verbose mode

  -h, --help
          Print help (see a summary with '-h')
```

```bash
usage: isopedia-splice-viz.py [-h] -i INPUT [-g GTF] [-t TEMPLATE] -o OUTPUT

Visualize isoforms from isopedia-anno-splice output

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file, the file is the output from isopedia splice with single query mode.
  -g GTF, --gtf GTF     Reference GTF file
  -t TEMPLATE, --template TEMPLATE
                        Templates HTML file
  -o OUTPUT, --output OUTPUT
                        Output file
```
</details>

### Example format of splice_bed:

```
id    chr1    pos1    chr2    pos2
```

**Note:** If you are using coordinates from GFF/GTF, please convert the 1-based to 0-based coordinates.

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
| SAMPLE columns      | Per-sample evidence, the format is: support_read:CPM:[read_start\|read_end\|strand,]*N   |



Visualization output example:
https://zhengxinchang.github.io/isopedia/ 

![isopedia-splice-view](./img/isopedia-splice-view.png)


# Computational resource usage


Isopedia 1.4.0 was tested on 1,007 long-read transcriptome datasets from SRA and ENCODE.

> Tested on AMD Ryzen 9 7940HX, 64 GB RAM.


| Step                       | Peak Mem (GB)  | Time(H:MM:SS) |
|----------------------------|------------------------ |---|
| isopedia merge             | 54.77 Gb                 | 28:26:08 |
| isopedia index               | 45.85 Gb              | 5:48:55  |
| isopedia isoform(280K GENCODE V49 basic)      | 9.19 Gb  |   4:52:18|



# How to install Isopedia

## Install from conda

```bash
conda install -c zhengxinchang isopedia
```

## Check out the latest release

https://github.com/zhengxinchang/isopedia/releases

Note that the `isopedia-<version>.linux.tar.gz` was compiled in Amazon Linux 2 with GCC 7.3, Glibc 2.26, and Binutils 2.29.1. It should work in most Linux distributions. However, if your Linux distribution cannot run it, you can still try to use the `isopedia-<version>.musl.tar.gz` (statically linked with musl). 

## From source code


Rust, cargo, and musl are required for building the project from source.

```bash
git clone https://github.com/zhengxinchang/isopedia.git
cd isopedia
cargo build --release
cargo build --release --target x86_64-unknown-linux-musl
```


# Annotate ORF with ORFannotate

Isopedia is designed to interoperate with [ORFannotate](https://github.com/egustavsson/ORFannotate), which can process Isopedia outputs to predict ORFs and UTRs, annotate CDS features in GTF files, and produce transcript-level summaries.

Example to use ORF annotate:

```
<place holder>
```


# Development Roadmap

- Built-in download for remote indexes

# Contact

* zhengxc93@gmail.com
* fritz.sedlazeck@bcm.edu