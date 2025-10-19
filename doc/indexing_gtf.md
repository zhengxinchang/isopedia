
# Indexing GTF Files

Isopedia has capabilities to index GTF files for isoform annotation. The main difference is using the `GTF` as input in the `isopedia profile` command instead of BAM/CRAM files.

Here is an example command to indexing GTF files:

```bash
isopedia profile  \
    -g /path/to/input.gtf \  
    -o /path/to/output.gtf_isoforms.gz \
    --tid \ # include transcript IDs in the output
    --gid   # include gene IDs in the output
```

`--tid` and `--gid` are optional flags to include transcript IDs and gene IDs in the output file, respectively. They are disabled by default but can be useful for downstream analysis.

After indexing, the output file can be used in the same way as the isoform profile files generated from BAM/CRAM files for `merge`, `index`, and `isoform` subcommands.

Particularly, when using the `isoform` subcommand, `--info` flag can be used to include additional information such as transcript IDs and gene IDs in the annotation output.

```bash
isopedia isoform \
    -i /path/to/index/ \  # path to the isopedia index directory
    -g /path/to/input.gtf \  # path to the input GTF file
    -o /path/to/output.anno.isoform.gz \  # path to the output annotation file
    --info  # include additional information in the output
```

By doing this, the output annotation file will contain additional fields in each sample column for transcript IDs and gene IDs.

## Note on CPM Field 

When querying against GTF based index, the `CPM` (Counts Per Million) field in the output annotation file will always be no sense, as each transcript in GTF file is represented as one isoform with count 1. Thus, the `CPM` field will not provide meaningful information in this context.

## Note on COUNT Field

As isopedia merge/compare transcripts by `FSM`(full splice match) status, the transcripts from GTF files that share same splice junctions will be merged together. Therefore, the `COUNT` field in the output annotation file may be greater than 1, indicating that multiple transcripts from the GTF file correspond to the same isoform in the index. The transcript ID and gene ID fields will list all the corresponding IDs separated by commas.