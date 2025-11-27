#!/usr/bin/env python3

"""
A script to profile GTF files for transcript statistics.
1. how many of different genes share the same splice site
2. how may of different genes overlaped in the genome
"""

import argparse
from collections import defaultdict
import gzip
import os
import sys

def stats_overlap(gtf_file):
    tx_vec = []
    total_gene_n = 0
    with (gzip.open(gtf_file, 'rt') if gtf_file.endswith('.gz') else open(gtf_file, 'r')) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')   
            if fields[2] == 'gene':
                total_gene_n += 1
                strand = fields[6]
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                info = fields[8]
                info_fields = {item.split(' ')[0]: item.split(' ')[1].strip('"') for item in info.split('; ') if item}
                gene_id = info_fields.get('gene_id', 'NA')
                gene_name = info_fields.get('gene_name', 'NA')
                transcript_id = info_fields.get('transcript_id', 'NA')
                
                tx_vec.append((chrom, start, end, strand, gene_id, transcript_id,gene_name))
    # sort by chrom and start
    tx_vec.sort(key=lambda x: (x[0], x[1]))

    # now profile overlaps
    overlap_dict = defaultdict(set)  # gene_id -> set of overlapping gene_ids
    for i in range(len(tx_vec)):
        chrom_i, start_i, end_i, strand_i, gene_id_i, transcript_id_i, gene_name_i = tx_vec[i]
        for j in range(i + 1, len(tx_vec)):
            chrom_j, start_j, end_j, strand_j, gene_id_j, transcript_id_j, gene_name_j = tx_vec[j]
            if chrom_i != chrom_j:
                break  # different chromosome
            if start_j > end_i:
                break  # no overlap
            # check overlap
            if not (end_i < start_j or end_j < start_i):
                overlap_dict[f"{gene_name_i}({gene_id_i},{strand_i})"].add(f"{gene_name_j}({gene_id_j},{strand_j})")
                overlap_dict[f"{gene_name_j}({gene_id_j},{strand_j})"].add(f"{gene_name_i}({gene_id_i},{strand_i})")
    # print stats
    for gene_id, overlaps in overlap_dict.items():
        print(f"Gene {gene_id} overlaps with {len(overlaps)} other genes: {', '.join(overlaps)}")
    print(f"Total number of genes: {total_gene_n}, number of genes with overlaps: {len(overlap_dict)}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Profile GTF file for transcript statistics.")
    parser.add_argument("gtf_file", help="Path to the GTF file (can be gzipped).")
    args = parser.parse_args()
    
    if not os.path.isfile(args.gtf_file):
        print(f"Error: File {args.gtf_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    
    stats_overlap(args.gtf_file)