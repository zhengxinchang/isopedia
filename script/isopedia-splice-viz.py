import gzip
import argparse
import sys
from collections import OrderedDict
import json

class IsoformRecord:
    
    def __init__(self,line,sample_names):
        fields = line.strip().split("\t")
        try:
            self.id = fields[0]
            self.chr1 = fields[1]
            self.pos1 = fields[2]
            self.chr2 = fields[3]
            self.pos2 = fields[4]
            self.total_evidence = int(fields[5])
            self.cpm = float(fields[6])
            self.match_sj_idx = int(fields[7])
            self.dist_to_match_1 = int(fields[8].split(",")[0])
            self.dist_to_match_2 = int(fields[8].split(",")[1])
            self.n_exon = int(fields[9])
            self.start_pos_left = int(fields[10])
            self.start_pos_right = int(fields[11])
            self.end_pos_left = int(fields[12])
            self.end_pos_right = int(fields[13])
            splice_sites = fields[14].split(",")
            self.splice_sites = []
            self.splice_junction = []
            self.all_sample_names = sample_names
            for sj_str in splice_sites:
                a,b = sj_str.split("-")
                a = int(a)
                b = int(b)
                self.splice_sites.append(a)
                self.splice_sites.append(b)
                self.splice_junction.append((a, b))
            self.positive_samples = OrderedDict()

            for sample_idx, sample in enumerate(fields[16:]):
                if sample == "0:0:NULL":
                    continue
                count,cpm,reads = sample.split(":")
                count = int(count)
                cpm = float(cpm)
                reads = reads.split(",")
                reads_vec = []
                for read in reads:
                    start, end, strand = read.split("|")
                    start = int(start)
                    end = int(end)
                    reads_vec.append({
                        "start": start,
                        "end": end,
                        "strand": strand
                    })
                    reads_vec.sort(key=lambda x: x["start"])

                self.positive_samples[self.all_sample_names[sample_idx]] = {
                    "sample": self.all_sample_names[sample_idx],
                    "count": count,
                    "cpm": cpm,
                    "reads": reads_vec
                }
        except Exception as e:
            print(f"Error processing data line: {e}", file=sys.stderr)
            exit(1)

        self.positive_sample_size = len(self.positive_samples)
    def get_json(self):
        return {
            "id": self.id,
            "chr1": self.chr1,
            "pos1": self.pos1,
            "chr2": self.chr2,
            "pos2": self.pos2,
            "total_evidence": self.total_evidence,
            "cpm": self.cpm,
            "match_sj_idx": self.match_sj_idx,
            "dist_to_match_1": self.dist_to_match_1,
            "dist_to_match_2": self.dist_to_match_2,
            "n_exon": self.n_exon,
            "start_pos_left": self.start_pos_left,
            "start_pos_right": self.start_pos_right,
            "end_pos_left": self.end_pos_left,
            "end_pos_right": self.end_pos_right,
            "splice_sites": self.splice_sites,
            "splice_junction": self.splice_junction,
            "positive_sample_size": self.positive_sample_size,
            "positive_samples": self.positive_samples
        }


class MetaTable:
    def __init__(self):
        self.attrs = []
        self.line_no = 0
        self.data = {}
        self.sample_names = []
        self.left_most = None
        self.right_most = None

    def add_record(self, line):
        if self.line_no == 0:
            # header line
            clean_line = MetaTable._trim(line)
            fields = clean_line.strip().split("\t")
            self.attrs = fields
            self.line_no += 1
        else:
            # data line
            clean_line = MetaTable._trim(line)
            fields = clean_line.strip().split("\t")
            self.data[fields[0]] = fields[1:]
            self.line_no += 1
            self.sample_names.append(fields[0])

    def get_attr_list(self):
        return self.attrs

    def get_attr(self,sample_name, attr_name):
        if sample_name in self.data:
            record = self.data[sample_name]
            if attr_name in self.attrs:
                return record[self.attrs.index(attr_name)]
        return None

    def _trim(s):
        return s.lstrip("##")
    
    def get_json(self):
        return {
            "attrs": self.attrs,
            "data": self.data,
            "all_sample_names": self.sample_names,
            "record_no": self.line_no
        }


class GTFRecord:
    def __init__(self, line):
        fields = line.strip().split("\t")
        self.feature = fields[2]
        self.chr = fields[0]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.strand = fields[6]
        self.attributes = GTF._parse_attributes(fields[8])
        self.transcript = OrderedDict()
        self.elements = []

    def add_transcript(self, start, end, attributes):
        trans_id = attributes.get('transcript_id', None)
        self.transcript[trans_id] = {
            'transcript_id': trans_id,
            'gene_id': attributes.get('gene_id', None),
            'start': start,
            'end': end,
            'attributes': attributes,
            'elements': []
        }

    def add_elements(self, start, end, attributes):
        gene_id = attributes.get('gene_id', None)
        transcript_id = attributes.get('transcript_id', None)

        self.transcript[transcript_id]['elements'].append({
            'element_id': attributes.get('element_id', None),
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'start': start,
            'end': end,
            'attributes': attributes
        })

    def get_json(self):
        # print(f"Gene {self.attributes.get('gene_id', None)} has {len(self.transcript)} transcripts", file=sys.stderr)
        return {
            "feature": self.feature,
            "gene_id": self.attributes.get("gene_id", None),
            "chr": self.chr,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "attributes": self.attributes,
            "transcripts": self.transcript,
        }


class GTF:

    def __init__(self,gtf_path,q_chr,q_start,q_end):

        self.gene_records = {}

        with open(gtf_path,'r') as infile:
            for line in infile:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                # print(len(fields) )
                if len(fields) < 9:
                    continue
                feature_type = fields[2]
                # print(feature_type)
                chrom = trim_chr(fields[0])
                start = int(fields[3])
                end = int(fields[4])

                if chrom != q_chr:
                    continue
                
                if feature_type == "gene":
                    
                    if GTF._overlap(start,end,q_start,q_end) > 0 and chrom == trim_chr(q_chr):
                        ## gene that overlap with the query

                        attrs = GTF._parse_attributes(fields[8])
                        gene_id = attrs.get('gene_id', None)
                        self.gene_records[gene_id] = GTFRecord(line)

                elif feature_type == "transcript":
                    attrs = GTF._parse_attributes(fields[8])
                    gene_id = attrs.get('gene_id', None)
                    if gene_id in self.gene_records.keys():
                        self.gene_records[gene_id].add_transcript(start, end, attrs)
                else:
                    attrs = GTF._parse_attributes(fields[8])
                    gene_id = attrs.get('gene_id', None)
                    if gene_id in self.gene_records.keys():
                        self.gene_records[gene_id].add_elements(start, end, attrs)



    def _overlap(a1,a2,b1,b2):
        return max(0, min(a2, b2) - max(a1, b1))

    def _parse_attributes(string):
        attributes = {}
        for attr in string.split(";"):
            attr = attr.strip()
            if attr == "":
                continue
            key_value = attr.split(" ")
            if len(key_value) >= 2:
                key = key_value[0]
                value = " ".join(key_value[1:]).strip('"')
                attributes[key] = value
        return attributes


    def get_json(self):
        gene_list = []
        for gene in self.gene_records.values():
            # print("aa",gene.get_json())
            gene_list.append(gene.get_json())
        return gene_list


def trim_chr(chrom):

    chr = chrom.strip("chr").upper()
    if chr == "MT":
        return "M"
    return chr


def main():
    args = parse_args()

    metatable = MetaTable()
    isoform_records = []
    sample_names = []
    current_id = None

    left_most = sys.maxsize
    right_most = 0
    with gzip.open(args.input, 'rt') as infile:
        for line in infile:
            if line.startswith("##"):
                metatable.add_record(line)
            elif line[0] == "#" and line[1] != "#":
                # header lines:
                sample_names = line.strip().split("\t")[16:]
                # print(f"Sample names: {sample_names}", file=sys.stderr)
            else:
                if len(sample_names) == 0:
                    print("Error: No header line found before data lines.", file=sys.stderr)
                    exit(1)

                isoform_record = IsoformRecord(line,sample_names)
                if current_id is None:
                    current_id = isoform_record.id
                if current_id != isoform_record.id:
                    print(f"Error: Only single query are supported, however multiple query were found in the input file.",file=sys.stderr)
                    exit(1)

                if isoform_record.start_pos_left < left_most:
                    left_most = isoform_record.start_pos_left
                if isoform_record.end_pos_right > right_most:   
                    right_most = isoform_record.end_pos_right
                isoform_records.append(isoform_record)

    chrom = isoform_record.chr1
    # print(chrom,left_most,right_most)
    gtf = GTF(gtf_path=args.gtf,q_chr=chrom,q_start=left_most,q_end=right_most)


    # print(gtf.get_json(),file=sys.stderr)

    out_json = {
        "meta": metatable.get_json(),
        'annotation': gtf.get_json(),
        "isoforms": [record.get_json() for record in isoform_records]
    }

    with open(args.output, 'w') as outfile:
        json.dump(out_json, outfile)

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize isoforms from isopedia-anno-splice output")
    parser.add_argument("-i","--input", required=True, help="Input file, the file is the output from isopedia-anno-splice with single query mode.")
    parser.add_argument("-g","--gtf", required=False, help="Reference GTF file")
    parser.add_argument("-o","--output", required=True, help="Output file")
    return parser.parse_args()

if __name__ == "__main__":
    main()