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


def main():
    args = parse_args()

    metatable = MetaTable()
    isoform_records = []
    sample_names = []
    current_id = None
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
                isoform_records.append(isoform_record)


    out_json = {
        "meta": metatable.get_json(),
        "isoforms": [record.get_json() for record in isoform_records]
    }

    with open(args.output, 'w') as outfile:
        json.dump(out_json, outfile, indent=4)

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize isoforms from isopedia-anno-splice output")
    parser.add_argument("-i","--input", required=True, help="Input file, the file is the output from isopedia-anno-splice with single query mode.")
    parser.add_argument("-g","--gtf", required=False, help="Reference GTF file")
    parser.add_argument("-o","--output", required=True, help="Output file")
    return parser.parse_args()

if __name__ == "__main__":
    main()