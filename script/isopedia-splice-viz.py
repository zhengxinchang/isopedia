import gzip
import argparse
import os

class Record:
    pass

class MetaTable:
    pass

def main():
    args = parse_args()

    with gzip.open(args.input, 'rt') as infile:
        pass 


    pass


def parse_args():
    parser = argparse.ArgumentParser(description="Visualize isoforms from isopedia-anno-splice output")
    parser.add_argument("--input", required=True, help="Input file, the file is the output from isopedia-anno-splice with single query mode.")
    parser.add_argument("--output", required=True, help="Output file")
    return parser.parse_args()

