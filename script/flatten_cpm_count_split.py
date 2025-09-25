#!/usr/bin/env python3
"""
Author: Zev Kronenberg(github account:zeeev)
Flatten a CPM:Count matrix into tall format with separate CPM and Count columns.

Output columns:
    Sample    Transcript    GeneName    CPM    Count
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import List, Optional, Tuple, TextIO

def read_header_and_first_data(fp: TextIO) -> Tuple[str, Optional[str]]:
    header_line = None
    for line in fp:
        if line.startswith("##"):
            continue
        header_line = line.rstrip("\n")
        break
    if header_line is None:
        raise ValueError("No main header line found (after skipping '##' lines).")

    first_data = None
    for line in fp:
        if line.startswith("##"):
            continue
        first_data = line.rstrip("\n")
        break
    return header_line, first_data

def find_sample_start_idx(header_cols: List[str], first_data_cols: Optional[List[str]]) -> int:
    if "FORMAT" in header_cols:
        idx = header_cols.index("FORMAT") + 1
        if idx < len(header_cols):
            return idx
    if first_data_cols is not None:
        for i, tok in enumerate(first_data_cols):
            if ":" in tok:
                return i
    for sentinel in ("attributes", "positive_count/sample_size", "min_read"):
        if sentinel in header_cols:
            i = header_cols.index(sentinel) + 1
            if i < len(header_cols):
                return i
    return max(0, len(header_cols) - 1)

def extract_gene_name(attributes: str) -> str:
    for field in attributes.split(";"):
        field = field.strip()
        if field.startswith("gene_name:"):
            return field.split(":", 1)[1]
    return ""

def split_cpm_count(cell: str) -> tuple[float, int]:
    if ":" not in cell:
        return 0.0, 0
    left, _, right = cell.partition(":")
    try:
        cpm = float(left)
    except Exception:
        cpm = 0.0
    try:
        cnt = int(float(right))
    except Exception:
        cnt = 0
    return cpm, cnt

def flatten(in_path: Path, out_path: Optional[Path] = None) -> None:
    out_fp: TextIO
    if out_path:
        out_fp = out_path.open("w", newline="")
        close_out = True
    else:
        import sys
        out_fp = sys.stdout
        close_out = False

    with in_path.open("r", errors="replace") as f:
        header_line, first_data_line = read_header_and_first_data(f)
        header_cols = header_line.lstrip("#").split("\t")
        first_data_cols = first_data_line.split("\t") if first_data_line is not None else None

        try:
            trans_idx = header_cols.index("trans_id")
        except ValueError:
            raise ValueError("Header must contain a 'trans_id' column.")

        attr_idx = header_cols.index("attributes") if "attributes" in header_cols else None

        sample_start_idx = find_sample_start_idx(header_cols, first_data_cols)
        sample_names = header_cols[sample_start_idx:]

        # Write new header
        out_fp.write("Sample\tTranscript\tGeneName\tCPM\tCount\n")

        def emit_row(parts: List[str]) -> None:
            if len(parts) < sample_start_idx:
                return
            transcript = parts[trans_idx] if trans_idx < len(parts) else ""
            attributes = parts[attr_idx] if (attr_idx is not None and attr_idx < len(parts)) else ""
            gene_name = extract_gene_name(attributes) if attributes else ""
            for s_col, s_name in enumerate(sample_names, start=sample_start_idx):
                cell = parts[s_col] if s_col < len(parts) else "0:0"
                cpm, cnt = split_cpm_count(cell)
                out_fp.write(f"{s_name}\t{transcript}\t{gene_name}\t{cpm}\t{cnt}\n")

        if first_data_cols is not None:
            emit_row(first_data_cols)
        for line in f:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            emit_row(parts)

    if close_out:
        out_fp.close()

def main():
    ap = argparse.ArgumentParser(description="Flatten to Sample, Transcript, GeneName, CPM, Count.")
    ap.add_argument("input", help="Input TSV path")
    ap.add_argument("-o", "--output", help="Output TSV path (default: stdout)")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise FileNotFoundError(in_path)

    out_path = Path(args.output) if args.output else None
    flatten(in_path, out_path)

if __name__ == "__main__":
    main()
