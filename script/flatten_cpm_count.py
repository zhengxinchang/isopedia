#!/usr/bin/env python3
"""
Author: Zev Kronenberg(github name: zeeev)
Flatten a CPM:Count matrix with per-sample trailing columns into tall/long format.

Input format assumptions:
- Leading lines starting with '##' are metadata and should be skipped.
- The first non-'##' line is a tab-delimited header.
- Columns include:
    - 'trans_id' (Transcript)
    - 'attributes' (contains 'gene_name:...' among semicolon-separated key:value pairs)
    - trailing sample columns whose cells look like '<CPM>:<Count>'
- Often a 'FORMAT' column appears right before the sample columns.

Output:
A TSV with columns:
    SAMPLE    Transcript    GeneName    CPM:Count

Usage:
    python flatten_cpm_count.py INPUT.tsv [-o OUTPUT.tsv]
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import List, Optional, Tuple, TextIO
import gzip

def read_header_and_first_data(fp: TextIO) -> Tuple[str, Optional[str]]:
    """
    Return (header_line, first_data_line_or_None), skipping lines starting with '##'.
    Leaves fp positioned right after the first data line.
    """
    header_line = None
    # Find header
    for line in fp:
        if line.startswith("##"):
            continue
        header_line = line.rstrip("\n")
        break
    if header_line is None:
        raise ValueError("No main header line found (after skipping '##' lines).")
    print(f"Header: {header_line}")
    # Find first data line
    first_data = None
    for line in fp:
        if line.startswith("##"):
            continue
        first_data = line.rstrip("\n")
        break
    return header_line, first_data

def find_sample_start_idx(header_cols: List[str], first_data_cols: Optional[List[str]]) -> int:
    """
    Determine where sample columns begin.

    Strategy:
    1) If a 'FORMAT' column exists, samples begin immediately after it.
    2) Else, if we have a first data row, find the first column whose cell looks like '<something>:<something>'.
       (Sample cells are 'CPM:Count'; header cells before that are normal fields.)
    3) Fallback to the last known metadata column index + 1 (conservative).
    """
    # 1) Prefer FORMAT sentinel if present
    if "FORMAT" in header_cols:
        idx = header_cols.index("FORMAT") + 1
        if idx < len(header_cols):
            return idx

    # 2) Heuristic from first data row
    if first_data_cols is not None:
        for i, tok in enumerate(first_data_cols):
            if ":" in tok:
                # Good enough: sample cells are stored as 'CPM:Count' strings
                return i

    # 3) Fallback: try after attributes if present
    for sentinel in ("attributes", "positive_count/sample_size", "min_read"):
        if sentinel in header_cols:
            i = header_cols.index(sentinel) + 1
            if i < len(header_cols):
                return i

    # Last resort
    return max(0, len(header_cols) - 1)

def extract_gene_name(attributes: str) -> str:
    """
    attributes like: 'transcript_id:ENST...;gene_id:ENSG...;gene_name:OR4F5'
    """
    for field in attributes.split(";"):
        field = field.strip()
        if field.startswith("gene_name:"):
            return field.split(":", 1)[1]
    return ""

def flatten(
    in_path: Path,
    out_path: Optional[Path] = None,
) -> None:
    out_fp: TextIO
    if out_path:
        out_fp = out_path.open("w", newline="")
        close_out = True
    else:
        import sys
        out_fp = sys.stdout
        close_out = False

    if in_path.suffix == ".gz":
        in_path_open = gzip.open
    else:
        in_path_open = open


    with in_path_open(in_path, "rt", errors="replace") as f:
        header_line, first_data_line = read_header_and_first_data(f)
        header_cols = header_line.lstrip("#").split("\t")
        first_data_cols = first_data_line.split("\t") if first_data_line is not None else None

        # Locate key columns
        try:
            trans_idx = header_cols.index("trans_id")
        except ValueError:
            raise ValueError("Header must contain a 'trans_id' column.")

        attr_idx = header_cols.index("attributes") if "attributes" in header_cols else None

        sample_start_idx = find_sample_start_idx(header_cols, first_data_cols)
        sample_names = header_cols[sample_start_idx:]

        # Write output header
        out_fp.write("SAMPLE\tTranscript\tGeneName\tCPM:Count\n")

        def emit_row(parts: List[str]) -> None:
            if len(parts) < sample_start_idx:
                return
            transcript = parts[trans_idx] if trans_idx < len(parts) else ""
            attributes = parts[attr_idx] if (attr_idx is not None and attr_idx < len(parts)) else ""
            gene_name = extract_gene_name(attributes) if attributes else ""
            # Emit one line per sample, preserving original CPM:Count token
            for s_col, s_name in enumerate(sample_names, start=sample_start_idx):
                cell = parts[s_col] if s_col < len(parts) else "0:0"
                out_fp.write(f"{s_name}\t{transcript}\t{gene_name}\t{cell}\n")

        # Emit first data row (already read)
        if first_data_cols is not None:
            emit_row(first_data_cols)

        # Stream the rest
        for line in f:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            emit_row(parts)

    if close_out:
        out_fp.close()

def main():
    ap = argparse.ArgumentParser(description="Flatten CPM:Count matrix to SAMPLE, Transcript, GeneName, CPM:Count.")
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
