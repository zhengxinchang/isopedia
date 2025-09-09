import pandas as pd
import sys

def dedup_meta_tbl(meta_tbl_path, output_path):
    # Read the metadata table
    df = pd.read_csv(meta_tbl_path, sep='\t')

    # Deduplicate based on 'unique_id', keeping the first occurrence
    dedup_df = df.drop_duplicates(subset='unique_id', keep='first')

    # Write the deduplicated table to a new file
    dedup_df.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dedup_meta_tbl.py <input_meta_tbl> <output_meta_tbl>")
        sys.exit(1)

    input_meta_tbl = sys.argv[1]
    output_meta_tbl = sys.argv[2]

    dedup_meta_tbl(input_meta_tbl, output_meta_tbl)