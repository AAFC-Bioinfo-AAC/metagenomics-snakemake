#!/usr/bin/env python3
"""
This script removes the host taxa from Bracken output files (as defined in the TAXA_FILTERS)
Then the samples are re-normalized using the total remaining read counts (recomputing *_frac columns from *_num columns).
"""

import pandas as pd
import sys

# Redirect stdout and stderr to the log file provided by Snakemake
if hasattr(snakemake, "log") and snakemake.log:
    sys.stdout = open(snakemake.log[0], "w")
    sys.stderr = sys.stdout

# Define taxa to filter at each level
TAXA_FILTERS = {
    "domain": ["Eukaryota"],         # keep only prokaryotes in domain tables
    "phylum": ["Chordata"],
    "genus": ["Bos", "Sus", "Homo"],
    "species": ["Bos taurus", "Bos indicus", "Sus scrofa", "Homo sapiens"]
}

def recompute_fracs_from_nums(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each sample, sets '<sample>..._frac' = '<sample>..._num' / sum('<sample>..._num').
    If the sum is zero, sets the fraction column to 0.0.
    """
    num_cols = [c for c in df.columns if c.endswith('_num')]
    for num_col in num_cols:
        base = num_col[:-4]  # strip '_num'
        frac_col = f"{base}_frac"
        if frac_col in df.columns:
            nums = pd.to_numeric(df[num_col], errors='coerce').fillna(0)
            total = nums.sum()
            df[frac_col] = (nums / total) if total > 0 else 0.0
    return df

def clean_bracken_data(df: pd.DataFrame, host_taxa: list) -> pd.DataFrame:
    # Drop host/unwanted taxa by exact match on 'name'
    df_filtered = df[~df['name'].isin(host_taxa)].copy()
    # Recompute each *_frac from its *_num AFTER filtering
    df_filtered = recompute_fracs_from_nums(df_filtered)
    return df_filtered

def main():
    # For each level defined in input
    for level in ("species", "genus", "phylum", "domain"):
        # Only proceed if the rule actually specified this level
        if not hasattr(snakemake.input, level):
            continue
        in_file = snakemake.input[level]
        out_file = snakemake.output[level]
        host_taxa = TAXA_FILTERS[level]

        print(f"Processing {level}: {in_file} -> {out_file}")

        df = pd.read_csv(in_file, sep='\t')
        cleaned = clean_bracken_data(df, host_taxa)
        cleaned.to_csv(out_file, sep='\t', index=False)

        # Sanity check: show min-max of frac-column sums across samples
        frac_cols = [c for c in cleaned.columns if c.endswith('_frac')]
        if frac_cols:
            col_sums = cleaned[frac_cols].sum(axis=0)
            print("Cleaned {}: frac-column sums range: {:.6f} to {:.6f}".format(
                level, col_sums.min(), col_sums.max()))
        else:
            print(f"Cleaned {level}: (no *_frac columns found)")

if __name__ == "__main__":
    main()