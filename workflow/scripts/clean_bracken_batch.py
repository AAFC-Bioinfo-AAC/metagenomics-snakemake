#!/usr/bin/env python3
import pandas as pd

# clean_bracken_batch.py (for Snakemake script:)
# Remove the host taxa from Bracken output files

# Define taxa to filter at each level
TAXA_FILTERS = {
    "phylum": ["Chordata"],
    "genus": ["Bos", "Sus", "Homo"],
    "species": ["Bos taurus", "Bos indicus", "Sus scrofa", "Homo sapiens"]
}

def clean_bracken_data(df, host_taxa):
    df_filtered = df[~df['name'].isin(host_taxa)].copy()
    frac_cols = ['fraction_total_reads']  # Only normalize this column
    frac_sums = df_filtered[frac_cols].sum()
    df_filtered[frac_cols] = df_filtered[frac_cols].div(frac_sums)
    return df_filtered

# You get snakemake.input, snakemake.output as dicts
for level in ("species", "genus", "phylum"):
    in_file = snakemake.input[level]
    out_file = snakemake.output[level]
    host_taxa = TAXA_FILTERS[level]
    print(f"Processing {level}: {in_file} -> {out_file}")
    if not in_file or in_file == '/':
        raise ValueError(f"Invalid input file: {in_file}")
    if not out_file or out_file == '/':
        raise ValueError(f"Invalid output file: {out_file}")
    df = pd.read_csv(in_file, sep='\t')
    cleaned = clean_bracken_data(df, host_taxa)
    cleaned.to_csv(out_file, sep='\t', index=False)
    print(f"Cleaned: {in_file} -> {out_file}")