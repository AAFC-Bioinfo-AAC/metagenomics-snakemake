#!/usr/bin/env python3

import pandas as pd
import os
from pathlib import Path

# Input directory where Bracken files are located
INPUT_DIR = Path("/gpfs/fs7/aafc/projects/J-002259_milk_microbiome_gut/data/final/Sow_milk_metagenomes/Kraken2_with_GTDB226_w_host/Bracken")

# Taxa to remove by file type
TAXA_FILTERS = {
    "phylum": ["Chordata"],
    "genus": ["Bos", "Sus", "Homo"],
    "species": ["Bos taurus", "Bos indicus", "Sus scrofa", "Homo sapiens"]
}

def clean_bracken_data(df, host_taxa):
    df_filtered = df[~df['name'].isin(host_taxa)].copy()
    frac_cols = [col for col in df.columns if col.endswith('_frac')]
    frac_sums = df_filtered[frac_cols].sum()
    df_filtered[frac_cols] = df_filtered[frac_cols].div(frac_sums)
    return df_filtered

def process_file(filepath, host_taxa):
    filename = filepath.stem
    output_path = filepath.parent / f"{filename}_cleaned.txt"
    df = pd.read_csv(filepath, sep='\t')
    cleaned_df = clean_bracken_data(df, host_taxa)
    cleaned_df.to_csv(output_path, sep='\t', index=False)
    print(f"? Cleaned: {filepath.name} ? {output_path.name}")

def main():
    txt_files = list(INPUT_DIR.glob("*.txt"))
    if not txt_files:
        print(f"? No .txt files found in: {INPUT_DIR}")
        return

    for file in txt_files:
        fname = file.name.lower()
        if "phylum" in fname:
            taxa = TAXA_FILTERS["phylum"]
        elif "genus" in fname:
            taxa = TAXA_FILTERS["genus"]
        elif "species" in fname:
            taxa = TAXA_FILTERS["species"]
        else:
            print(f"? Skipping unrecognized file: {file.name}")
            continue

        process_file(file, taxa)

if __name__ == "__main__":
    main()
