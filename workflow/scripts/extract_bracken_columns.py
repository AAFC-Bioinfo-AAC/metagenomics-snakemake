#!/usr/bin/env python3
"""
This script extracts the relevant columns from the merged Bracken output files after the
relative abundance had been recalculated by dividing the reads for each taxon by the total
bacteria + archaea reads for that sample.

At each Taxonomic level there are three output files:
1. Raw counts
2. Bracken fractions (default, not adjusted)
3. Recalculated abundance based on Bracken Domain read totals for bacteria + archaea (adjusted)
"""
import pandas as pd
import re

def extract_sample(col):
    # Only strip from sample columns (not metadata columns)
    match = re.match(r"^([A-Za-z0-9\-]+)_bracken.*", col)
    return match.group(1) if match else col

input_tables = {
    "species": snakemake.input.species_adjusted,
    "genus": snakemake.input.genus_adjusted,
    "phylum": snakemake.input.phylum_adjusted
}
output_raw = {
    "species": snakemake.output.species_raw,
    "genus": snakemake.output.genus_raw,
    "phylum": snakemake.output.phylum_raw
}
output_rel = {
    "species": snakemake.output.species_rel,
    "genus": snakemake.output.genus_rel,
    "phylum": snakemake.output.phylum_rel
}
output_rel_adjusted = {
    "species": snakemake.output.species_rel_recalc,
    "genus": snakemake.output.genus_rel_recalc,
    "phylum": snakemake.output.phylum_rel_recalc
}

for level in ["species", "genus", "phylum"]:
    df = pd.read_csv(input_tables[level], sep='\t')
    meta_cols = list(df.columns[:3])

    # --- RAW ABUNDANCE ---
    num_cols = meta_cols + [c for c in df.columns if c.endswith('_num')]
    out_num = df[num_cols].copy()
    out_num.columns = [extract_sample(c) if c not in meta_cols else c for c in out_num.columns]
    out_num.to_csv(output_raw[level], index=False)

    # --- BRACKEN FRACTIONS (DEFAULT, not adjusted) ---
    frac_cols = meta_cols + [c for c in df.columns if c.endswith('_frac') and not c.endswith('_frac_adjusted')]
    out_frac = df[frac_cols].copy()
    out_frac.columns = [extract_sample(c) if c not in meta_cols else c for c in out_frac.columns]
    out_frac.to_csv(output_rel[level], index=False)

    # --- ADJUSTED FRACTIONS ---
    adj_cols = meta_cols + [c for c in df.columns if c.endswith('_frac_adjusted')]
    out_adj = df[adj_cols].copy()
    out_adj.columns = [extract_sample(c) if c not in meta_cols else c for c in out_adj.columns]
    out_adj.to_csv(output_rel_adjusted[level], index=False)

    print(f"Extracted tables for {level}: raw={output_raw[level]}, rel={output_rel[level]}, rel_adjusted={output_rel_adjusted[level]}")