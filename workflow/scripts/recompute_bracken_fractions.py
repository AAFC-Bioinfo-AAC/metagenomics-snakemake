#!/usr/bin/env python3
"""
Efficiently recalculates *_frac_adjusted columns in species, genus, and phylum tables,
by dividing *_num by each sample's total (Bacteria + Archaea) reads from the domain table.
Adds new *_frac_adjusted columns in one operation.
"""

import pandas as pd
import numpy as np
import sys
import re

if hasattr(snakemake, "log") and snakemake.log:
    sys.stdout = open(snakemake.log[0], "w")
    sys.stderr = sys.stdout

domain_path = snakemake.input["domain"]
inputs = {
    "species": snakemake.input["species"],
    "genus": snakemake.input["genus"],
    "phylum": snakemake.input["phylum"],
}
outputs = {
    "species": snakemake.output["species_adjusted"],
    "genus": snakemake.output["genus_adjusted"],
    "phylum": snakemake.output["phylum_adjusted"],
}

domain_df = pd.read_csv(domain_path, sep='\t')

def extract_sample_id(col):
    # Extracts the full sample ID before _bracken
    # e.g., test_SUB2008RP9D70 or H5CK2004RP4D70 from *_bracken.*.txt_num
    m = re.match(r"^(.+?)_bracken", col)
    return m.group(1) if m else None

# Build a {sample_id: domain_col} map for *_num in domain table
domain_num_cols = [c for c in domain_df.columns if c.endswith('_num')]
domain_col_for_sample = {extract_sample_id(c):c for c in domain_num_cols if extract_sample_id(c)}

print("DEBUG: Detected sample IDs from domain table:")
for sample_id in domain_col_for_sample.keys():
    print(f"  - {sample_id}")

# Compute prokaryote (Bacteria+Archaea) total for each sample_id
prokaryote_totals = {}
for sample_id, col in domain_col_for_sample.items():
    total = float(domain_df.loc[domain_df['name'] == 'Bacteria', col].values.sum()) + \
            float(domain_df.loc[domain_df['name'] == 'Archaea', col].values.sum())
    prokaryote_totals[sample_id] = total
    print(f"DEBUG: {sample_id} -> prokaryote total = {total}")

def adjust_table(infile, outfile, prokaryote_totals, level):
    df = pd.read_csv(infile, sep='\t')
    sample_num_cols = [c for c in df.columns if c.endswith('_num')]
    new_cols = {}
    adjusted_colnames = []  # To track order for sum QC

    for col in sample_num_cols:
        sample_id = extract_sample_id(col)
        total = prokaryote_totals.get(sample_id, None)
        new_col = col.replace("_num", "_frac_adjusted")
        adjusted_colnames.append(new_col)
        if total is not None and total > 0:
            new_cols[new_col] = df[col].astype(float) / total
        else:
            # Create a Series of zeros with the same index as df
            new_cols[new_col] = pd.Series(0.0, index=df.index)
            print(f"WARNING: No matching prokaryote total for {col} (sample_id={sample_id}), all zeros.")

    # Add all new columns at once for performance
    adjusted_df = pd.DataFrame(new_cols, index=df.index)
    df = pd.concat([df, adjusted_df], axis=1)

    # QC: check sum per fractional column
    sums = df[adjusted_colnames].sum()
    print(f"Saved recalculated abundances for {level} -> {outfile}")
    print("-- ADJUSTED frac sum per sample (should all be ~= 1.0) --")
    print(sums.to_string())
    if not np.allclose(sums, 1.0, atol=1e-4):
        print("WARNING: Some adjusted fractions do not sum to 1.0 within tolerance (1e-4)!")
    else:
        print("All adjusted columns sum to 1.0 within tolerance.")

    df.to_csv(outfile, sep='\t', index=False)

for level in ["species", "genus", "phylum"]:
    adjust_table(inputs[level], outputs[level], prokaryote_totals, level)