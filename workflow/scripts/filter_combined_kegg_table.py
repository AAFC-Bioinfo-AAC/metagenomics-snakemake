'''
    Filename: aggregate_minpath_pathways.py
    Author: Devin Holman, Katherine James-Gzyl, AI
    Date created: 2025/09/26
    Snakemake version: 9.9.0
    python version: 3.10

    Script is filtering out non-microbial pathways from the combined KEGG pathway abundance table.
    The exclusion list can be edited as needed and is located in reasources/KEGG_BRITE_pathway_exclusion_file.txt
'''

import pandas as pd

# Snakemake gives you paths directly
combined_file = str(snakemake.input.combined)
exclude_file = str(snakemake.input.exclude_list)
output_file = str(snakemake.output)   # no [0] needed

# Load files
df = pd.read_csv(combined_file, sep="\t")
exclude = pd.read_csv(exclude_file, sep="\t")

# Ensure IDs are strings to avoid mismatch (e.g., "00073" vs 73)
exclude_ids = exclude['Pathway_ID'].astype(str).tolist()

# Filter out unwanted pathways
filtered = df[~df['Pathway'].astype(str).isin(exclude_ids)]

# Save
filtered.to_csv(output_file, sep="\t", index=False)

