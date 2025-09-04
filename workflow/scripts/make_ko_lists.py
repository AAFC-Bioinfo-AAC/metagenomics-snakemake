'''
    Filename: make_ko_list.py
    Author: Katherine James-Gzyl
    Date created: 2025/09/03
    Snakemake version: 9.9.0
    python version: 3.10

    Script generates KO lists from the abundance table
      1. Extract KO IDs from the abundance table
      2. Save raw KO list
      3. Fix KO list formatting
      4. Save fixed KO list
'''
import pandas as pd
from pathlib import Path

# Redirect logging
if hasattr(snakemake, "log") and snakemake.log:
    sys.stdout = open(snakemake.log[0], "w")
    sys.stderr = sys.stdout

# Input/output from Snakemake
abundance_file = snakemake.input["abundance"]
ko_list_raw = snakemake.output["ko_list_raw"]
ko_list_fixed = snakemake.output["ko_list_fixed"]

# Read KO abundance table
df = pd.read_csv(abundance_file, sep="\t", dtype=str)

# Extract and filter KO column
kos = df["KO"]
kos = kos[~kos.isin(["No_KO", "KO"])]           # Exclude 'No_KO' and header
kos = kos.str.replace("ko:", "", regex=False)   # Remove 'ko:' prefix

unique_kos = sorted([k for k in kos if pd.notnull(k) and k.strip()])

# Write one-per-line KO list (raw)
Path(ko_list_raw).write_text("\n".join(unique_kos) + "\n")

if not unique_kos:
    print("KO list is empty. Skipping MinPath and aggregation.")

# Write two-column KO list for MinPath (fixed)
with open(ko_list_fixed, "w") as out:
    for i, ko in enumerate(unique_kos, 1):
        out.write(f"dummy{i}\t{ko}\n")