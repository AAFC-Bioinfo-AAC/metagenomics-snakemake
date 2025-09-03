'''
    Filename: gene_ko_abundance.py
    Author: Devin Holman and Katherine James-Gzyl
    Date created: 2025/09/03
    Snakemake version: 9.9.0
    python version: 3.8 #double check

    Script processes DIAMOND output for each sample to generate a KO-annotated abundance table
      1. Map genes to KEGG Orthology IDS
      2. Count gene occurrences and calculate lengths
      3. Normalize counts to RPK (reads per kilobase) and CPM (counts per million)
      4. Save the results as a tab delimited table
'''

import pandas as pd
import sys

# Redirect logging
if hasattr(snakemake, "log") and snakemake.log:
    sys.stdout = open(snakemake.log[0], "w")
    sys.stderr = sys.stdout

# Input/output from Snakemake
diamond_output = snakemake.input["diamond"]
ko_genes_list = snakemake.input["ko_genes"]
output_file = snakemake.output["abundance"]
with open(snakemake.input["read_count"]) as f:
    read_count = int(f.read().strip())

# Map each gene to its KO (KEGG Orthology) ID
ko_mapping = {}
with open(ko_genes_list, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            ko, gene = parts
            ko_mapping[gene.strip()] = ko.strip()

# Parse the DIAMOND output to build two dictionaries:
#   - gene_counts: number of times each gene appears
#   - gene_lengths: nucleotide length of each gene
gene_counts = {}
gene_lengths = {}
with open(diamond_output, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        gene = parts[1]
        gene_length_aa = int(parts[2])
        gene_length_nt = (gene_length_aa * 3) + 3
        gene_lengths.setdefault(gene, gene_length_nt)
        gene_counts[gene] = gene_counts.get(gene, 0) + 1

# Loop takes each gene, looks up its KO ID, length, and count, then calculates two normalizations
#   - RPK (Reads per Kilobase): adjusts counts by gene length.
#   - CPM (Counts per Million): adjusts counts by sequencing depth.

output_data = []
for gene, count in gene_counts.items():
    ko = ko_mapping.get(gene, "No_KO")
    length_kb = gene_lengths[gene] / 1000
    rpk = count / length_kb if length_kb > 0 else 0
    cpm = (count / read_count * 1e6) if read_count > 0 else 0
    output_data.append((gene, ko, count, rpk, cpm, read_count))

df = pd.DataFrame(output_data, columns=["Gene", "KO", "Abundance", "RPK", "Copies_Per_Million_Reads", "Read_Count"])
df.to_csv(output_file, sep='\t', index=False)