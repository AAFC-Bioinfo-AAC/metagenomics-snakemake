# extract_bracken_columns.py (for Snakemake script:)

import pandas as pd
import re

input_tables = {
    "species": snakemake.input.species_table,
    "genus": snakemake.input.genus_table,
    "phylum": snakemake.input.phylum_table
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

for level in ["species", "genus", "phylum"]:
    df = pd.read_csv(input_tables[level], sep='\t')

    raw_suffix = f"_bracken.{level}.report.txt_num"
    rel_suffix = f"_bracken.{level}.report.txt_frac"

    # Raw columns
    num_columns = list(df.columns[:3]) + [c for c in df.columns if c.endswith(raw_suffix)]
    abundance_num = df[num_columns]
    abundance_num.columns = [
        re.sub(re.escape(raw_suffix) + '$', '', c)
        for c in abundance_num.columns
    ]
    abundance_num.to_csv(output_raw[level], index=False, quoting=0)

    # Relative columns
    frac_columns = list(df.columns[:3]) + [c for c in df.columns if c.endswith(rel_suffix)]
    abundance_frac = df[frac_columns]
    abundance_frac.columns = [
        re.sub(re.escape(rel_suffix) + '$', '', c)
        for c in abundance_frac.columns
    ]
    abundance_frac.to_csv(output_rel[level], index=False, quoting=0)