# extract_bracken_columns.py (for Snakemake script:)

import pandas as pd
import re

#!/usr/bin/env python3
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

    # Accept both variants with/without '.cleaned' in the filename
    raw_suffixes = [
        f"_bracken.{level}.report.cleaned.txt_num",
        f"_bracken.{level}.report.txt_num"
    ]
    rel_suffixes = [
        f"_bracken.{level}.report.cleaned.txt_frac",
        f"_bracken.{level}.report.txt_frac"
    ]

    # Find columns that end with any suffix (raw)
    num_columns = list(df.columns[:3]) + [
        c for c in df.columns if any(c.endswith(suf) for suf in raw_suffixes)
    ]
    # Remove suffix (with or without .cleaned) for output columns
    abundance_num = df[num_columns]
    abundance_num.columns = [
        re.sub(f"(_bracken\.{level}\.report(\.cleaned)?\.txt_num)$", "", c)
        for c in abundance_num.columns
    ]
    abundance_num.to_csv(output_raw[level], index=False, quoting=0)

    # Find columns ending with frac-suffix
    frac_columns = list(df.columns[:3]) + [
        c for c in df.columns if any(c.endswith(suf) for suf in rel_suffixes)
    ]
    abundance_frac = df[frac_columns]
    abundance_frac.columns = [
        re.sub(f"(_bracken\.{level}\.report(\.cleaned)?\.txt_frac)$", "", c)
        for c in abundance_frac.columns
    ]
    abundance_frac.to_csv(output_rel[level], index=False, quoting=0)