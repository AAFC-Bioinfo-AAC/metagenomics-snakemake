'''
    Filename: kegg_category_summary.py
    Author: Devin Holman and Katherine James-Gzyl
    Date created: 2025/09/10
    Snakemake version: 9.9.0
    python version: 3.10

    The BRITE hierarchy file from the KEGG database is used to summarize the MinPath-confirmed pathways into higher-level categories
    for each sample using the sample_aggregated_minpath.tsv file.
      1. Parse BRITE hierarchy (ko00001.keg)
      2. Annotate aggregated MinPath output (sample_aggregated_minpath.tsv)
      3. Confirm KEGG IDs are 5-digit strings
      4. Save the results as a tab delimited table
'''
import pandas as pd
import re
import sys
import traceback

# Snakemake output logging (redirect stdout & stderr)
if "snakemake" in globals():
    if hasattr(snakemake, "log") and snakemake.log:
        sys.stdout = open(snakemake.log[0], "w")
        sys.stderr = sys.stdout

try:
    # Input/output from Snakemake
    minpath_table = snakemake.input["minpath_table"]
    brite_file = snakemake.input["BRITE_hierarchy"]
    output_file = snakemake.output[0]

    print(f"[INFO] minpath_table: {minpath_table}")
    print(f"[INFO] brite_file: {brite_file}")
    print(f"[INFO] output_file: {output_file}")

    # Parse BRITE hierarchy file to build a mapping of pathway IDs to categories
    pathway_categories = {}
    top_category, sub_category = None, None

    try:
        with open(brite_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("A"):
                    top_category = line[1:].strip()
                elif line.startswith("B"):
                    sub_category = line[1:].strip()
                elif line.startswith("C"):
                    # Example line:
                    # C 00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
                    match = re.match(r"C\s+(\d{5})\s+(.*?)\s*\[", line)
                    if match:
                        pid, pname = match.groups()
                        pathway_categories[pid] = {
                            "name": pname.strip(),
                            "top_category": top_category,
                            "sub_category": sub_category
                        }
        print(f"[INFO] Parsed {len(pathway_categories)} pathways from BRITE hierarchy.")
    except Exception as e:
        print(f"[ERROR] Failed parsing BRITE file '{brite_file}': {e}")
        traceback.print_exc()
        sys.exit(1)

    # Helper lookup functions
    def get_name(pid): 
        return pathway_categories.get(pid, {}).get("name", "Unknown")

    def get_top_category(pid): 
        return pathway_categories.get(pid, {}).get("top_category", "Unknown")

    def get_sub_category(pid): 
        return pathway_categories.get(pid, {}).get("sub_category", "Unknown")

    # Annotate MinPath output
    try:
        df = pd.read_csv(minpath_table, sep="\t")
        print(f"[INFO] Loaded {minpath_table}: {df.shape[0]} rows, columns: {', '.join(df.columns)}")
    except Exception as e:
        print(f"[ERROR] Could not read input MinPath table '{minpath_table}': {e}")
        traceback.print_exc()
        sys.exit(1)

    if "Pathway" not in df.columns:
        print(f"[ERROR] Input file {minpath_table} does not contain a 'Pathway' column. Columns present: {df.columns.to_list()}")
        sys.exit(1)

    # Ensure KEGG IDs are 5-digit strings
    df["Pathway"] = df["Pathway"].astype(str).str.zfill(5)

    # Add KEGG annotations
    df["Pathway_Name"] = df["Pathway"].apply(get_name)
    df["Top_Category"] = df["Pathway"].apply(get_top_category)
    df["Sub_Category"] = df["Pathway"].apply(get_sub_category)

    unmatched = df[df["Pathway_Name"] == "Unknown"]["Pathway"].unique()
    if len(unmatched) > 0:
        print(f"[WARNING] {len(unmatched)} pathways in MinPath not found in BRITE hierarchy: {', '.join(unmatched)}")

   
    # Save annotated file
    try:
        df.to_csv(output_file, sep="\t", index=False)
        print(f"[INFO] Wrote annotated pathways to {output_file}")
    except Exception as e:
        print(f"[ERROR] Could not write output file '{output_file}': {e}")
        traceback.print_exc()
        sys.exit(1)

except Exception as e:
    print("[FATAL ERROR] Unhandled exception occurred!")
    traceback.print_exc()
    sys.exit(1)