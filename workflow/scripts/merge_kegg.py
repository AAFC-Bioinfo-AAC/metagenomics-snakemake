'''
    Filename: merge_kegg.py
    Author: Devin Holman, Katherine James-Gzyl and AI
    Date created: 2025/12/11
    Snakemake version: 9.9.0
    python version: 3.10

    Script combines the ko_pathway_abundance_with_category.tsv, aggregated_minpath.tsv, and gene_ko_abundance.tsv into the following combined files:
        1) Pathways_categorized_CPM.tsv 
            column names are: Pathway,Pathway_Name,Top_Category,Sub_Category,sample1,sample2,...
            rows are unique Pathways
        2) Pathways_no_categorization_CPM.tsv
            column names are: Pathway,sample1,sample2,...
            rows are unique Pathways
        3) KEGG_gene_hits_raw.tsv
            column names are: Gene,KO,sample1,sample2,...
            rows are copies per million reads values for each Gene-KO pair 
        4) KO_CPM.tsv
            column names are: KO,sample1,sample2,...
            rows are unique KOs with summed Copies_Per_Million_Reads values across Genes per sample
        5) Read_counts_per_sample.tsv
            column names are: Sample,Read_Count
            rows are samples with their corresponding read counts extracted from the gene_ko_abundance.tsv files
'''
import csv
import sys
from collections import OrderedDict, defaultdict

# Print all stdout & stderr to snakemake log file
if snakemake.log:
    log_file = snakemake.log[0] if isinstance(snakemake.log, (list, tuple)) else snakemake.log
    log_handle = open(log_file, "w")
    sys.stdout = log_handle
    sys.stderr = log_handle

print("Starting merge_kegg_results...")
print(f"Number of input files for pathways with category: {len(snakemake.input.pathways_cat)}")
print(f"Number of input files for pathways without category: {len(snakemake.input.pathways_nocat)}")
print(f"Number of gene_ko files: {len(snakemake.input.gene_ko)}")


def merge_files(id_cols, info_cols, metric_col, files_by_sample, out_file):
    """
    Merge per-sample TSVs into a wide matrix.

    Parameters
    ----------
    id_cols : list[str]
        Column(s) that uniquely identify a row (e.g. ["Pathway"] or ["Gene", "KO"])
    info_cols : list[str]
        Extra annotation columns to preserve (e.g. ["Top_Category"])
    metric_col : str
        Column containing the numerical metric to merge (e.g. "total_cpm")
    files_by_sample : dict[str, str]
        Mapping {sample_name: filepath}
    out_file : str
        Output TSV filename
    """
    files_by_sample = dict(files_by_sample)
    if not files_by_sample:
        print(f"No files provided for output {out_file}")
        return

    all_samples = sorted(files_by_sample.keys())
    print(f"\n=== Output: {out_file} ===")
    print(f"Found {len(all_samples)} samples.")

    all_keys = OrderedDict()      # key -> info tuple
    data = defaultdict(dict)      # key -> {sample: metric_value (string)}

    for sample, filepath in files_by_sample.items():
        print(f"  Processing: {filepath} as sample {sample}")

        with open(filepath, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Empty or malformed file: {filepath}")

            # Allow for whitespace in headers
            header_map = {h.strip(): h for h in reader.fieldnames}
            required = id_cols + info_cols + [metric_col]

            for req in required:
                if req not in header_map:
                    raise ValueError(
                        f"Missing required column '{req}' in {filepath}. "
                        f"Found: {reader.fieldnames}"
                    )

            for row in reader:
                key = tuple(row[header_map[col]] for col in id_cols)

                if key not in all_keys:
                    all_keys[key] = tuple(
                        row[header_map[col]] for col in info_cols
                    )

                # keep metric as string (no rounding / conversion)
                data[key][sample] = row[header_map[metric_col]]

    # Write wide matrix (no Combined column)
    with open(out_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        header = id_cols + info_cols + all_samples
        writer.writerow(header)

        for key, info in all_keys.items():
            row = list(key) + list(info)
            for sample in all_samples:
                row.append(data[key].get(sample, "0"))
            writer.writerow(row)

    print(f"  Wrote {out_file}")


def merge_gene_ko_files_and_readcounts(files_by_sample, ko_cpm_out, readcounts_out):
    """
    From per-sample *_gene_ko_abundance.tsv files, produce:

      1) ko_cpm.tsv:
         - Rows: KO
         - Columns: samples
         - Values: sum of Copies_Per_Million_Reads per KO per sample

      2) read_counts_per_sample.tsv:
         - Columns: Sample, Read_Count
    """
    files_by_sample = dict(files_by_sample)
    if not files_by_sample:
        print(f"No files provided for KO merge into {ko_cpm_out}")
        return

    all_samples = sorted(files_by_sample.keys())
    print(f"\n=== Output: {ko_cpm_out} & {readcounts_out} ===")
    print(f"Found {len(all_samples)} samples.")

    all_kos = OrderedDict()
    data = defaultdict(lambda: defaultdict(float))  # KO -> sample -> CPM
    sample_to_readcount = {}

    for sample, filepath in files_by_sample.items():
        print(f"  Processing: {filepath} as sample {sample}")
        got_readcount = False

        with open(filepath, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Empty or malformed file: {filepath}")

            for row in reader:
                KO = row["KO"].strip()

                # Grab per-sample read count from the first row we see it
                if not got_readcount and "Read_Count" in row:
                    sample_to_readcount[sample] = row["Read_Count"]
                    got_readcount = True

                # Skip non-KO entries
                if KO.lower() == "no_ko":
                    continue

                cpm_str = row["Copies_Per_Million_Reads"]
                if cpm_str in ("", None):
                    continue

                try:
                    cpm = float(cpm_str)
                except ValueError:
                    # Skip weird non-numeric values
                    continue

                all_kos[KO] = None
                data[KO][sample] += cpm

    # --- ko_cpm.tsv ---
    with open(ko_cpm_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["KO"] + all_samples)

        for KO in all_kos:
            row = [KO]
            for sample in all_samples:
                value = data[KO].get(sample, 0.0)
                row.append(f"{value:.6f}")
            writer.writerow(row)

    print(f"  Wrote {ko_cpm_out}")

    # --- read_counts_per_sample.tsv ---
    with open(readcounts_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Sample", "Read_Count"])
        for sample in all_samples:
            rc = sample_to_readcount.get(sample, "NA")
            writer.writerow([sample, rc])

    print(f"  Wrote {readcounts_out}")


# ================== Entry point when called via Snakemake ==================

# List of sample names (same order as in expand(..., sample=SAMPLE_NAMES))
samples = list(snakemake.params.samples)

# Build {sample: file} mappings from the input lists
pathways_cat_by_sample = dict(zip(samples, snakemake.input.pathways_cat))
pathways_nocat_by_sample = dict(zip(samples, snakemake.input.pathways_nocat))
gene_ko_by_sample = dict(zip(samples, snakemake.input.gene_ko))

# Pathways with categories
merge_files(
    id_cols=["Pathway"],
    info_cols=["Pathway_Name", "Top_Category", "Sub_Category"],
    metric_col="total_cpm",
    files_by_sample=pathways_cat_by_sample,
    out_file=str(snakemake.output.pathways_categorized),
)

# Pathways without categories
merge_files(
    id_cols=["Pathway"],
    info_cols=[],
    metric_col="total_cpm",
    files_by_sample=pathways_nocat_by_sample,
    out_file=str(snakemake.output.pathways_nocat),
)

# Gene–KO table (raw, no combined/aggregation)
merge_files(
    id_cols=["Gene", "KO"],
    info_cols=[],
    metric_col="Copies_Per_Million_Reads",
    files_by_sample=gene_ko_by_sample,
    out_file=str(snakemake.output.gene_hits_raw),
)

# KO CPM + read counts
merge_gene_ko_files_and_readcounts(
    files_by_sample=gene_ko_by_sample,
    ko_cpm_out=str(snakemake.output.ko_cpm),
    readcounts_out=str(snakemake.output.read_counts_per_sample),
)

print("\nAll merges completed.")
if snakemake.log:
    log_handle.close()