# workflow/scripts/helpers.py

def bowtie2_index_files(sample, SAMPLE_ASSEMBLY):
    """Return the 6 Bowtie2 index files for a given sample."""
    return [
        f"{SAMPLE_ASSEMBLY}/{sample}_assembly.bt2.{i}.bt2" for i in range(1, 5)
    ] + [
        f"{SAMPLE_ASSEMBLY}/{sample}_assembly.bt2.rev.{i}.bt2" for i in (1, 2)
    ]

def parse_filtered_samples(filtered_list_file):
    """Read filtered samples file and return as list of sample names."""
    with open(filtered_list_file) as f:
        return [line.strip() for line in f if line.strip()]
