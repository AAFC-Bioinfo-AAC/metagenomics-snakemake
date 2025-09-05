'''
    Filename: kgg.smk
    Author: Katherine James-Gzyl
    Date created: 2025/08/25
    Snakemake version: 9.9.0
'''
rule merge_read_pairs:
    input: 
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz"
    shell:
        r"""
        mkdir -p "$(dirname {output.merged})"

        zcat {input.R1} {input.R2} | gzip > {output.merged}
        """
rule kegg_diamond:
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz",
        diamond_db = f"{KEGG_DIAMOND}/prokaryotes.pep.dmnd"
    output:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_kegg_diamond.log"
    conda:
        "../envs/diamond.yaml"
    params:
        pigz_threads=4,
        sensitivity=config.get("kegg_diamond", {}).get("sensitivity", "--sensitive"),
        max_target_num=config.get("kegg_diamond", {}).get("max-target-seqs", 1),
        out_file_format=config.get("kegg_diamond", {}).get("outfmt", "6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore")

    threads: config.get("kegg_diamond", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.diamond})"

        # per job temp folder
        tmpbase="${{TMPDIR:-/tmp}}"
        jobtmp="$(mktemp -d "$tmpbase/diamond_{wildcards.sample}_XXXXXX")" || {{ echo "Failed to create job tmpdir $jobtmp" >> "{log}"; exit 1; }}
        echo "Using job tmpdir base: $jobtmp" >> "{log}"
        piz_temp="$jobtmp/piz_tmp"
        diamond_temp="$jobtmp/diamond_tmp"
        mkdir -p "$piz_temp" "$diamond_temp"

        cleanup() {{
             # Extra sanity check: never allow jobtmp == tmpbase
            if [[ "$jobtmp" == "$tmpbase" ]]; then
                echo "Sanity check failed: jobtmp is exactly tmpbase—refusing to delete." >> "{log}"
                return
            fi
            
            if [[ -n "$jobtmp" && -d "$jobtmp" && "$jobtmp" == "$tmpbase"/diamond_* ]]; then
                rm -rf "$jobtmp"
            fi
        }}
        trap cleanup EXIT

        # Set pigz threads from Snakemake param
        pigz_threads={params.pigz_threads}
        diamond_threads=$(( {threads} - pigz_threads ))
        # Ensure no thread count is 0
        [ "$diamond_threads" -lt 1 ] && diamond_threads=1

        # Decompress first, then run DIAMOND
        pigz -dc -p "$pigz_threads" {input.merged} > "$piz_temp/{wildcards.sample}.fastq" 2>> {log} || {{ echo "Failed to decompress {input.merged}" >> "{log}"; exit 1; }}

        diamond blastx \
        -d {input.diamond_db} \
        -q "$piz_temp/{wildcards.sample}.fastq" \
        -o {output.diamond} \
        --tmpdir "$diamond_temp" \
        {params.sensitivity} \
        {params.max_target_num} \
        {params.out_file_format} \
        --threads "$diamond_threads"
        """
rule count_reads:
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz"
    output:
        read_count = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_read_count.txt"
    conda:
        "../envs/python3.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.read_count})"

        pigz -dc {input.merged} | wc -l | awk '{{print int($1/4)}}' > {output.read_count}
        """
rule gene_ko_abundance:
    input:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8",
        read_count = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_read_count.txt",
        ko_genes_list = f"{KEGG_KO}/ko_genes.list"
    output:
       abundance = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_gene_ko_abundance.tsv", 
    log:
        f"{LOG_DIR}/kegg/{{sample}}_gene_ko_abundance.log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/gene_ko_abundance.py"
rule make_ko_lists:
    input:
        abundance = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_gene_ko_abundance.tsv"
    output:
        ko_list_raw = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_list_raw.txt",
        ko_list_fixed = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_list_fixed.txt"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_make_ko_lists.log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/make_ko_lists.py"
rule minpath:
    input:
        ko_list_fixed = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_list_fixed.txt"
    output:
        minpath_out= f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_minpath_output.txt"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_minpath.log"
    conda:
        "../envs/minpath.yaml"
    shell:
        r"""
        echo "MinPath version is 1.6" >> {log}
        if [[ ! -s {input.ko_list_fixed} ]]; then
            echo "[`date '+%Y-%m-%d %H:%M:%S'`] KO list is empty for {wildcards.sample}. Skipping MinPath." | tee -a {log}
            touch {output.minpath_out}
        else
            python {MINPATH_SCRIPT} -ko {input.ko_list_fixed} -report {output.minpath_out} &>> {log}
        fi
        """
rule aggregate_minpath_pathways:
    input:
        minpath_out= f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_minpath_output.txt",
        abundance = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_gene_ko_abundance.tsv",
        ko_pathway_list = f"{KEGG_KO}/ko_pathway.list"
    output:
        minpath_table= f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_aggregated_minpath.tsv"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_aggregate_minpath.log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/aggregate_minpath_pathways.py"
