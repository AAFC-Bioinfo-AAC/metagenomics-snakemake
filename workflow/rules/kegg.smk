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
rule make_kegg_diamond_db:
    input:
        prokaryotes_fasta = f"{KEGG_FASTA}/prokaryotes.pep.gz"
    output:
        done_kegg_diamond = f"{LOG_DIR}/prokaryotes_db_done.txt",
        diamond_db = f"{KEGG_DIAMOND}/prokaryotes.pep.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads: config.get("make_kegg_diamond_db", {}).get("threads", 8)
    shell:
        r"""
        DIAMOND_DB="{KEGG_DIAMOND}/prokaryotes.pep.dmnd"
        KEGG_PROK_FASTA="{input.prokaryotes_fasta}"
        THREADS={threads}
        DONE="{output.done_kegg_diamond}"

        mkdir -p $(dirname {output.done_kegg_diamond})

        if [[ ! -f "$DIAMOND_DB" ]]; then
            echo "DIAMOND database not found. Creating from $KEGG_PROK_FASTA"
            
            if [[ ! -f "$KEGG_PROK_FASTA" ]]; then
                echo "Error: FASTA file not found at $KEGG_PROK_FASTA"
                exit 1
            fi

            zcat "$KEGG_PROK_FASTA" | diamond makedb --in - -d "${KEGG_DIAMOND}/prokaryotes.pep" --threads "$THREADS"
            echo "DIAMOND database created successfully."
        else
            echo "DIAMOND database already exists. Skipping creation."
        fi

        date > "$DONE"
        """
rule kegg_diamond:
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz",
        diamond_db = f"{KEGG_DIAMOND}/prokaryotes.pep.dmnd"
    output:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8",
        temp_fastq = temp(f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_tmp.fastq")
    log:
        f"{LOG_DIR}/kegg/{{sample}}_kegg_diamond.log"
    conda:
        "../envs/diamond.yaml"
    params:
        pigz_threads=2,
        sensitivity=config.get("kegg_diamond", {}).get("sensitivity", ""), #default is not setting sensitivity, and default is designed for finding hits of >60% identity and short read alignment. Its sensitivity is between --fast and --mid-sensitive
        max_target_num=config.get("kegg_diamond", {}).get("max-target-seqs", 1),
        out_file_format=config.get("kegg_diamond", {}).get("outfmt", "6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore")
    threads: config.get("kegg_diamond", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.diamond})"

        pigz_threads={params.pigz_threads}
        diamond_threads=$(( {threads} - pigz_threads ))
        [[ "$diamond_threads" -lt 1 ]] && diamond_threads=1

        pigz -dc -p "$pigz_threads" {input.merged} > {output.temp_fastq} 2>> {log} || {{ echo "Failed to decompress {input.merged}" >> {log}; exit 1; }}

        tmpbase="${{TMPDIR:-/tmp}}"
        diamond_tmp="$(mktemp -d "$tmpbase/diamond_{wildcards.sample}_XXXXXX")" || {{ echo "Failed to create DIAMOND temp dir" >> {log}; exit 1; }}
        echo "Using DIAMOND tmp dir: ${{diamond_tmp}}" >> {log}

        cleanup() {{
            if [[ -n "${{diamond_tmp:-}}" && -d "${{diamond_tmp}}" ]]; then
                echo "Cleaning up DIAMOND tmp dir: ${{diamond_tmp}}" >> {log}
                rm -rf -- "${{diamond_tmp}}"
            fi
        }}
        trap cleanup EXIT

        diamond blastx \
            -d {input.diamond_db} \
            -q {output.temp_fastq} \
            -o {output.diamond} \
            --tmpdir "$diamond_tmp" \
            {params.sensitivity} \
            --max-target-seqs {params.max_target_num} \
            --outfmt {params.out_file_format} \
            --threads "$diamond_threads" \
            2>> {log}
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
    params:
        minpath_script = "workflow/scripts/MinPath/MinPath.py"

    shell:
        r"""
        echo "MinPath version is 1.6" >> {log}
        if [[ ! -s {input.ko_list_fixed} ]]; then
            echo "[`date '+%Y-%m-%d %H:%M:%S'`] KO list is empty for {wildcards.sample}. Skipping MinPath." | tee -a {log}
            touch {output.minpath_out}
        else
            python {params.minpath_script} -ko {input.ko_list_fixed} -report {output.minpath_out} &>> {log}
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
rule kegg_category_mapping:
    input:
        minpath_table = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_aggregated_minpath.tsv",
        BRITE_hierarchy = f"{KEGG_BRITE_HIERARCHY}/ko00001.keg"
    output:
        f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_pathway_abundance_with_category.tsv"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_kegg_category_mapping.log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/kegg_category_summary.py"
rule kegg_category_sampleID:
    input:
        f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_pathway_abundance_with_category.tsv"
    output:
        f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_pathway_abundance_with_category_sampleID.tsv"
    shell:
        """
        awk -v id={wildcards.sample} 'BEGIN {{OFS="\\t"}} NR==1 {{print "sampleID", $0; next}} {{print id, $0}}' {input} > {output}
        """
rule combine_kegg_category_tables:
    input:
        expand(f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_pathway_abundance_with_category_sampleID.tsv", sample=SAMPLES)
    output:
        f"{KEGG_OUTPUT_DIR}/combined_ko_pathway_abundance_with_category.tsv"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -q -n +2 {input} >> {output}
        """ 
rule filter_combined_kegg_table:
    input:
        combined = f"{KEGG_OUTPUT_DIR}/combined_ko_pathway_abundance_with_category.tsv",
        exclude_list = f"{KEGG_CUSTOM_LIST}/KEGG_BRITE_pathway_exclusion_file.txt"
    output:
        f"{KEGG_OUTPUT_DIR}/combined_ko_pathway_abundance_with_category_filtered.tsv"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/filter_combined_kegg_table.py"  

