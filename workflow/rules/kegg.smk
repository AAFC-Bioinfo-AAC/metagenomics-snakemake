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
    log:
        f"{LOG_DIR}/{kegg}/kegg_merge.log"
    shell:
        r"""

        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.merged})"

        "cat {input.R1} {input.R2} > {output.merged}"
        echo "Merged {input.R1} and {input.R2} into {output.merged}" >> {log}
        """

rule: kegg_diamond
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz",
        diamond_db = f"{KEGG_DIR}/prokaryotes.pep.dmnd"

    output:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_kegg_diamond.log"
    conda:
        "../envs/diamond.yaml"
    params:
        pigz_threads=4
    threads: config.get("kegg_diamond", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.diamond})"

        # per job temp folder
        tmpbase="${TMPDIR:-/tmp}"
        jobtmp="$tmpbase/diamond_{wildcards.sample}_$RANDOM"
        mkdir -p "$jobtmp" || { echo "Failed to create job tmpdir $jobtmp" >> "{log}"; exit 1; }
        echo "Using job tmpdir base: $jobtmp" >> "{log}"
        piz_temp="$jobtmp/piz_tmp"
        diamond_temp="$jobtmp/diamond_tmp"
        mkdir -p "$piz_temp" "$diamond_temp"

        cleanup() {
            if [[ -n "$jobtmp" && -d "$jobtmp" && "$jobtmp" == "$tmpbase"/diamond_* ]]; then
                rm -rf "$jobtmp"
            fi
        }
        trap cleanup EXIT

        # Set pigz threads from Snakemake param
        pigz_threads={params.pigz_threads}
        diamond_threads=$(( {threads} - pigz_threads ))
        # Ensure no thread count is 0
        [ "$diamond_threads" -lt 1 ] && diamond_threads=1

        # Decompress first, then run DIAMOND
        pigz -dc -p "$pigz_threads" {input.merged} > "$piz_temp/{wildcards.sample}.fastq" 2>> {log} || { echo "Failed to decompress {input.merged}" >> "{log}"; exit 1; }

        diamond blastx \
        -d {input.diamond_db} \
        -q "$piz_temp/{wildcards.sample}.fastq" \
        -o {output.diamond} \
        --tmpdir "$diamond_temp" \
        --sensitive \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --threads "$diamond_threads"
        """
rule count_reads:
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz"
    output:
        read_count = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_read_count.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.read_count})"

        "zcat {input.merged} | wc -l | awk '{{print int($1/4)}}' > {output.read_count}"
        """

rule gene_ko_abundance:
    input:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8",
        read_count = f"{KEGG_OUTPUT_DIR}/read_counts/{{sample}}_read_count.txt",
        ko_genes_list = KEGG_DIR #/gpfs/fs7/aafc/common/reference/kegg/genes/ko/ko_genes.list
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
        minpath = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_minpath_output.txt"
    log:
        f"{LOG_DIR}/kegg/{{sample}}_minpath.log"
    conda:
        "../envs/minpath.yaml"
    script:
        "../scripts/minpath.py"

rule: kegg_minpath
    input:
        merged = f"{MERGED_R1_R2}/{{sample}}_merged.fastq.gz",
        ko_genes_list = KEGG_DIR,
        ko_pathway_list =KEGG_DIR,
    
    output:
        diamond = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_diamond_output.m8",
        ko = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_gene_ko_abundance.tsv",
        ko_list_raw = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_list_raw.txt",
        ko_list_fixed = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_list_fixed.txt",
        minpath = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_minpath_output.txt",
        pathway = f"{KEGG_OUTPUT_DIR}/{{sample}}/{{sample}}_ko_pathway_abundance.tsv"
    log:
        f"{LOG_DIR}/{kegg}/{{sample}}_kegg_minpath.log"
    conda:
        "../envs/diamond.yaml"
    scripts:
        "../scripts/MinPath.py" #need to move in to scripts
    params:
    shell:
        r"""
        # read count 
        read_count=$(zcat {input.merged} | wc -l)
        read_count=$((read_count/4))
    
        # Gene-to-KO Mapping and Abundance


        
        
        """

