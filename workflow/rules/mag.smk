'''
    Filename: mag.smk
    Author: Katherine James-Gzyl
    Date created: 2025/09/23
    Snakemake version: 9.9.0
'''
rule megahit_assembly:
    input: 
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_megahit.log"
    conda:
        "../envs/megahit.yaml"
    params:
        min_contig_length = config.get("megahit", {}).get("min_contig_length", 1000), 
        out_prefix = config.get("megahit", {}).get("out_prefix", "final")
    threads: config.get("megahit", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.assembly})"

        tmpbase="${{TMPDIR:-/tmp}}"
        run_dir="$(mktemp -d "$tmpbase/magahit_run_{wildcards.sample}_XXXXXX")" || {{ echo "Failed to create MEGAHIT run dir" >> {log}; exit 1; }}
        tmp_dir="$(mktemp -d "$tmpbase/magahit_tmp_{wildcards.sample}_XXXXXX")" || {{ echo "Failed to create MEGAHIT temp dir" >> {log}; exit 1; }}
        echo "Using MEGAHIT run dir: ${{run_dir}}" >> {log}
        echo "Using MEGAHIT tmp dir: ${{tmp_dir}}" >> {log}

        cleanup() {{
            if [[ -n "${{tmp_dir:-}}" && -d "${{tmp_dir}}" ]]; then
                echo "Cleaning up MEGAHIT tmp dir: ${{tmp_dir}}" >> {log}
                rm -rf -- "${{tmp_dir}}"
            fi

            if [[ -n "${{run_dir:-}}" && -d "${{run_dir}}" ]]; then
                echo "Cleaning up MEGAHIT run dir: ${{run_dir}}" >> {log}
                rm -rf -- "${{run_dir}}"
            fi
        }}
        trap cleanup EXIT

        megahit \
          -1 {input.R1} \
          -2 {input.R2} \
          -t {threads} \
          --min-contig-len {params.min_contig_length} \
          --out-dir "$run_dir" --force \
          --out-prefix {params.out_prefix} \
          --tmp-dir "$tmp_dir" \
          >> {log} 2>&1

        # MEGAHIT writes {params.out_prefix}.contigs.fa inside out-dir
        src_contigs="$run_dir/{params.out_prefix}.contigs.fa"
        if [[ ! -s "$src_contigs" ]]; then
          echo "Expected contigs not found: $src_contigs; creating empty placeholder file" >> {log}
          : > "{output.assembly}"
        
          exit 2
        fi

        # Atomic place into final destination (handles cross-filesystem safely)
        dest="{output.assembly}"
        tmp_dest="${{dest}}.tmp.$$"
        cp --preserve=mode,timestamps "$src_contigs" "$tmp_dest"
        mv -f "$tmp_dest" "$dest"

        echo "Placed contigs to: $dest" >> {log}
        """
checkpoint filter_nonempty_assemblies:
    input:
        expand(f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa", sample=SAMPLE_NAMES)
    output:
         f"{SAMPLE_ASSEMBLY}/samples_with_contigs.txt"
    run:
        import os
        valid_samples = []
        for infile in input:
            if os.path.isfile(infile) and os.path.getsize(infile) > 0:
                sample = os.path.basename(infile).replace("_assembly.contigs.fa", "")
                valid_samples.append(sample)
        with open(output[0], "w") as out:
            out.write("\n".join(valid_samples) + "\n")

rule index_assembly:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    output:
        expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.bt2.{{suffix}}",
            sample=["{sample}"],
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_bowtie2_index.log"
    conda:
        "../envs/bowtie2.yaml"
    threads: config.get("index_assembly", {}).get("threads", 8)
    params:
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly.bt2"
    shell:
        r"""
        set -euo pipefail
        bowtie2-build --threads {threads} --quiet {input.assembly} {params.index_base} >> {log} 2>&1
        """
rule map_reads_to_assembly:
    input:
        index = expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.bt2.{{suffix}}",
            sample=["{sample}"],
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        ),
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        bam = f"{SAMPLE_ASSEMBLY}/{{sample}}.bam"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_bowtie2_mapping.log"
    conda:
        "../envs/bowtie2.yaml"
    threads: config.get("map_reads", {}).get("threads", 16)
    params:
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly.bt2",
        max_mem_per_thread = config.get("map_reads", {}).get("max_mem_per_thread", "4G")
    shell:
        r"""
        set -euo pipefail

        # 80/20 split (rounded)
        t_bowtie2=$(( ({threads} * 80) / 100 ))
        t_sort=$(( {threads} - t_bowtie2 ))

        [ $t_bowtie2 -lt 1 ] && t_bowtie2=1
        [ $t_sort -lt 1 ] && t_sort=1

        echo "thread splitting: bowtie2=$t_bowtie2, samtools=$t_sort" >> {log}

        bowtie2 -x {params.index_base} -1 {input.R1} -2 {input.R2} \
            --threads $t_bowtie2 2>> {log} \
        | samtools sort -m {params.max_mem_per_thread} --threads $t_sort \
            -O BAM -o {output.bam} - 2>> {log}
        """






        