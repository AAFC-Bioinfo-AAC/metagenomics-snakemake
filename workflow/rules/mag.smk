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
          > {log} 2>&1

        # MEGAHIT writes {params.out_prefix}.contigs.fa inside out-dir
        src_contigs="$run_dir/{params.out_prefix}.contigs.fa"
        if [[ ! -s "$src_contigs" ]]; then
          echo "Expected contigs not found: $src_contigs" >> {log}
          exit 2
        fi

        # Atomic place into final destination (handles cross-filesystem safely)
        dest="{output.assembly}"
        tmp_dest="${{dest}}.tmp.$$"
        cp --preserve=mode,timestamps "$src_contigs" "$tmp_dest"
        mv -f "$tmp_dest" "$dest"

        echo "Placed contigs to: $dest" >> {log}
        """