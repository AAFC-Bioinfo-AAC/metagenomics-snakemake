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

        src_contigs="$run_dir/{params.out_prefix}.contigs.fa"
        dest="{output.assembly}"

        if [[ ! -s "$src_contigs" ]]; then
            echo "No contigs produced; creating empty placeholder: $dest" >> {log}
            : > "$dest"
            touch "${{dest}}.EMPTY"
            echo "Marker file created: ${{dest}}.EMPTY" >> {log}
            exit 0
        fi

        tmp_dest="${{dest}}.tmp.$$"
        cp --preserve=mode,timestamps "$src_contigs" "$tmp_dest"
        mv -f "$tmp_dest" "$dest"

        echo "Placed contigs to: $dest" >> {log}
        """
checkpoint filter_assemblies:
    input:
        expand(f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa", sample=SAMPLES)
    output:
        f"{SAMPLE_ASSEMBLY}/passed_checkpoint_assemblies.txt"
    params:
        min_len_for_stats = config.get("assembly_filter", {}).get("min_len_for_stats", 2000),
        min_total_bp      = config.get("assembly_filter", {}).get("min_total_bp", 100000),
        min_contigs       = config.get("assembly_filter", {}).get("min_contigs", 100),
        min_fasta_bytes   = config.get("assembly_filter", {}).get("min_fasta_bytes", 1),
        metrics_tsv       = f"{SAMPLE_ASSEMBLY}/samples_with_contigs.metrics.tsv"
    run:
        import os, gzip

        def fasta_stats_ge_len(path, minlen):
            """Return (total_bp_ge_min, n_contigs_ge_min, n_contigs_all).
               Robust to empty/missing; supports .gz."""
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                return 0, 0, 0

            opener = gzip.open if path.endswith(".gz") else open
            total = 0
            n_ge = 0
            n_all = 0
            cur_len = 0
            seen = False
            with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seen:
                            # finalize previous contig
                            n_all += 1
                            if cur_len >= minlen:
                                n_ge += 1
                                total += cur_len
                        cur_len = 0
                        seen = True
                    else:
                        cur_len += len(line.strip())
                # finalize last contig
                if seen:
                    n_all += 1
                    if cur_len >= minlen:
                        n_ge += 1
                        total += cur_len
            return total, n_ge, n_all

        passed = []
        # More descriptive column names
        rows = [(
            "sample",
            "fasta_file_bytes",
            f"total_bp_ge{params.min_len_for_stats}bp",
            f"num_contigs_ge{params.min_len_for_stats}bp",
            "num_contigs_total",
            "passed_filter",
            "filter_failure_reason"
        )]

        for infile in input:
            sample = os.path.basename(infile).replace("_assembly.contigs.fa", "").replace(".gz","")
            fasta_bytes = os.path.getsize(infile) if os.path.exists(infile) else 0
            total_bp, num_contigs_passing_filter, num_contigs_total = fasta_stats_ge_len(infile, params.min_len_for_stats)

            ok = True
            reasons = []
            if fasta_bytes < params.min_fasta_bytes:
                ok = False
                reasons.append(f"file_size<{params.min_fasta_bytes}B")
            if total_bp < params.min_total_bp:
                ok = False
                reasons.append(f"total_bp_≥{params.min_len_for_stats}bp<{params.min_total_bp}")
            if num_contigs_passing_filter < params.min_contigs:
                ok = False
                reasons.append(f"num_contigs_≥{params.min_len_for_stats}bp<{params.min_contigs}")

            if ok:
                passed.append(sample)
            
            rows.append((
                sample,
                str(fasta_bytes),
                str(total_bp),
                str(num_contigs_passing_filter),
                str(num_contigs_total),
                "PASS" if ok else "FAIL",
                ",".join(reasons) if reasons else "OK"
            ))

        with open(output[0], "w") as out:
            out.write("\n".join(passed) + ("\n" if passed else ""))

        # Optional sidecar metrics
        try:
            os.makedirs(os.path.dirname(params.metrics_tsv), exist_ok=True)
            with open(params.metrics_tsv, "w") as m:
                m.write("\t".join(rows[0]) + "\n")
                for r in rows[1:]:
                    m.write("\t".join(r) + "\n")
        except Exception:
            pass

rule index_assembly:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    output:
        expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.{{suffix}}",
            sample=["{sample}"],
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_bowtie2_index.log"
    conda:
        "../envs/bowtie2.yaml"
    threads: config.get("index_assembly", {}).get("threads", 8)
    params:
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly"
    shell:
        r"""
        set -euo pipefail
        bowtie2-build --threads {threads} --quiet {input.assembly} {params.index_base} >> {log} 2>&1
        """
rule map_reads_to_assembly:
    input:
        index = expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.{{suffix}}",
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
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly",
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
rule sample_depth_file:
    input:
        bam = f"{SAMPLE_ASSEMBLY}/{{sample}}.bam"
    output:
        depth = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/{{sample}}_depth.txt"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_depth.log"
    conda:
        "../envs/metabat2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.depth})"
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam} 2>> {log}
        """
rule metabat2_binning:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
        depth    = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/{{sample}}_depth.txt"
    output:
        bins_dir     = directory(f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/bins"),
        unbinned_dir = directory(f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/unbinned")
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_metabat2.log"
    conda:
        "../envs/metabat2.yaml"
    threads: config.get("metabat2", {}).get("threads", 8)
    params:
        min_contig_length = config.get("metabat2", {}).get("min_contig_length", 2000)
    shell:
        r"""
        set -euo pipefail
        
        tmpbase="${{TMPDIR:-/tmp}}"
        metabat2_dir="$(mktemp -d "$tmpbase/metabat2_{wildcards.sample}_XXXXXX")" || {{
            echo "Failed to create metabat2 dir" >> {log}; exit 1;
        }}
        echo "Using metabat2 dir: ${{metabat2_dir}}" >> {log}

        cleanup() {{
            if [[ -n "${{metabat2_dir:-}}" && -d "${{metabat2_dir}}" ]]; then
                echo "Cleaning up metabat2 dir: ${{metabat2_dir}}" >> {log}
                rm -rf -- "${{metabat2_dir}}"
            fi
        }}
        trap cleanup EXIT

        # Run MetaBAT2 in scratch
        metabat2 --inFile {input.assembly} \
            --outFile "${{metabat2_dir}}/{wildcards.sample}.bin" \
            --abdFile {input.depth} \
            --numThreads {threads} \
            --minContig {params.min_contig_length} \
            --unbinned 2>> {log}

        # Ensure destination dirs exist
        mkdir -p "{output.bins_dir}" "{output.unbinned_dir}"

        # Copy binned FASTAs and their BinInfo.txt files, preserving timestamps
        for f in "${{metabat2_dir}}/{wildcards.sample}.bin."*[0-9].fa; do
            [[ -f "$f" ]] && cp --preserve=mode,timestamps "$f" "{output.bins_dir}/"
        done
        # Copy the single BinInfo summary file
        if [[ -f "${{metabat2_dir}}/{wildcards.sample}.bin.BinInfo.txt" ]]; then
            cp --preserve=mode,timestamps "${{metabat2_dir}}/{wildcards.sample}.bin.BinInfo.txt" "{output.bins_dir}/"
        fi

         # Copy special bins (tooShort, lowDepth, unbinned)
        for f in "${{metabat2_dir}}/{wildcards.sample}.bin."{{tooShort,lowDepth,unbinned}}.fa; do
            [[ -f "$f" ]] && cp --preserve=mode,timestamps "$f" "{output.unbinned_dir}/"
        done

        echo "MetaBAT2 outputs copied to {output.bins_dir} and {output.unbinned_dir}" >> {log}
        """
rule checkm2_bins:
    input:
        bins_dir = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/bins",
        checkm2_db = f"{CHECKM2_DB}"
    output:
        checkm2_dir = directory(f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/checkm2/"),
        checkm2_summary = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/checkm2/quality_report.tsv"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_checkm2.log"
    conda:
        "../envs/checkm2.yaml"
    threads: config.get("checkm2", {}).get("threads", 4)
    params:
        memory_usage=config.get("checkm2", {}).get("memory_usage", "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.checkm2_dir})"
        checkm2 predict \
            --threads {threads} \
            {params.memory_usage} \
            -x fa \
            --database_path {input.checkm2_db} \
            --input {input.bins_dir} \
            --output-directory {output.checkm2_dir} 2>> {log}
        """






        