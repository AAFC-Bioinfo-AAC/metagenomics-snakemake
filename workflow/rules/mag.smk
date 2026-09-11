'''
    Filename: mag.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2025/09/11
    Snakemake version: 9.20.0

'''

# -------------------------------------------------------------------
# Functions for checkpoint
# -------------------------------------------------------------------
def get_filtered_samples(wildcards):
    ckpt_output = checkpoints.filter_assemblies.get().output.passed_samples
    return parse_filtered_samples(ckpt_output)


def get_filtered_list_file(wildcards):
    return checkpoints.filter_assemblies.get().output.passed_samples


# -------------------------------------------------------------------
# Checkpoint to filter assemblies based on quality metrics
# Occurs after assembly but before binning
# -------------------------------------------------------------------
checkpoint filter_assemblies:
    input:
        assemblies = expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
            sample=SAMPLES
        ),
        gate = f"{LOG_DIR}/envs/conda_gate_mag.txt"
    output:
        # CHANGED: Both files written by the checkpoint are declared outputs.
        passed_samples = f"{SAMPLE_ASSEMBLY}/passed_checkpoint_assemblies.txt",
        metrics = f"{SAMPLE_ASSEMBLY}/samples_with_contigs.metrics.tsv"
    params:
        min_len_for_stats = config.get("assembly_filter", {}).get("min_len_for_stats", 2000),
        min_total_bp = config.get("assembly_filter", {}).get("min_total_bp", 100000),
        min_contigs = config.get("assembly_filter", {}).get("min_contigs", 100),
        min_fasta_bytes = config.get("assembly_filter", {}).get("min_fasta_bytes", 1)
    run:
        import gzip
        import os

        def fasta_stats_ge_len(path, minlen):
            """Return total bp >= minlen, contigs >= minlen and all contigs."""
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
                            n_all += 1
                            if cur_len >= minlen:
                                n_ge += 1
                                total += cur_len
                        cur_len = 0
                        seen = True
                    else:
                        cur_len += len(line.strip())

                if seen:
                    n_all += 1
                    if cur_len >= minlen:
                        n_ge += 1
                        total += cur_len

            return total, n_ge, n_all

        passed = []
        rows = [(
            "sample",
            "fasta_file_bytes",
            f"total_bp_ge{params.min_len_for_stats}bp",
            f"num_contigs_ge{params.min_len_for_stats}bp",
            "num_contigs_total",
            "passed_filter",
            "filter_failure_reason"
        )]

        for infile in input.assemblies:
            sample = (
                os.path.basename(infile)
                .replace("_assembly.contigs.fa", "")
                .replace(".gz", "")
            )
            fasta_bytes = os.path.getsize(infile) if os.path.exists(infile) else 0
            total_bp, num_contigs_passing_filter, num_contigs_total = (
                fasta_stats_ge_len(infile, params.min_len_for_stats)
            )

            ok = True
            reasons = []

            if fasta_bytes < params.min_fasta_bytes:
                ok = False
                reasons.append(f"file_size<{params.min_fasta_bytes}B")
            if total_bp < params.min_total_bp:
                ok = False
                reasons.append(
                    f"total_bp_≥{params.min_len_for_stats}bp<{params.min_total_bp}"
                )
            if num_contigs_passing_filter < params.min_contigs:
                ok = False
                reasons.append(
                    f"num_contigs_≥{params.min_len_for_stats}bp<{params.min_contigs}"
                )

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

        os.makedirs(os.path.dirname(output.passed_samples), exist_ok=True)

        with open(output.passed_samples, "w", encoding="utf-8") as out:
            out.write("\n".join(passed) + ("\n" if passed else ""))

        # CHANGED: A metrics-write failure now fails the checkpoint instead of
        # being silently ignored.
        with open(output.metrics, "w", encoding="utf-8") as metrics:
            metrics.write("\t".join(rows[0]) + "\n")
            for row in rows[1:]:
                metrics.write("\t".join(row) + "\n")


# -------------------------------------------------------------------
# Pre-warm conda environments for MAG rules after the checkpoint
# -------------------------------------------------------------------
MAG_ENVS = ["bowtie2", "metabat2", "checkm2"]

localrules: prewarm_mag_env, prewarm_mag_gate


rule prewarm_mag_env:
    output:
        f"{LOG_DIR}/envs/mag/{{env}}.prewarmed"
    conda:
        "../envs/{env}.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})"
        touch {output:q}
        """


rule prewarm_mag_gate:
    input:
        expand(f"{LOG_DIR}/envs/mag/{{env}}.prewarmed", env=MAG_ENVS)
    output:
        touch(f"{LOG_DIR}/envs/conda_gate_mag.txt")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})"
        touch {output:q}
        """


# -------------------------------------------------------------------
# Assemble cleaned reads with MEGAHIT
# -------------------------------------------------------------------
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
    threads:
        config.get("megahit", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        # CHANGED: All generated-file and log directories are created explicitly.
        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {output.assembly:q})"

        tmpbase="${{TMPDIR:-/tmp}}"
        run_dir="$(mktemp -d "$tmpbase/megahit_run_XXXXXX")" || {{
            echo "Failed to create MEGAHIT run directory" >> {log:q}
            exit 1
        }}
        tmp_dir="$(mktemp -d "$tmpbase/megahit_tmp_XXXXXX")" || {{
            echo "Failed to create MEGAHIT temporary directory" >> {log:q}
            rm -rf -- "$run_dir"
            exit 1
        }}

        echo "Using MEGAHIT run directory: $run_dir" >> {log:q}
        echo "Using MEGAHIT temporary directory: $tmp_dir" >> {log:q}

        cleanup() {{
            if [[ -n "${{tmp_dir:-}}" && -d "$tmp_dir" ]]; then
                echo "Cleaning up MEGAHIT temporary directory: $tmp_dir" >> {log:q}
                rm -rf -- "$tmp_dir"
            fi
            if [[ -n "${{run_dir:-}}" && -d "$run_dir" ]]; then
                echo "Cleaning up MEGAHIT run directory: $run_dir" >> {log:q}
                rm -rf -- "$run_dir"
            fi
        }}
        trap cleanup EXIT

        out_prefix={params.out_prefix:q}

        megahit \
            -1 {input.R1:q} \
            -2 {input.R2:q} \
            -t {threads} \
            --min-contig-len {params.min_contig_length} \
            --out-dir "$run_dir" \
            --force \
            --out-prefix "$out_prefix" \
            --tmp-dir "$tmp_dir" \
            >> {log:q} 2>&1

        src_contigs="$run_dir/${{out_prefix}}.contigs.fa"
        dest={output.assembly:q}

        # CHANGED: Remove markers created by older versions of this rule.
        # The empty assembly itself is sufficient for the filtering checkpoint.
        rm -f -- "${{dest}}.EMPTY"

        if [[ ! -s "$src_contigs" ]]; then
            echo "No contigs were produced; creating empty assembly: $dest" >> {log:q}
            : > "$dest"
            exit 0
        fi

        tmp_dest="${{dest}}.tmp.$$"
        cp --preserve=mode,timestamps "$src_contigs" "$tmp_dest"
        mv -f -- "$tmp_dest" "$dest"

        echo "Placed contigs at: $dest" >> {log:q}
        """


# -------------------------------------------------------------------
# Build the Bowtie2 index for each retained assembly
# -------------------------------------------------------------------
rule index_assembly:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    output:
        # CHANGED: Force and declare the predictable large-index .bt2l format.
        temp(
            expand(
                f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.{{suffix}}",
                sample=["{sample}"],
                suffix=[
                    "1.bt2l",
                    "2.bt2l",
                    "3.bt2l",
                    "4.bt2l",
                    "rev.1.bt2l",
                    "rev.2.bt2l"
                ]
            )
        )
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_bowtie2_index.log"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config.get("index_assembly", {}).get("threads", 8)
    params:
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly"
    shadow:
        "shallow"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {params.index_base:q})"

        bowtie2-build \
            --large-index \
            --threads {threads} \
            --quiet \
            {input.assembly:q} \
            {params.index_base:q} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Map cleaned reads to each sample assembly
# -------------------------------------------------------------------
rule map_reads_to_assembly:
    input:
        index = expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.{{suffix}}",
            sample=["{sample}"],
            suffix=[
                "1.bt2l",
                "2.bt2l",
                "3.bt2l",
                "4.bt2l",
                "rev.1.bt2l",
                "rev.2.bt2l"
            ]
        ),
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        bam = f"{SAMPLE_ASSEMBLY}/{{sample}}.bam"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_bowtie2_mapping.log"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config.get("map_reads", {}).get("threads", 16)
    params:
        index_base = lambda wildcards: f"{SAMPLE_ASSEMBLY}/{wildcards.sample}_assembly",
        max_mem_per_thread = config.get("map_reads", {}).get("max_mem_per_thread", "4G")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {output.bam:q})"

        # CHANGED: Require enough threads for the two concurrently running tools.
        total_threads={threads}
        if (( total_threads < 2 )); then
            echo "map_reads_to_assembly requires at least 2 threads" >> {log:q}
            exit 1
        fi

        # Allocate about 20% of the total threads to samtools sort. The samtools
        # --threads value counts additional threads beyond its main thread.
        sort_total=$(( total_threads / 5 ))
        (( sort_total < 1 )) && sort_total=1
        t_bowtie2=$(( total_threads - sort_total ))
        t_sort_extra=$(( sort_total - 1 ))

        echo "Thread allocation: bowtie2=$t_bowtie2, samtools_total=$sort_total" >> {log:q}

        bowtie2 \
            -x {params.index_base:q} \
            -1 {input.R1:q} \
            -2 {input.R2:q} \
            --threads "$t_bowtie2" \
            2>> {log:q} \
        | samtools sort \
            -m {params.max_mem_per_thread:q} \
            --threads "$t_sort_extra" \
            -O BAM \
            -o {output.bam:q} \
            - \
            2>> {log:q}
        """


# -------------------------------------------------------------------
# Calculate per-contig depth for MetaBAT2
# -------------------------------------------------------------------
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
        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {output.depth:q})"

        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth:q} \
            {input.bam:q} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Bin contigs with MetaBAT2
# -------------------------------------------------------------------
rule metabat2_binning:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
        depth = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/{{sample}}_depth.txt"
    output:
        bins_dir = directory(f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/bins"),
        unbinned_dir = directory(f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/unbinned")
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_metabat2.log"
    conda:
        "../envs/metabat2.yaml"
    threads:
        config.get("metabat2", {}).get("threads", 8)
    params:
        min_contig_length = config.get("metabat2", {}).get("min_contig_length", 2000)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"

        tmpbase="${{TMPDIR:-/tmp}}"
        metabat2_dir="$(mktemp -d "$tmpbase/metabat2_XXXXXX")" || {{
            echo "Failed to create MetaBAT2 temporary directory" >> {log:q}
            exit 1
        }}
        echo "Using MetaBAT2 temporary directory: $metabat2_dir" >> {log:q}

        cleanup() {{
            if [[ -n "${{metabat2_dir:-}}" && -d "$metabat2_dir" ]]; then
                echo "Cleaning up MetaBAT2 temporary directory: $metabat2_dir" >> {log:q}
                rm -rf -- "$metabat2_dir"
            fi
        }}
        trap cleanup EXIT

        sample={wildcards.sample:q}
        bin_prefix="$metabat2_dir/${{sample}}.bin"

        metabat2 \
            --inFile {input.assembly:q} \
            --outFile "$bin_prefix" \
            --abdFile {input.depth:q} \
            --numThreads {threads} \
            --minContig {params.min_contig_length} \
            --unbinned \
            >> {log:q} 2>&1

        mkdir -p {output.bins_dir:q} {output.unbinned_dir:q}

        # Copy numbered genome bins.
        for file in "$bin_prefix."*[0-9].fa; do
            [[ -f "$file" ]] && \
                cp --preserve=mode,timestamps "$file" {output.bins_dir:q}/
        done

        # Copy the MetaBAT2 bin-information summary when present.
        if [[ -f "${{bin_prefix}}.BinInfo.txt" ]]; then
            cp --preserve=mode,timestamps \
                "${{bin_prefix}}.BinInfo.txt" \
                {output.bins_dir:q}/
        fi

        # Copy special outputs that are not genome bins.
        for file in "${{bin_prefix}}."{{tooShort,lowDepth,unbinned}}.fa; do
            [[ -f "$file" ]] && \
                cp --preserve=mode,timestamps "$file" {output.unbinned_dir:q}/
        done

        bin_count=$(find {output.bins_dir:q} -maxdepth 1 -type f -name '*.fa' | wc -l)
        echo "MetaBAT2 produced $bin_count genome bin(s)" >> {log:q}
        echo "MetaBAT2 outputs copied to {output.bins_dir} and {output.unbinned_dir}" >> {log:q}
        """


# -------------------------------------------------------------------
# Estimate MAG completeness and contamination with CheckM2
# -------------------------------------------------------------------
rule checkm2_bins:
    input:
        bins_dir = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/bins",
        checkm2_db = f"{CHECKM2_DB}"
    output:
        # CHANGED: Do not declare both a directory and a file within it.
        checkm2_summary = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/checkm2/quality_report.tsv",
        # CHANGED: Record whether CheckM2 ran or no bins were available.
        checkm2_status = f"{SAMPLE_ASSEMBLY}/metabat2/{{sample}}/checkm2/status.tsv"
    log:
        f"{LOG_DIR}/individual_assemblies/{{sample}}_checkm2.log"
    conda:
        "../envs/checkm2.yaml"
    threads:
        config.get("checkm2", {}).get("threads", 4)
    params:
        memory_usage = config.get("checkm2", {}).get("memory_usage", "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {output.checkm2_summary:q})"

        # CHANGED: An empty MetaBAT2 result is a valid biological outcome. Do
        # not call CheckM2 with an empty input directory.
        if ! find {input.bins_dir:q} -maxdepth 1 -type f -name '*.fa' -print -quit \
            | grep -q .; then
            echo "No genome bins were produced; CheckM2 was not run" >> {log:q}
            : > {output.checkm2_summary:q}
            printf 'sample\tstatus\tdetails\n%s\tNO_BINS\tMetaBAT2 produced no genome bins\n' \
                {wildcards.sample:q} > {output.checkm2_status:q}
            exit 0
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        checkm2_run_dir="$(mktemp -d "$tmpbase/checkm2_XXXXXX")" || {{
            echo "Failed to create CheckM2 temporary directory" >> {log:q}
            exit 1
        }}

        cleanup() {{
            if [[ -n "${{checkm2_run_dir:-}}" && -d "$checkm2_run_dir" ]]; then
                echo "Cleaning up CheckM2 temporary directory: $checkm2_run_dir" >> {log:q}
                rm -rf -- "$checkm2_run_dir"
            fi
        }}
        trap cleanup EXIT

        checkm2 predict \
            --threads {threads} \
            {params.memory_usage} \
            -x fa \
            --database_path {input.checkm2_db:q} \
            --input {input.bins_dir:q} \
            --output-directory "$checkm2_run_dir" \
            >> {log:q} 2>&1

        if [[ ! -s "$checkm2_run_dir/quality_report.tsv" ]]; then
            echo "CheckM2 did not produce a non-empty quality_report.tsv" >> {log:q}
            exit 1
        fi

        checkm2_dest="$(dirname {output.checkm2_summary:q})"
        cp -a "$checkm2_run_dir"/. "$checkm2_dest"/

        printf 'sample\tstatus\tdetails\n%s\tCHECKM2_COMPLETED\tquality_report.tsv created\n' \
            {wildcards.sample:q} > {output.checkm2_status:q}
        """
