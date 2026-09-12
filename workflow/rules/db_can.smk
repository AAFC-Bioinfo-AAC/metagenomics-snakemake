'''
    Filename: db_can.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2026/09/11
    Snakemake version: 9.20.0

'''

import os


# -------------------------------------------------------------------
# Functions for the non-empty-assembly checkpoint
# -------------------------------------------------------------------

def get_filtered_samples_dbcan(wildcards):
    """Return samples whose assembly file exists and is non-empty."""
    checkpoint_output = checkpoints.nonempty_assemblies.get().output[0]
    return parse_filtered_samples(checkpoint_output)

def bwa_mapping_threads():
    try:
        value = int(config.get("bwa_mem_mapping", {}).get("threads", 12))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "bwa_mem_mapping.threads must be an integer of at least 2."
        ) from exc

    if value < 2:
        raise ValueError(
            "bwa_mem_mapping.threads must be at least 2 so BWA and SAMtools "
            "each receive a thread."
        )
    return value


# -------------------------------------------------------------------
# Checkpoint to retain samples with a non-empty assembly
# -------------------------------------------------------------------

checkpoint nonempty_assemblies:
    input:
        assemblies = expand(
            f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
            sample=SAMPLE_NAMES
        ),
        gate = f"{LOG_DIR}/envs/conda_gate_dbcan.txt"
    output:
        f"{SAMPLE_ASSEMBLY}/nonempty_assemblies.txt"
    run:
        valid_samples = []
        suffix = "_assembly.contigs.fa"

        for assembly in input.assemblies:
            if os.path.isfile(assembly) and os.path.getsize(assembly) > 0:
                filename = os.path.basename(assembly)
                sample = (
                    filename[:-len(suffix)]
                    if filename.endswith(suffix)
                    else os.path.splitext(filename)[0]
                )
                valid_samples.append(sample)

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w", encoding="utf-8") as handle:
            handle.write(
                "\n".join(valid_samples) + ("\n" if valid_samples else "")
            )


# -------------------------------------------------------------------
# Pre-warm conda environments for the dbCAN workflow after the checkpoint
# -------------------------------------------------------------------

DBCAN_ENVS = ["pyrodigal", "dbcan", "bwa"]

localrules: prewarm_dbcan_env, prewarm_dbcan_gate


rule prewarm_dbcan_env:
    output:
        f"{LOG_DIR}/envs/dbcan/{{env}}.prewarmed"
    conda:
        "../envs/{env}.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})"
        touch {output:q}
        """


rule prewarm_dbcan_gate:
    input:
        expand(
            f"{LOG_DIR}/envs/dbcan/{{env}}.prewarmed",
            env=DBCAN_ENVS
        )
    output:
        touch(f"{LOG_DIR}/envs/conda_gate_dbcan.txt")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})"
        touch {output:q}
        """


# -------------------------------------------------------------------
# Gene prediction
# -------------------------------------------------------------------

rule pyrodigal:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    output:
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        nuc_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}.cds"
    log:
        f"{LOG_DIR}/dbcan/prodigal/{{sample}}.log"
    threads:
        config.get("pyrodigal", {}).get("threads", 8)
    conda:
        "../envs/pyrodigal.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.gff:q})"
        : > {log:q}

        pyrodigal \
            -p meta \
            -i {input.assembly:q} \
            -d {output.nuc_file:q} \
            -a {output.faa:q} \
            -f gff \
            -o {output.gff:q} \
            -j {threads} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Optional run_dbCAN analyses
#
# These rules deliberately remain separate so users can request the most
# comprehensive analysis they require. easy_CGC repeats CAZyme annotation and
# easy_substrate repeats both CAZyme annotation and CGC identification.
# -------------------------------------------------------------------

rule cazyme_annotation:
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_cazyme/overview.tsv"
    params:
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_cazyme"
    log:
        f"{LOG_DIR}/dbcan/cazyme_annotation/{{sample}}.log"
    threads:
        config.get("cazyme_annotation", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan CAZyme_annotation \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


rule cgc_calling:
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul/overview.tsv",
        cgc = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul/cgc.gff",
        cgc_table = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul/cgc_standard_out.tsv"
    params:
        # The existing _pul name is retained for compatibility with the
        # current workflow/Snakefile targets.
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul"
    log:
        f"{LOG_DIR}/dbcan/cgc_calling/{{sample}}.log"
    threads:
        config.get("cgc_calling", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan easy_CGC \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --input_gff {input.gff:q} \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


rule substrate_prediction:
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/overview.tsv",
        cgc_table = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/cgc_standard_out.tsv",
        substrate = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/substrate_prediction.tsv"
    params:
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan"
    log:
        f"{LOG_DIR}/dbcan/substrate_prediction/{{sample}}.log"
    threads:
        config.get("substrate_prediction", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan easy_substrate \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --input_gff {input.gff:q} \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Read mapping and gene-depth calculation
# -------------------------------------------------------------------

rule bwa_mem_mapping:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        bam = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam",
        bai = temp(f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam.bai"),
        bwa_idx = temp(
            multiext(
                f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa"
            )
        )
    log:
        f"{LOG_DIR}/dbcan/bwa_mem_mapping/{{sample}}.log"
    threads:
        bwa_mapping_threads()
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.bam:q})"
        : > {log:q}

        sort_threads=$(( {threads} / 4 ))
        if (( sort_threads < 1 )); then
            sort_threads=1
        fi
        bwa_threads=$(( {threads} - sort_threads ))

        echo "[$(date)] Total threads: {threads}; BWA: $bwa_threads; SAMtools: $sort_threads" \
            >> {log:q}

        bwa index {input.assembly:q} >> {log:q} 2>&1

        bwa mem \
            -t "$bwa_threads" \
            {input.assembly:q} \
            {input.R1:q} \
            {input.R2:q} \
            2>> {log:q} \
        | samtools sort \
            -@ "$sort_threads" \
            -o {output.bam:q} \
            - \
            >> {log:q} 2>&1

        samtools index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """


rule dbcan_depth:
    input:
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        bam = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam.bai"
    output:
        depth_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    log:
        f"{LOG_DIR}/dbcan/dbcan_depth/{{sample}}.log"
    params:
        overlap_base_ratio = config.get("dbcan_depth", {}).get(
            "overlap_base_ratio", 0.2
        ),
        mapping_quality = config.get("dbcan_depth", {}).get(
            "mapping_quality", 30
        ),
        identity = config.get("dbcan_depth", {}).get("identity", 0.98)
    threads:
        config.get("dbcan_depth", {}).get("threads", 16)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.depth_file:q})"
        : > {log:q}

        dbcan_utils cal_coverage \
            -g {input.gff:q} \
            -i {input.bam:q} \
            -o {output.depth_file:q} \
            -t {threads} \
            --overlap_base_ratio {params.overlap_base_ratio} \
            --mapping_quality {params.mapping_quality} \
            --identity {params.identity} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Abundance calculation
# -------------------------------------------------------------------

rule get_abundances_rpm:
    input:
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/overview.tsv",
        cgc_table = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/cgc_standard_out.tsv",
        substrate = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/substrate_prediction.tsv",
        depth_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    output:
        cazy_fam_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/fam_abund.out",
        cazy_subfam_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/subfam_abund.out",
        EC_number = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/EC_abund.out",
        substrate_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/fam_substrate_abund.out",
        cgc_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_abund.out",
        substrate_ho = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_substrate_PUL_homology.out",
        cgc_substrate_voting = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_substrate_majority_voting.out"
    params:
        dbcan_dir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan",
        ab_dir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund",
        marker_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/did_not_run_get_abundances_rpm.txt"
    log:
        f"{LOG_DIR}/dbcan/get_abundances_rpm/{{sample}}.log"
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.ab_dir:q}
        : > {log:q}
        rm -f -- {params.marker_file:q}

        # Return success only when a tabular file contains at least one
        # non-empty row after its header.
        has_data_rows() {{
            local input_file="$1"
            [[ -s "$input_file" ]] && \
                awk 'NR > 1 && NF {{ found = 1; exit }} END {{ exit !found }}' \
                    "$input_file"
        }}

        require_output() {{
            local output_file="$1"
            if [[ ! -e "$output_file" ]]; then
                echo "ERROR: Expected output was not created: $output_file" \
                    >> {log:q}
                exit 1
            fi
        }}

        record_skip() {{
            printf '%s\n' "$1" >> {params.marker_file:q}
        }}

        (
            cd {params.ab_dir:q}

            if has_data_rows {input.overview:q}; then
                dbcan_utils fam_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cazy_fam_ab:q}
                require_output {output.cazy_subfam_ab:q}
                require_output {output.EC_number:q}

                dbcan_utils fam_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ab:q}
            else
                : > {output.cazy_fam_ab:q}
                : > {output.cazy_subfam_ab:q}
                : > {output.EC_number:q}
                : > {output.substrate_ab:q}
                record_skip \
                    "Family-level abundance calculations were skipped because overview.tsv contained no annotation rows."
            fi

            if has_data_rows {input.cgc_table:q}; then
                dbcan_utils CGC_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cgc_ab:q}
            else
                : > {output.cgc_ab:q}
                record_skip \
                    "CGC abundance calculation was skipped because cgc_standard_out.tsv contained no CGC rows."
            fi

            if has_data_rows {input.substrate:q}; then
                dbcan_utils CGC_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ho:q}
                require_output {output.cgc_substrate_voting:q}
            else
                : > {output.substrate_ho:q}
                : > {output.cgc_substrate_voting:q}
                record_skip \
                    "CGC substrate abundance calculations were skipped because substrate_prediction.tsv contained no prediction rows."
            fi
        )

        """
